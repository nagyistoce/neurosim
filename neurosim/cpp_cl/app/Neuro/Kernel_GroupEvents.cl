/*
  TODO:
    - bug. Fails if GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER is not 16 for WG sizes more than 1 and
      if GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT is set.
    - detect events with identical time and fuse them into a single event (need to sync with update
      phase on that)
    - instead of using atomics on output histogram it could use same method of computing histogram
      as it does for local one.
    - tik-tok buffer size can be reduced by 1/3 since 2nd key space is not used. Key and value
      buffers could be separate ones.
    - enable radix greater than 4.
    - implement scaling with GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG. May be useful for 
      narrowing WG size for variants.
    - for multi-WF per WG try to explore best-performance barreir placement.
*/

  /*Bariers and fences. Enable in case of errors.*/
#define GROUP_EVENTS_LOCAL_BARRIER_WF                     //barrier(CLK_LOCAL_MEM_FENCE);
#define GROUP_EVENTS_LOCAL_MEM_FENCE_WF                   //mem_fence(CLK_LOCAL_MEM_FENCE);

  /*Local memory accessors*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  #define HISTOGRAM(i)                                    *(cache + \
                                                          GROUP_EVENTS_WF_CACHE_SIZE_WORDS*\
                                                          GROUP_EVENTS_WG_SIZE_WF + i)
  #define OFFSET(i)                                       *((int __local*)(cache + i))
  /*Alternative form: __global unsigned short *x = (unsigned short __global*)(test+idx);*/
#endif

#if (GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT)
#if GROUP_EVENTS_ELEMENTS_PER_WI != 4
  #error GROUP_EVENTS_ELEMENTS_PER_WI has to be 4
#endif

uint localPrefixSum
(
          uint  pData, 
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
          uint  wi_id, 
          uint  wf_id, 
#else
          uint  wi_id, 
#endif 
          uint* totalSum, 
  __local uint  cache[]
){
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  uint wfOffset  = wf_id*GROUP_EVENTS_WF_CACHE_SIZE_WORDS;
  cache[wfOffset + wi_id] = 0;
  wfOffset += (GROUP_EVENTS_WF_SIZE_WI - GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET);
  cache[wfOffset + wi_id] = pData;
  uint idx = wfOffset + 2*wi_id + 1;
#else
  cache[wi_id] = 0;
  cache[wi_id + GROUP_EVENTS_WF_SIZE_WI - 
    GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] = pData;
  uint idx = 2*wi_id + (GROUP_EVENTS_WF_SIZE_WI + 1) - 
    GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET;
#endif

	GROUP_EVENTS_LOCAL_BARRIER_WF
  
  /*Prefix sum*/
#if (GROUP_EVENTS_OPTIMIZATION_2LEVEL_REDUCE)
  uint u0, u1, u2;
  u0 = cache[idx-3];
  u1 = cache[idx-2];
  u2 = cache[idx-1];
  atomic_add( &cache[idx], u0+u1+u2 );
  GROUP_EVENTS_LOCAL_MEM_FENCE_WF

  u0 = cache[idx-12];
  u1 = cache[idx-8];
  u2 = cache[idx-4];
  atomic_add( &cache[idx], u0+u1+u2 );
  GROUP_EVENTS_LOCAL_MEM_FENCE_WF

  u0 = cache[idx-48];
  u1 = cache[idx-32];
  u2 = cache[idx-16];
  atomic_add( &cache[idx], u0+u1+u2 );
  GROUP_EVENTS_LOCAL_MEM_FENCE_WF
  
  cache[idx-1] += cache[idx-2];
  GROUP_EVENTS_LOCAL_MEM_FENCE_WF
#else
  cache[idx] += cache[idx-1];
  mem_fence(CLK_LOCAL_MEM_FENCE);
  cache[idx] += cache[idx-2];			
  mem_fence(CLK_LOCAL_MEM_FENCE);
  cache[idx] += cache[idx-4];
  mem_fence(CLK_LOCAL_MEM_FENCE);
  cache[idx] += cache[idx-8];
  mem_fence(CLK_LOCAL_MEM_FENCE);
  cache[idx] += cache[idx-16];
  mem_fence(CLK_LOCAL_MEM_FENCE);
  cache[idx] += cache[idx-32];
  mem_fence(CLK_LOCAL_MEM_FENCE);
  cache[idx-1] += cache[idx-2];
  mem_fence(CLK_LOCAL_MEM_FENCE);
#endif

	GROUP_EVENTS_LOCAL_BARRIER_WF
  
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
	wfOffset--;
  *totalSum = cache[wfOffset + GROUP_EVENTS_WF_SIZE_WI];
	uint addValue = cache[wfOffset + wi_id];
#else
	*totalSum = cache[GROUP_EVENTS_WF_SIZE_WI*2 - 1 - 
    GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
	uint addValue = cache[wi_id + GROUP_EVENTS_WF_SIZE_WI - 1 - 
    GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
#endif
  
	return addValue;
}



/*
  2 scan, 2 exchange.
  If GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI > 256 may result in overflow of allocated
  4 partitions per 32 bit word. Best if performed by 64-wide WFs independetly (no sync).
*/
void sort4Bits1
(
#if (GROUP_EVENTS_DEBUG_ENABLE)
  __global uint  *gm_debug_device,
#endif
#if GROUP_EVENTS_VALUES_MODE
  uint dataValues[GROUP_EVENTS_ELEMENTS_PER_WI], 
#endif
  uint dataKeys[GROUP_EVENTS_ELEMENTS_PER_WI], 
  int startBit, 
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  int wi_id, 
  int wf_id, 
#else
  int wi_id, 
#endif
  __local uint* cache
){
  /*Compute static variables*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  uint wfOffset  = wf_id*GROUP_EVENTS_WF_CACHE_SIZE_WORDS;
  uint dstAddr = wfOffset + GROUP_EVENTS_ELEMENTS_PER_WI*wi_id;
#else
  uint dstAddr = GROUP_EVENTS_ELEMENTS_PER_WI*wi_id;
#endif

  /*ibit = 0, 2, ... GROUP_EVENTS_HISTOGRAM_BIN_BITS-1*/
	for(uint ibit = startBit; ibit < GROUP_EVENTS_HISTOGRAM_BIN_BITS+startBit; ibit += 2)
	{
    /*WI: Right-shift and store my current 2-bit value starting from startBit inside a uint4*/
    uint4 b = (uint4)((dataKeys[0]>>(ibit)) & 0x3, 
                      (dataKeys[1]>>(ibit)) & 0x3, 
                      (dataKeys[2]>>(ibit)) & 0x3, 
                      (dataKeys[3]>>(ibit)) & 0x3);

		/*WI: partition a word: each key type (0,1,2,3) is allocated a chunk of bits. Store keys within
      bit chunks.*/

    /*For 64 WIs and 4 keys per WI the max possible sum is 256. Thus, allocate the following bit 
      chunks: 0:[7-0], 1:[15-8], 2:[23-16], 3:[31-24]. Pack using mult 8*(0,1,2,3) = (0x1, 0x100, 
      0x10000,  0x1000000).*/
    uint key4;
		uint sKeyPacked[GROUP_EVENTS_ELEMENTS_PER_WI];
    sKeyPacked[0] = 1<<(8*b.x);
    sKeyPacked[1] = 1<<(8*b.y);
    sKeyPacked[2] = 1<<(8*b.z);
    sKeyPacked[3] = 1<<(8*b.w);
    /*Sum the keys within their bit chunks and scan them.*/
    key4 = sKeyPacked[0] + sKeyPacked[1] + sKeyPacked[2] + sKeyPacked[3];

    /*WG: scan the keys within their bit chunks and get total sums*/
		uint sumPacked;
    uint rankPacked = localPrefixSum
    (
      key4, 
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
      wi_id, 
      wf_id, 
#else
      wi_id, 
#endif 
      &sumPacked, 
      cache 
    );

		GROUP_EVENTS_LOCAL_BARRIER_WF

		uint newOffset[GROUP_EVENTS_ELEMENTS_PER_WI] = { 0,0,0,0 };
		{
      /*WI: scan my private keys within their bit chunks*/
			{
				uint sum4 = 0;
        uint tmp = sKeyPacked[0];
        sKeyPacked[0] = sum4;
        sum4 += tmp;
        tmp = sKeyPacked[1];
        sKeyPacked[1] = sum4;
        sum4 += tmp;
        tmp = sKeyPacked[2];
        sKeyPacked[2] = sum4;
        sum4 += tmp;
        tmp = sKeyPacked[3];
        sKeyPacked[3] = sum4;
        sum4 += tmp;
			}
      
      /*WI: scan the chunks of sum: */
      
      /*For a WG with 64 WIs and 4 keys per WI: bits [7-0] are all zeros, bits [15-8] hold sum of 0s, 
        bits [23-16] hold sum of 0s and 1s, bits [31-24] hold sum of 0s, 1s and 2s.*/
      uint sumScanned = (sumPacked << 8) + 
                        (sumPacked << 8*2) + 
                        (sumPacked << 8*3);

      /*WI: compute offsets*/
			{
				/*WI: compute offsets for my keys within their bit chunks. Sum provides offset between
          global key groups (0,1,2,3). Rank provides offset between private key groups*/
        uint sumPlusRank = sumScanned + rankPacked;
        
        /*WI: compute private offsets within my private key groups*/
				{	uint ie = b.x;
          sKeyPacked[0] += sumPlusRank;
					newOffset[0] = (sKeyPacked[0] >> (ie*8)) & 0xFF;
				}
				{	uint ie = b.y;
					sKeyPacked[1] += sumPlusRank;
					newOffset[1] = (sKeyPacked[1] >> (ie*8)) & 0xFF;
				}
				{	uint ie = b.z;
					sKeyPacked[2] += sumPlusRank;
					newOffset[2] = (sKeyPacked[2] >> (ie*8)) & 0xFF;
				}
				{	uint ie = b.w;
					sKeyPacked[3] += sumPlusRank;
					newOffset[3] = (sKeyPacked[3] >> (ie*8)) & 0xFF;
				}
			}
		}

		GROUP_EVENTS_LOCAL_BARRIER_WF

    /*WF: exchange keys and values*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
    newOffset[0] += wfOffset;
    newOffset[1] += wfOffset;
    newOffset[2] += wfOffset;
    newOffset[3] += wfOffset;
#endif
    cache[newOffset[0]] = dataKeys[0];
    cache[newOffset[1]] = dataKeys[1];
    cache[newOffset[2]] = dataKeys[2];
    cache[newOffset[3]] = dataKeys[3];
    GROUP_EVENTS_LOCAL_BARRIER_WF

    dataKeys[0] = cache[dstAddr+0];
    dataKeys[1] = cache[dstAddr+1];
    dataKeys[2] = cache[dstAddr+2];
    dataKeys[3] = cache[dstAddr+3];
    GROUP_EVENTS_LOCAL_BARRIER_WF
    
#if GROUP_EVENTS_VALUES_MODE
    cache[newOffset[0]] = dataValues[0];
    cache[newOffset[1]] = dataValues[1];
    cache[newOffset[2]] = dataValues[2];
    cache[newOffset[3]] = dataValues[3];
    GROUP_EVENTS_LOCAL_BARRIER_WF

    dataValues[0] = cache[dstAddr+0];
    dataValues[1] = cache[dstAddr+1];
    dataValues[2] = cache[dstAddr+2];
    dataValues[3] = cache[dstAddr+3];
    GROUP_EVENTS_LOCAL_BARRIER_WF
#endif 
	}
}
#endif



/*
  group_events
  
  WG loads global offsets for its source events. These offsets define placement of this WG 
  events in the global for all WGs data structure of globally sorted events. For this:
  each WF iterates over source events from a current time slot located in GM in chunks of 
  WF x GROUP_EVENTS_ELEMENTS_PER_WI (each WI picks GROUP_EVENTS_ELEMENTS_PER_WI items):
    - loads each chunk
    - sorts within a chunk by bit-masked keys, specifically by GROUP_EVENTS_HISTOGRAM_BIN_BITS 
      bits in keys starting from bit GROUP_EVENTS_HISTOGRAM_BIT_SHIFT
    - computes local offsets (AKA pointers, histogram) of these keys using scan
    - computes global offsets for events in current chunk based on local offets and WG-level 
      global offsets (AKA scan, histogram, pointers) loaded at the beginning
    - stores keys (and optionally values) at global offsets
    - increments WG-level global offsets before proceeding to the next chunk
  
  Result: globally sorted events (sorted by GROUP_EVENTS_HISTOGRAM_BIN_BITS bits in keys starting 
          from bit position GROUP_EVENTS_HISTOGRAM_BIT_SHIFT
*/
__kernel
void group_events
(
#if (GROUP_EVENTS_DEBUG_ENABLE)
  __global      uint          *gm_debug_host,
  __global      uint          *gm_debug_device,
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  __global      uint          *gm_error_code,
#endif
#if GROUP_EVENTS_VALUES_MODE && (GROUP_EVENTS_RELOCATE_VALUES || GROUP_EVENTS_REPLACE_KEY)
  __global      uint          *gm_event_targets,
  __global      uint          *gm_event_delays,
  __global      uint          *gm_event_weights,
#endif
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
	__global      uint          *gm_histogram_out,
#endif
#if (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 0)
  __global      uint          *gm_source_event_counts,
#endif
  __global      uint          *gm_source_events,
  __global      uint          *gm_destination_events,
  __global      uint          *gm_offsets,
                uint          step
){
  __local uint cache[GROUP_EVENTS_CACHE_SIZE_WORDS];
  
	__local uint lmlocalHistogramReference[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS];
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT && GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT)
	__local uint lmGroupEventsHistogramOut[GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT];
#endif

	uint wi_id = get_local_id(0);
	uint wg_id = get_group_id(0);
  
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  uint local_wf_id  = (wi_id/GROUP_EVENTS_WF_SIZE_WI);
  uint wi_id_wf_scope = (wi_id%GROUP_EVENTS_WF_SIZE_WI);
#else
  #define local_wf_id             0
  #define wi_id_wf_scope          wi_id
#endif

#if (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 0)
  /*Load event totals*/
  /*step is a current time slot, which is due for execution: */
  if(wi_id_wf_scope == 0)
  {
    /* Get WG total events in a temporary LM location*/
    uint index_ptr = GROUP_EVENTS_TIME_SLOTS * wg_id + step;
    cache[0] = *(gm_source_event_counts + index_ptr);
    /*Get global total events in a temporary LM location*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    uint offset = 
      step*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE+1) + 
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;
    cache[1] = gm_offsets[offset];
#endif
  }

  /*Load global bin bit offsets for data allocated to this WG */
	if(wi_id < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*WI: load histogram item for my WG from a bin mapped to my wi_id*/
    uint offset = 
      /*Time slot*/
      step*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE+1) + 
      /*WG offset*/
      wg_id + 
      /*WI pitch*/
      wi_id*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

    lmlocalHistogramReference[wi_id] = gm_offsets[offset];
	}

  /*Init histogram for the next scan/sort depending where it is placed: LDS or GM/L1,L2*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
#if ((GROUP_EVENTS_WG_SIZE_WI < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*\
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT)) || !GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS)
  {
#if !GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
    uint ptr = wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT);
#endif
    for
    (
      uint i = wi_id; 
      i < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); 
      i += GROUP_EVENTS_WG_SIZE_WI
    ){
#if !GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
      gm_histogram_out[ptr + i] = 0; 
#else
      lmGroupEventsHistogramOut[i] = 0;
#endif
    }
  }
#else
  if(wi_id < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))
  {
#if !GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
    gm_histogram_out[wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + wi_id] = 0; 
#else
      lmGroupEventsHistogramOut[wi_id] = 0;
#endif
  }
#endif
#endif

  /*Get base data address for WG and time slot*/
  uint base_address = 
    /*Event data buffers*/
    wg_id * GROUP_EVENTS_TIME_SLOTS * GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE +
    /*Current event data buffer*/
    step * GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE;

  /*Broadcast total events to all WIs in a WF*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  barrier(CLK_LOCAL_MEM_FENCE);
#else
  GROUP_EVENTS_LOCAL_BARRIER_WF
#endif
  uint total_wg_synaptic_events = cache[0];
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  uint total_synaptic_events = cache[1];
#endif
  
  /*Compute limit address*/
  uint limit_address = base_address + total_wg_synaptic_events;
  
  /*Compute total chunks for this WG.*/
  uint totalSynapticEventChunksPerWi = (uint)
    ceil(((float)total_wg_synaptic_events)/(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI));
    
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Compute total chunks for grid*/
  uint total_synaptic_event_chunks = (uint)
    ceil(((float)total_synaptic_events)/(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI));

  /*Compute total chunks per WG for output histogram*/
  uint wg_chunk_size = (uint)
    ceil(((float)total_synaptic_event_chunks)/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE);

  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);
#endif

#elif (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 1)

  /* Get totals for the events in a temporary LM location*/
  if(wi_id == 0)
  {
    cache[0] = *(gm_offsets + GROUP_EVENTS_HISTOGRAM_BIN_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS);
  }

  /*Load global bin bit offsets for data allocated to this WG */
	if( wi_id < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS) )
	{
    /*WI: load histogram item for my WG from a bin mapped to my wi_id (bin size is NUM_WGS)*/
    uint offset = 
    /*WG offset*/
    wg_id + 
    /*WI pitch*/
    wi_id*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

    lmlocalHistogramReference[wi_id] = gm_offsets[offset];
	}
  
  /*Init histogram for the next scan/sort depending where it is placed: LDS or GM/L1,L2*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
#if ((GROUP_EVENTS_WG_SIZE_WI < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*\
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT)) || !GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS)
  {
#if !GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
    uint ptr = wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT);
#endif
    for
    (
      uint i = wi_id; 
      i < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); 
      i += GROUP_EVENTS_WG_SIZE_WI
    ){
#if !GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
      gm_histogram_out[ptr + i] = 0; 
#else
      lmGroupEventsHistogramOut[i] = 0;
#endif
    }
  }
#else
  if(wi_id < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))
  {
#if !GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
    gm_histogram_out[wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + wi_id] = 0; 
#else
      lmGroupEventsHistogramOut[wi_id] = 0;
#endif
  }
#endif
#endif

  /*Broadcast total events to all WIs in the group*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
  barrier(CLK_LOCAL_MEM_FENCE);
#else
  GROUP_EVENTS_LOCAL_BARRIER_WF
#endif
  uint limit_address = cache[0];

  /*Compute total chunks for grid*/
  uint total_synaptic_event_chunks = (uint)
    ceil(((float)limit_address)/(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI));
  
  /*Compute total chunks for WG*/
  uint chunks_per_wg = (uint)
    ceil(((float)total_synaptic_event_chunks)/GROUP_EVENTS_GRID_SIZE_WG);
  
  /*Compute WG start and end chunk*/
  uint chunk_start = wg_id*chunks_per_wg;
  uint chunk_end = chunk_start + chunks_per_wg;
  if(chunk_start > total_synaptic_event_chunks){chunk_start = total_synaptic_event_chunks;}
  if(chunk_end > total_synaptic_event_chunks){chunk_end = total_synaptic_event_chunks;}
  
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Compute total chunks per WG for output histogram*/
  uint wg_chunk_size = (uint)
    ceil(((float)total_synaptic_event_chunks)/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE);
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);
#endif
#endif

  GROUP_EVENTS_LOCAL_BARRIER_WF
  
  /*WI: iterate through my blocks*/
#if (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 0)
	/*TODO: for last iterations with WFs > 1 the empty WFs could skip some stuff.*/
  for(uint i = 0; i < totalSynapticEventChunksPerWi; i++)
  {
    uint addr = base_address + i*(GROUP_EVENTS_ELEMENTS_PER_WI*GROUP_EVENTS_WG_SIZE_WI) +
      wi_id*GROUP_EVENTS_ELEMENTS_PER_WI;

    /*WI: load data for my block: keys are data elements, 
    values are pointers to original address of these data elements*/
    uint dataKeys[GROUP_EVENTS_ELEMENTS_PER_WI];
#if GROUP_EVENTS_VALUES_MODE
    uint dataValues[GROUP_EVENTS_ELEMENTS_PER_WI];
#endif
    
		for(uint i=0; i < GROUP_EVENTS_ELEMENTS_PER_WI; i++)
    {
      uint gm_address = addr + i;
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
			dataKeys[i] = ( gm_address < limit_address )? gm_source_events[gm_address] : 0xffffffff;
#else
			dataKeys[i] = gm_source_events[gm_address];
#endif
#if GROUP_EVENTS_VALUES_MODE
      dataValues[i] = gm_address;
#endif
    }
    
#elif (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 1)

	for(uint i = chunk_start; i < chunk_end; i++)
  {
    uint addr = 
      /*WI within chunk*/
      wi_id*GROUP_EVENTS_ELEMENTS_PER_WI + 
      /*Chunk for this WG*/
      i*GROUP_EVENTS_ELEMENTS_PER_WI*GROUP_EVENTS_WG_SIZE_WI;
    
    /*WI: load data for my block: keys are data elements, 
    values are pointers to original address of these data elements*/
    uint dataKeys[GROUP_EVENTS_ELEMENTS_PER_WI];
#if GROUP_EVENTS_VALUES_MODE
    uint dataValues[GROUP_EVENTS_ELEMENTS_PER_WI];
#endif
		for(uint i=0; i<GROUP_EVENTS_ELEMENTS_PER_WI; i++)
    {
      uint gm_address = addr + i;
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
			dataKeys[i] = ( gm_address < limit_address )? gm_source_events[gm_address + 
        GROUP_EVENTS_SOURCE_KEY_OFFSET*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE] : 0xffffffff;
#else
			dataKeys[i] = gm_source_events[gm_address + GROUP_EVENTS_SOURCE_KEY_OFFSET*
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE];
#endif
#if GROUP_EVENTS_VALUES_MODE == 1
      dataValues[i] = gm_address;
#elif GROUP_EVENTS_VALUES_MODE == 2
      dataValues[i] = gm_source_events[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + 
        gm_address];
#endif
    }
#endif /*END: GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE*/

#if (GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT)
    /*
      Sort data based on 4 bits starting with start_bit.
      In order to coalesce the memory writes as much as
      possible the keys are sorted in the LDS. During the sorting the kernel
      performs scattered writes into LDS, but this memory is designed
      for efficient random access. Global memory hosted in DRAM is far
      less suited to scattered writes than LDS
    */
    sort4Bits1
    (
#if (GROUP_EVENTS_DEBUG_ENABLE)
      gm_debug_device, 
#endif
#if GROUP_EVENTS_VALUES_MODE
      dataValues, 
#endif
      dataKeys, 
#if GROUP_EVENTS_ENABLE_STEP_SHIFT
      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT*step,
#else
      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT, 
#endif
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
      wi_id_wf_scope,
      local_wf_id,
#else
      wi_id, 
#endif
      cache
    );
#endif/*GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT*/

    /*Get bin keys for dataKeys based on 4 bits starting with start_bit*/
		uint keys[GROUP_EVENTS_ELEMENTS_PER_WI];
		for(uint i=0; i<GROUP_EVENTS_ELEMENTS_PER_WI; i++)
#if GROUP_EVENTS_ENABLE_STEP_SHIFT
			keys[i] = (dataKeys[i] >> GROUP_EVENTS_HISTOGRAM_BIT_SHIFT*step)&
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#else
			keys[i] = (dataKeys[i] >> GROUP_EVENTS_HISTOGRAM_BIT_SHIFT)&
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#endif

    GROUP_EVENTS_LOCAL_BARRIER_WF
    
    /*Initialize counts to zero. Reuse cache[]*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
#if (GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS == GROUP_EVENTS_WF_SIZE_WI)
    cache[GROUP_EVENTS_WF_CACHE_SIZE_WORDS*local_wf_id + wi_id_wf_scope] = 0;
#elif ((GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS > GROUP_EVENTS_WF_SIZE_WI) || \
      !GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS)
    for
    (
      uint i = wi_id_wf_scope; 
      i < GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS; 
      i+= GROUP_EVENTS_WF_SIZE_WI
    ){
      cache[GROUP_EVENTS_WF_CACHE_SIZE_WORDS*local_wf_id + i] = 0;
    }
#else
    if(wi_id_wf_scope < GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS)
    {
      cache[GROUP_EVENTS_WF_CACHE_SIZE_WORDS*local_wf_id + wi_id_wf_scope] = 0;
    }
#endif

    /*Count key occurence. 
    Every GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER threads share a single counter buffer*/
    uint setIdx = wi_id_wf_scope/GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER;
    for(uint i = 0; i < GROUP_EVENTS_ELEMENTS_PER_WI; i++)
    {
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
      if( addr + i < limit_address )
#endif
      {
        atomic_inc(&cache[GROUP_EVENTS_WF_CACHE_SIZE_WORDS*local_wf_id + 
          setIdx*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + keys[i]]);
      }
    }
    
#else /*(GROUP_EVENTS_WG_SIZE_WF == 1)*/

#if (GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS == GROUP_EVENTS_WF_SIZE_WI)
    cache[wi_id] = 0;
#elif ((GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS > GROUP_EVENTS_WF_SIZE_WI) || \
      !GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS)
    for
    (
      uint i = wi_id; 
      i < GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS; 
      i+= GROUP_EVENTS_WF_SIZE_WI
    ){
      cache[i] = 0;
    }
#else
    if(wi_id < GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS)
    {
      cache[wi_id] = 0;
    }
#endif

    /*Count key occurence. 
    Every GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER threads share a single counter buffer*/
    uint setIdx = wi_id/GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER;
    for(uint i = 0; i < GROUP_EVENTS_ELEMENTS_PER_WI; i++)
    {
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
      if( addr + i < limit_address )
#endif
      {
        atomic_inc(&cache[setIdx*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + keys[i]]);
      }
    }
#endif

    GROUP_EVENTS_LOCAL_BARRIER_WF
    
    /*Compute bin histogram across bin buffers. Each WI computes a total for assigned to it bin*/
    uint myHistogram = 0;
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
    uint hIdx = wi_id_wf_scope - GROUP_EVENTS_HISTOGRAM_TOTAL_BINS;
    /*TODO: need the case if it is larger than WF size.*/
    if(wi_id_wf_scope < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)
    {
      /*Run across bin buffers*/
      for(uint i = 0; i < GROUP_EVENTS_WF_SIZE_WI/GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER; i++)
      {
        /*Add bin value to a total for this bin*/
        myHistogram += cache[GROUP_EVENTS_WF_CACHE_SIZE_WORDS*local_wf_id + 
          i*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id_wf_scope]; 
      }
      /*WI: modify my ptr for histogram initialization*/
      hIdx = GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id_wf_scope;
      
      /*Store histogram for later offsetting.*/
      if(local_wf_id < GROUP_EVENTS_WG_SIZE_WF-1)
      {
        HISTOGRAM(local_wf_id*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id_wf_scope) = myHistogram;
      }
    }
#else
    uint hIdx = wi_id - GROUP_EVENTS_HISTOGRAM_TOTAL_BINS;
    /*TODO: need the case if it is larger than WF size.*/
    if( wi_id < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)
    {
      /*Run across bin buffers*/
      for(uint i = 0; i < GROUP_EVENTS_WF_SIZE_WI/GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER; i++)
      {
        /*Add bin value to a total for this bin*/
        myHistogram += cache[i*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id]; 
      }
      /*WI: modify my ptr for histogram initialization*/
      hIdx = GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id;
    }
#endif

    GROUP_EVENTS_LOCAL_BARRIER_WF
    
    /*Initialize prefix sum with histogram data. Reuse cache[]*/
    /*TODO: need the case if it is larger than WF size.*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
    hIdx += local_wf_id*GROUP_EVENTS_WF_CACHE_SIZE_WORDS;
    
    if(wi_id_wf_scope < 2*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)
    {
      cache[hIdx] = myHistogram;
    }
#else
    if(wi_id < 2*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)
    {
      cache[hIdx] = myHistogram;
    }
#endif

    GROUP_EVENTS_LOCAL_BARRIER_WF

    /*Perform prefix sum on the histogram*/
#if (GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT)
    if(wi_id_wf_scope < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)
    {
#if (GROUP_EVENTS_OPTIMIZATION_2LEVEL_REDUCE)

      cache[hIdx] = cache[hIdx-1];
      mem_fence(CLK_LOCAL_MEM_FENCE);
      
      uint u0, u1, u2;

      /*WI: Add 4 values to the right including mine, stride 1*/
      u0 = cache[hIdx-3];
      u1 = cache[hIdx-2];
      u2 = cache[hIdx-1];
      atomic_add( &cache[hIdx], u0 + u1 + u2 );
      GROUP_EVENTS_LOCAL_MEM_FENCE_WF
      
      /*WI: Add 4 values to the right including mine, stride 4*/
      u0 = cache[hIdx-12];
      u1 = cache[hIdx-8];
      u2 = cache[hIdx-4];
      atomic_add( &cache[hIdx], u0 + u1 + u2 );
      GROUP_EVENTS_LOCAL_MEM_FENCE_WF
      
#else /*!GROUP_EVENTS_OPTIMIZATION_2LEVEL_REDUCE*/

      cache[hIdx] = cache[hIdx-1];
      mem_fence(CLK_LOCAL_MEM_FENCE);
      cache[hIdx] += cache[hIdx-1];
      mem_fence(CLK_LOCAL_MEM_FENCE);
      cache[hIdx] += cache[hIdx-2];
      mem_fence(CLK_LOCAL_MEM_FENCE);
      cache[hIdx] += cache[hIdx-4];
      mem_fence(CLK_LOCAL_MEM_FENCE);
      cache[hIdx] += cache[hIdx-8];
      mem_fence(CLK_LOCAL_MEM_FENCE);
#endif
    }
#endif /*END: GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT*/

    /*Barier*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
    barrier(CLK_LOCAL_MEM_FENCE);
#else
    GROUP_EVENTS_LOCAL_BARRIER_WF
#endif

    /*Offset prefix sum*/
    /*OFFSET() is the same as cache[], but changed to int type*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
    if((local_wf_id > 0) && (wi_id_wf_scope < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS))
    {
      uint offset = HISTOGRAM(wi_id_wf_scope);
#if (GROUP_EVENTS_WG_SIZE_WF > 2)
      if(local_wf_id > 1)
      {
        offset += HISTOGRAM(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id_wf_scope);
      }
#if (GROUP_EVENTS_WG_SIZE_WF > 3)
      if(local_wf_id > 2)
      {
        offset += HISTOGRAM(2*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id_wf_scope);
      }
#endif
#endif
#if (GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT)
      OFFSET(local_wf_id*GROUP_EVENTS_WF_CACHE_SIZE_WORDS + GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + 
        wi_id_wf_scope) -= offset;
#else
      cache[local_wf_id*GROUP_EVENTS_WF_CACHE_SIZE_WORDS + GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + 
        wi_id_wf_scope] += offset;
#endif
    }
    GROUP_EVENTS_LOCAL_BARRIER_WF
#else
    GROUP_EVENTS_LOCAL_BARRIER_WF
#endif

    /*Store data based on global and local offsets*/
    /*WI: Iterate through my data*/
    for(uint ie = 0; ie < GROUP_EVENTS_ELEMENTS_PER_WI; ie++)
    {
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
      if( addr + ie < limit_address )
#endif
      {
        uint bin = keys[ie];
        uint myIdx = lmlocalHistogramReference[bin];
      
#if (GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT)
        /*WI: compute global ptr for current item = global current bin ptr for this WG + item 
          ID relative to the first item in the current bin. Since items are sorted the prefix
          sum of bin totals is used to offset local item ID and get its relative ID. For a WG
          with multiple WFs the prefix sum of a WF has to be offset with bin counts of all 
          preceeding WFs.*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
        myIdx += GROUP_EVENTS_ELEMENTS_PER_WI*wi_id_wf_scope + ie;
        myIdx -= OFFSET(local_wf_id*GROUP_EVENTS_WF_CACHE_SIZE_WORDS + 
          GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + bin);
#else
        myIdx += GROUP_EVENTS_ELEMENTS_PER_WI*wi_id + ie;
        myIdx -= cache[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + bin];
#endif
#else
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
        myIdx += atomic_dec(&cache[local_wf_id*
          GROUP_EVENTS_WF_CACHE_SIZE_WORDS + GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + bin])-1;
#else
        myIdx += atomic_dec(&cache[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + bin])-1;
#endif
#endif
        /*TODO: optimize*/
#if GROUP_EVENTS_VALUES_MODE
#if GROUP_EVENTS_RELOCATE_VALUES
        gm_destination_events[myIdx] = dataKeys[ie];
        uint value_ptr = dataValues[ie];
        uint value1 = gm_event_delays[value_ptr];
        uint value2 = gm_event_weights[value_ptr];
        gm_destination_events[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + myIdx] = value1;
        gm_destination_events[2*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + myIdx] = value2;
#else /*!GROUP_EVENTS_RELOCATE_VALUES*/
#if GROUP_EVENTS_REPLACE_KEY
        uint new_key_ptr = dataValues[ie];
        uint new_key = gm_event_targets[new_key_ptr];
        gm_destination_events[myIdx] = new_key;
        gm_destination_events[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE+myIdx] = dataValues[ie];
#else /*!GROUP_EVENTS_REPLACE_KEY*/
        gm_destination_events[myIdx] = dataKeys[ie];
        gm_destination_events[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE+myIdx] = dataValues[ie];
#endif /*END: GROUP_EVENTS_REPLACE_KEY*/
#endif  /*END GROUP_EVENTS_RELOCATE_VALUES*/
#else /*!GROUP_EVENTS_VALUES_MODE*/
        gm_destination_events[myIdx] = dataKeys[ie];
#endif /*END GROUP_EVENTS_VALUES_MODE*/

#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
        /*Compute output histogram key for the next scan/sort*/
#if GROUP_EVENTS_REPLACE_KEY
#if GROUP_EVENTS_ENABLE_STEP_SHIFT
        uint hist_out_bin = (new_key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT*(step+1)) & 
          GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT;
#else
        uint hist_out_bin = (new_key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT) & 
          GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT;
#endif
#else /*!GROUP_EVENTS_REPLACE_KEY*/
#if GROUP_EVENTS_ENABLE_STEP_SHIFT
        uint hist_out_bin = (dataKeys[ie]>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT*(step+1)) & 
          GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT;
#else
        uint hist_out_bin = (dataKeys[ie]>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT) & 
          GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT;
#endif
#endif
        uint hist_out = 
        /*WG, which is going to work on this data element in the next scan/sort*/
        myIdx/wg_chunk_size + 
        /*Bin pointer*/
        hist_out_bin*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
        
        /*Increment counter for this neuron: it goes into a radix bin and future WG sub-bin*/
#if GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT
        atomic_inc(&lmGroupEventsHistogramOut[hist_out]);
#else
        atomic_inc(&gm_histogram_out[wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
          GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + hist_out]);
#endif
#endif
      }
    }
    
    /*Barier*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
    barrier(CLK_LOCAL_MEM_FENCE);
#else
    GROUP_EVENTS_LOCAL_BARRIER_WF
#endif
    
    /*Increment offset for next block of data*/
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
		if( wi_id_wf_scope < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
		{
      atomic_add(&lmlocalHistogramReference[wi_id_wf_scope], myHistogram);
    }
#else
		if( wi_id_wf_scope < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
		{
      lmlocalHistogramReference[wi_id_wf_scope] += myHistogram;
    }
#endif
	}
  
  /*Store histogram for the next scan/sort in GM exactly as it is in LM: radix bin is divided into 
    future WG sub-bins*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT && GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT)
#if ((GROUP_EVENTS_WG_SIZE_WI < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*\
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT)) || !GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS)
  {
    uint ptr = wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT);
    for
    (
      uint i = wi_id; 
      i < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); 
      i += GROUP_EVENTS_WG_SIZE_WI
    ){
        gm_histogram_out[ptr + i] = lmGroupEventsHistogramOut[i]; 
    }
  }
#else
  if(wi_id < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))
  {
    gm_histogram_out[wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + wi_id] = lmGroupEventsHistogramOut[wi_id]; 
  }
#endif
#endif
}
