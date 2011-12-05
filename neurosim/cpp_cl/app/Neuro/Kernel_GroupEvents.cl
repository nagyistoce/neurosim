/*
  TODO:

*/



#if (GROUP_EVENTS_LOCAL_SORT_ENABLE)
uint unpack4Key( uint key, int keyIdx ){ return (key>>(keyIdx*8)) & 0xff;}



uint bit8Scan(uint v){return (v<<8) + (v<<16) + (v<<24);}



uint localPrefixSum( uint pData, uint lIdx, uint* totalSum, __local uint sorterSharedMemory[], int wgSize /*64 or 128*/ )
{
	{	//	Set data
		sorterSharedMemory[lIdx] = 0;
		sorterSharedMemory[lIdx+wgSize] = pData;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	{	//	Prefix sum
		int idx = 2*lIdx + (wgSize+1);
#if defined(GROUP_EVENTS_USE_2LEVEL_REDUCE)
		if( lIdx < 64 )
		{
			uint u0, u1, u2;
			u0 = sorterSharedMemory[idx-3];
			u1 = sorterSharedMemory[idx-2];
			u2 = sorterSharedMemory[idx-1];
			atomic_add( &sorterSharedMemory[idx], u0+u1+u2 );			
			mem_fence(CLK_LOCAL_MEM_FENCE);

			u0 = sorterSharedMemory[idx-12];
			u1 = sorterSharedMemory[idx-8];
			u2 = sorterSharedMemory[idx-4];
			atomic_add( &sorterSharedMemory[idx], u0+u1+u2 );			
			mem_fence(CLK_LOCAL_MEM_FENCE);

			u0 = sorterSharedMemory[idx-48];
			u1 = sorterSharedMemory[idx-32];
			u2 = sorterSharedMemory[idx-16];
			atomic_add( &sorterSharedMemory[idx], u0+u1+u2 );			
			mem_fence(CLK_LOCAL_MEM_FENCE);
			if( wgSize > 64 )
			{
				sorterSharedMemory[idx] += sorterSharedMemory[idx-64];
				mem_fence(CLK_LOCAL_MEM_FENCE);
			}

			sorterSharedMemory[idx-1] += sorterSharedMemory[idx-2];
			mem_fence(CLK_LOCAL_MEM_FENCE);
		}
#else
		if( lIdx < 64 )
		{
			sorterSharedMemory[idx] += sorterSharedMemory[idx-1];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			sorterSharedMemory[idx] += sorterSharedMemory[idx-2];			
			mem_fence(CLK_LOCAL_MEM_FENCE);
			sorterSharedMemory[idx] += sorterSharedMemory[idx-4];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			sorterSharedMemory[idx] += sorterSharedMemory[idx-8];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			sorterSharedMemory[idx] += sorterSharedMemory[idx-16];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			sorterSharedMemory[idx] += sorterSharedMemory[idx-32];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			if( wgSize > 64 )
			{
				sorterSharedMemory[idx] += sorterSharedMemory[idx-64];
				mem_fence(CLK_LOCAL_MEM_FENCE);
			}

			sorterSharedMemory[idx-1] += sorterSharedMemory[idx-2];
			mem_fence(CLK_LOCAL_MEM_FENCE);
		}
#endif
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	*totalSum = sorterSharedMemory[wgSize*2-1];
	uint addValue = sorterSharedMemory[lIdx+wgSize-1];
	return addValue;
}



/*2 scan, 2 exchange*/
#if GROUP_EVENTS_VALUES_MODE
void sort4Bits1(uint dataKeys[4], uint dataValues[4], int startBit, int wi_id, __local uint* lmSortData)
#else
void sort4Bits1(uint dataKeys[4], int startBit, int wi_id, __local uint* lmSortData)
#endif
{
  /*ibit = 0, 2, ... GROUP_EVENTS_HISTOGRAM_BIN_BITS-1*/
	for(uint ibit=0; ibit<GROUP_EVENTS_HISTOGRAM_BIN_BITS; ibit+=2)
	{
    /*Right-shift and store current 2-bit value starting from startBit inside uint4*/
    uint4 b = (uint4)((dataKeys[0]>>(startBit+ibit)) & 0x3, 
                      (dataKeys[1]>>(startBit+ibit)) & 0x3, 
                      (dataKeys[2]>>(startBit+ibit)) & 0x3, 
                      (dataKeys[3]>>(startBit+ibit)) & 0x3);

		uint key4;
		uint sKeyPacked[4] = { 0, 0, 0, 0 };
		{
			/*2^ of (0,1,2,3)*8 = (0x1, 0x100, 0x10000,  0x1000000)*/
      sKeyPacked[0] |= 1<<(8*b.x);
			sKeyPacked[1] |= 1<<(8*b.y);
			sKeyPacked[2] |= 1<<(8*b.z);
			sKeyPacked[3] |= 1<<(8*b.w);

			key4 = sKeyPacked[0] + sKeyPacked[1] + sKeyPacked[2] + sKeyPacked[3];
		}

		uint rankPacked;
		uint sumPacked;
		{
			rankPacked = localPrefixSum( key4, wi_id, &sumPacked, lmSortData, GROUP_EVENTS_WG_SIZE_WI );
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		uint newOffset[4] = { 0,0,0,0 };
		{
			uint sumScanned = bit8Scan( sumPacked );

			uint scannedKeys[4];
			scannedKeys[0] = 1<<(8*b.x);
			scannedKeys[1] = 1<<(8*b.y);
			scannedKeys[2] = 1<<(8*b.z);
			scannedKeys[3] = 1<<(8*b.w);
			{	//	4 scans at once
				uint sum4 = 0;
				for(int ie=0; ie<4; ie++)
				{
					uint tmp = scannedKeys[ie];
					scannedKeys[ie] = sum4;
					sum4 += tmp;
				}
			}

			{
				uint sumPlusRank = sumScanned + rankPacked;
				{	uint ie = b.x;
					scannedKeys[0] += sumPlusRank;
					newOffset[0] = unpack4Key( scannedKeys[0], ie );
				}
				{	uint ie = b.y;
					scannedKeys[1] += sumPlusRank;
					newOffset[1] = unpack4Key( scannedKeys[1], ie );
				}
				{	uint ie = b.z;
					scannedKeys[2] += sumPlusRank;
					newOffset[2] = unpack4Key( scannedKeys[2], ie );
				}
				{	uint ie = b.w;
					scannedKeys[3] += sumPlusRank;
					newOffset[3] = unpack4Key( scannedKeys[3], ie );
				}
			}
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		{
			lmSortData[newOffset[0]] = dataKeys[0];
			lmSortData[newOffset[1]] = dataKeys[1];
			lmSortData[newOffset[2]] = dataKeys[2];
			lmSortData[newOffset[3]] = dataKeys[3];

			barrier(CLK_LOCAL_MEM_FENCE);

			uint dstAddr = 4*wi_id;
			dataKeys[0] = lmSortData[dstAddr+0];
			dataKeys[1] = lmSortData[dstAddr+1];
			dataKeys[2] = lmSortData[dstAddr+2];
			dataKeys[3] = lmSortData[dstAddr+3];

			barrier(CLK_LOCAL_MEM_FENCE);
      
#if GROUP_EVENTS_VALUES_MODE
			lmSortData[newOffset[0]] = dataValues[0];
			lmSortData[newOffset[1]] = dataValues[1];
			lmSortData[newOffset[2]] = dataValues[2];
			lmSortData[newOffset[3]] = dataValues[3];

			barrier(CLK_LOCAL_MEM_FENCE);

			dataValues[0] = lmSortData[dstAddr+0];
			dataValues[1] = lmSortData[dstAddr+1];
			dataValues[2] = lmSortData[dstAddr+2];
			dataValues[3] = lmSortData[dstAddr+3];

			barrier(CLK_LOCAL_MEM_FENCE);
#endif 
		}
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
          
  TODO: 
    - detect events with identical time and fuse them into a single event (need to sync with update
      phase on that)
    - currently WG consits of a single WF. Potential direction to overcome this limitation:
      a)  WG takes work of 2 or more WGs, each WF independently iterates over its data. No sync
          is needed, but can perform poorely for uneven data sizes. However, this can be
          metigated by pushing more allocations to each WF for better work balance.
      b)  WG splits work by 2 or more and allocates each big chunk to a separate WF. This would 
          require a sync between WFs.
    - lots of barrier(CLK_LOCAL_MEM_FENCE) in the code. Are they skipped by compiler at WF level?
    - could benefit from vectorisation of GM structures.
    - instead of using atomics on output histogram it could use same method of computing histogram
      as it does for local one.
    - tik-tok buffer size can be reduced by 1/3 since 2nd key space is not used. Key and value
      buffers could be separate ones.
*/
__kernel
/*__attribute__((reqd_work_group_size(GROUP_EVENTS_WG_SIZE_WI,1,1)))*/
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
	__global      uint          *gm_original_values,
#endif
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
	__global      uint          *gm_histogram_out,
#endif
  __global      uint          *gm_source_events,
  __global      uint          *gm_destination_events,
  __global      uint          *gm_offsets,
                uint          step
){
	__local uint lmSortData[GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI+16];
	__local uint lmlocalHistogramReference[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS];
	__local uint lmlocalHistogramScratchPad[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*2];
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
	__local uint lmGroupEventsHistogramOut
    [GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT];
#endif

	uint wi_id = get_local_id(0);
	uint wg_id = get_group_id(0);

#if (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 0)
  /*Number of the current time slot, which is due for execution: */
  uint current_time_slot = (step%GROUP_EVENTS_TIME_SLOTS);
  
  /* Get totals for the events in a temporary LM location*/
  if(wi_id == 0)
  {
    uint index_ptr = GROUP_EVENTS_TIME_SLOTS * wg_id + current_time_slot;
    lmlocalHistogramScratchPad[0] = *(gm_source_events + index_ptr);
  }

  /*Load global bin bit offsets for data allocated to this WG */
	if( wi_id < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS) )
	{
    /*WI: load histogram item for my WG from a bin mapped to my wi_id (bin size is NUM_WGS)*/
    uint offset = 
    /*Time slot*/
    current_time_slot*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE+1) + 
    /*WG offset*/
    wg_id + 
    /*WI pitch*/
    wi_id*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

    lmlocalHistogramReference[wi_id] = gm_offsets[offset];
    
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    /*Get global total events*/
    if(wi_id == 0)
    {
      uint offset = 
      current_time_slot*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE+1) + 
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;
      
      lmlocalHistogramScratchPad[1] = gm_offsets[offset];
    }
#endif
	}
  
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Init histogram for the next scan/sort*/
  for(uint i=wi_id; i<(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); 
    i+= GROUP_EVENTS_WG_SIZE_WI)
  {
    lmGroupEventsHistogramOut[i] = 0;
  }
#endif

  /*Get base data address for WG and time slot*/
  uint base_address = 
    /*Event totals buffers*/
    GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS * GROUP_EVENTS_TIME_SLOTS + 
    /*Event data buffers*/
    wg_id * GROUP_EVENTS_TIME_SLOTS * 
    (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE * GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS) +
    /*Current event data buffer*/
    current_time_slot * 
    (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE * GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS);

  /*Broadcast total events to all WIs in the group*/
  barrier(CLK_LOCAL_MEM_FENCE);
  uint total_wg_synaptic_events = lmlocalHistogramScratchPad[0];
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  uint total_synaptic_events = lmlocalHistogramScratchPad[1];
#endif

  /*Compute limit address*/
  uint limit_address = base_address + total_wg_synaptic_events*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  
  /*Compute total chunks for this WG.*/
  uint total_wg_synaptic_event_chunks = total_wg_synaptic_events/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
    
  if(total_wg_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI) < 
    total_wg_synaptic_events)
  {
    total_wg_synaptic_event_chunks++;
  }
  
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Compute total chunks for grid*/
  uint total_synaptic_event_chunks = total_synaptic_events/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI) < total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  /*Compute total chunks per WG for output histogram*/
  uint wg_chunk_size = total_synaptic_event_chunks/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  if(wg_chunk_size*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);
#endif
#elif (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 1)

  /* Get totals for the events in a temporary LM location*/
  if(wi_id == 0)
  {
    lmlocalHistogramScratchPad[0] = *(gm_offsets + 
      GROUP_EVENTS_HISTOGRAM_BIN_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS);
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
  
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Init histogram for the next scan/sort*/
  for(uint i=wi_id; i<(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); 
    i+= GROUP_EVENTS_WG_SIZE_WI)
  {
    lmGroupEventsHistogramOut[i] = 0;
  }
#endif

  /*Broadcast total events to all WIs in the group*/
  barrier(CLK_LOCAL_MEM_FENCE);
  uint total_synaptic_events = lmlocalHistogramScratchPad[0];

  /*Compute limit address*/
  uint limit_address = (total_synaptic_events)*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  
  /*Compute total chunks for grid*/
  uint total_synaptic_event_chunks = total_synaptic_events/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI) < total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks for WG*/
  uint chunks_per_wg = total_synaptic_event_chunks/GROUP_EVENTS_GRID_SIZE_WG;
  if(chunks_per_wg*GROUP_EVENTS_GRID_SIZE_WG  < total_synaptic_event_chunks)
  {
    chunks_per_wg++;
  }
  
  /*Compute WG start and end chunk*/
  uint chunk_start = wg_id*chunks_per_wg;
  uint chunk_end = chunk_start + chunks_per_wg;
  if(chunk_start > total_synaptic_event_chunks){chunk_start = total_synaptic_event_chunks;}
  if(chunk_end > total_synaptic_event_chunks){chunk_end = total_synaptic_event_chunks;}
  
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Compute total chunks per WG for output histogram*/
  uint wg_chunk_size = total_synaptic_event_chunks/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  if(wg_chunk_size*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);
#endif
#endif

  
  /*TODO: verify if this is really needed here to prevent from racing on lmlocalHistogramScratchPad 
  borrowed temporary above*/
  barrier(CLK_LOCAL_MEM_FENCE);

  /*WI: iterate through my blocks*/
#if (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 0)
	for(uint i = 0; i < total_wg_synaptic_event_chunks; i++)
  {
    uint addr = base_address + 
      wi_id*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*GROUP_EVENTS_ELEMENTS_PER_WI + i*
      (GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*GROUP_EVENTS_ELEMENTS_PER_WI*GROUP_EVENTS_WG_SIZE_WI);
#elif (GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE == 1)
	for(uint i = chunk_start; i < chunk_end; i++)
  {
    uint addr = 
      /*WI within chunk*/
      wi_id*GROUP_EVENTS_ELEMENTS_PER_WI + 
      /*Chunk for this WG*/
      i*GROUP_EVENTS_ELEMENTS_PER_WI*GROUP_EVENTS_WG_SIZE_WI;
    addr *= GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
#endif

		uint myHistogram = 0;

    /*WI: load data for my block: keys are data elements, 
    values are pointers to original address of these data elements*/
    uint dataKeys[GROUP_EVENTS_ELEMENTS_PER_WI];
#if GROUP_EVENTS_VALUES_MODE
    uint dataValues[GROUP_EVENTS_ELEMENTS_PER_WI];
#endif
		for(uint i=0; i<GROUP_EVENTS_ELEMENTS_PER_WI; i++)
    {
      uint gm_address = addr + i*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
			dataKeys[i] = ( gm_address < limit_address )? 
        gm_source_events[gm_address + GROUP_EVENTS_SOURCE_KEY_OFFSET] : 0xffffffff;
#else
			dataKeys[i] = gm_source_events[gm_address + GROUP_EVENTS_SOURCE_KEY_OFFSET];
#endif
#if GROUP_EVENTS_VALUES_MODE == 1
      dataValues[i] = gm_address;
#elif GROUP_EVENTS_VALUES_MODE == 2
      dataValues[i] = gm_source_events[gm_address+1];
#endif
    }

#if (GROUP_EVENTS_LOCAL_SORT_ENABLE)
    /*
      Sort data based on 4 bits starting with start_bit.
      In order to coalesce the memory writes as much as
      possible the keys are sorted in the LDS. During the sorting the kernel
      performs scattered writes into LDS, but this memory is designed
      for efficient random access. Global memory hosted in DRAM is far
      less suited to scattered writes than LDS
    */
    /*TODO: review method sort4Bits1, add handling of values*/
#if GROUP_EVENTS_ENABLE_STEP_SHIFT
#if GROUP_EVENTS_VALUES_MODE
		sort4Bits1(dataKeys, dataValues, GROUP_EVENTS_HISTOGRAM_BIT_SHIFT*step, wi_id, lmSortData);
#else
		sort4Bits1(dataKeys, GROUP_EVENTS_HISTOGRAM_BIT_SHIFT*step, wi_id, lmSortData);
#endif

#else /*!GROUP_EVENTS_ENABLE_STEP_SHIFT*/

#if GROUP_EVENTS_VALUES_MODE
		sort4Bits1(dataKeys, dataValues, GROUP_EVENTS_HISTOGRAM_BIT_SHIFT, wi_id, lmSortData);
#else
		sort4Bits1(dataKeys, GROUP_EVENTS_HISTOGRAM_BIT_SHIFT, wi_id, lmSortData);
#endif
#endif/*GROUP_EVENTS_ENABLE_STEP_SHIFT*/
#endif/*GROUP_EVENTS_LOCAL_SORT_ENABLE*/

    /*Get bin keys for dataKeys based on 4 bits starting with start_bit*/
		uint keys[GROUP_EVENTS_ELEMENTS_PER_WI];
		for(uint i=0; i<GROUP_EVENTS_ELEMENTS_PER_WI; i++)
#if GROUP_EVENTS_ENABLE_STEP_SHIFT
			keys[i] = (dataKeys[i]>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT*step)&
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#else
			keys[i] = (dataKeys[i]>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT)&
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#endif

    /*Create local offsets*/
		{
			/*Initialize offsets to zero*/
      if( wi_id < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
			{
				lmlocalHistogramScratchPad[wi_id] = 0;
			}
      /*Initialize counts to zero*/
			lmSortData[wi_id] = 0;
      
			barrier(CLK_LOCAL_MEM_FENCE);

      uint setIdx = wi_id/GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER;

      /*Count key occurence. 
      Every GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER threads share a single counter buffer*/
			for(uint i=0; i<GROUP_EVENTS_ELEMENTS_PER_WI; i++)
      {
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
				if( addr + i*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS < limit_address )
#endif
        {
          atomic_inc( &lmSortData[(setIdx)*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + keys[i]] );
        }
			}
      
			barrier(CLK_LOCAL_MEM_FENCE);

      /*Compute bin totals across bin buffers*/
			uint hIdx = GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id;
			/*Each WI computes a total for assigned to it bin*/
      if( wi_id < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
			{
				uint sum = 0;
        /*Run across bin buffers*/
				for(uint i=0; i<GROUP_EVENTS_WG_SIZE_WI/GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER; i++)
				{
					/*Add bin value to a total for this bin*/
          sum += lmSortData[i*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id]; 
				}
        /*WI: stores total for assigned to it bin*/
				myHistogram = sum;
				lmlocalHistogramScratchPad[hIdx] = sum;
			}
      
			barrier(CLK_LOCAL_MEM_FENCE);

      /*Perform scan*/
#if (GROUP_EVENTS_LOCAL_SORT_ENABLE)
#if defined(GROUP_EVENTS_USE_2LEVEL_REDUCE)
			if( wi_id < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
			{
				/*Shift to the right by one value*/
        lmlocalHistogramScratchPad[hIdx] = lmlocalHistogramScratchPad[hIdx-1];
				mem_fence(CLK_LOCAL_MEM_FENCE);

        /*WI: Add 4 values to the right including mine, stride 1*/
				uint u0, u1, u2;
				u0 = lmlocalHistogramScratchPad[hIdx-3];
				u1 = lmlocalHistogramScratchPad[hIdx-2];
				u2 = lmlocalHistogramScratchPad[hIdx-1];
				atomic_add( &lmlocalHistogramScratchPad[hIdx], u0 + u1 + u2 );
				mem_fence(CLK_LOCAL_MEM_FENCE);
        
        /*WI: Add 4 values to the right including mine, stride 4*/
				u0 = lmlocalHistogramScratchPad[hIdx-12];
				u1 = lmlocalHistogramScratchPad[hIdx-8];
				u2 = lmlocalHistogramScratchPad[hIdx-4];
				atomic_add( &lmlocalHistogramScratchPad[hIdx], u0 + u1 + u2 );
				mem_fence(CLK_LOCAL_MEM_FENCE);
			}
#else
			if( wi_id < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
			{
				lmlocalHistogramScratchPad[hIdx] = lmlocalHistogramScratchPad[hIdx-1];
				mem_fence(CLK_LOCAL_MEM_FENCE);
				lmlocalHistogramScratchPad[hIdx] += lmlocalHistogramScratchPad[hIdx-1];
				mem_fence(CLK_LOCAL_MEM_FENCE);
				lmlocalHistogramScratchPad[hIdx] += lmlocalHistogramScratchPad[hIdx-2];
				mem_fence(CLK_LOCAL_MEM_FENCE);
				lmlocalHistogramScratchPad[hIdx] += lmlocalHistogramScratchPad[hIdx-4];
				mem_fence(CLK_LOCAL_MEM_FENCE);
				lmlocalHistogramScratchPad[hIdx] += lmlocalHistogramScratchPad[hIdx-8];
				mem_fence(CLK_LOCAL_MEM_FENCE);
			}
#endif
			barrier(CLK_LOCAL_MEM_FENCE);
#endif
		}
    
    /*Store data based on global and local offsets*/
		{
			/*WI: Iterate through my data*/
      for(uint ie=0; ie<GROUP_EVENTS_ELEMENTS_PER_WI; ie++)
			{
				uint bin = keys[ie];
				uint wg_offset_in_global_scope = lmlocalHistogramReference[bin];
#if defined(GROUP_EVENTS_CHECK_BOUNDARY)
				if( addr + ie*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS < limit_address )
#endif
				{
#if (GROUP_EVENTS_LOCAL_SORT_ENABLE)
          uint data_id_in_wg_scope = GROUP_EVENTS_ELEMENTS_PER_WI*wi_id + ie;
          uint myIdx = data_id_in_wg_scope - 
            lmlocalHistogramScratchPad[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + bin];
#else
  #error(option with GROUP_EVENTS_LOCAL_SORT_ENABLE=0 has bugs)
          uint myIdx = 
            atomic_dec(&lmlocalHistogramScratchPad[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + bin])-1;
#endif
          myIdx += wg_offset_in_global_scope;
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
          /*Compute WG, which is going to work on this data element in the next scan/sort*/
          uint hist_out_wg = myIdx/wg_chunk_size;
#endif
          myIdx *= GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;

#if GROUP_EVENTS_VALUES_MODE
#if GROUP_EVENTS_RELOCATE_VALUES
          gm_destination_events[myIdx] = dataKeys[ie];
          uint value_ptr = dataValues[ie];
          uint value1 = gm_original_values[value_ptr+1];
          uint value2 = gm_original_values[value_ptr+2];
          gm_destination_events[myIdx+1] = value1;
          gm_destination_events[myIdx+2] = value2;
#else /*!GROUP_EVENTS_RELOCATE_VALUES*/
#if GROUP_EVENTS_REPLACE_KEY
          uint new_key_ptr = dataValues[ie] + GROUP_EVENTS_REPLACEMENT_KEY_OFFSET;
          uint new_key = gm_original_values[new_key_ptr];
          gm_destination_events[myIdx] = new_key;
          gm_destination_events[myIdx+1] = dataValues[ie];
#else /*!GROUP_EVENTS_REPLACE_KEY*/
          gm_destination_events[myIdx] = dataKeys[ie];
          gm_destination_events[myIdx+1] = dataValues[ie];
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
          /*Increment counter for this neuron.*/
          /*TODO: verify in ISA code that this results in non-returning efficient atomic*/
          atomic_inc(&lmGroupEventsHistogramOut[hist_out_wg + 
            hist_out_bin*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE]);
#endif
        }
			}
		}

		barrier(CLK_LOCAL_MEM_FENCE);

    /*Increment offset for next block of data*/
		if( wi_id < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS )
		{
			lmlocalHistogramReference[wi_id] += myHistogram;
		}
    
		barrier(CLK_LOCAL_MEM_FENCE);
	}
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  /*Store histogram for the next scan/sort. Offset wg_id, pitch GROUP_EVENTS_WG_SIZE_WI: 
  each WI stores an item to to global bin.*/
  for(uint j=0; j<(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE); j++)
	{
#if (GROUP_EVENTS_WG_SIZE_WI < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT)
    for(uint i=wi_id; i<(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); i+=GROUP_EVENTS_WG_SIZE_WI)
    {
      uint ptr = 
        wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
        j + i*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
      gm_histogram_out[ptr] = 
        lmGroupEventsHistogramOut[j + i*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE]; 
    }
#else
    if(wi_id <(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))
    {
      uint ptr = 
        wg_id*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
        j + wi_id*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
      gm_histogram_out[ptr] = 
        lmGroupEventsHistogramOut[j + wi_id*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE]; 
    }
#endif
  }
#endif
}

