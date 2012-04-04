
/*
  GOALS:

*/



#define SCAN_LOCAL_BARRIER_WF                   //barrier(CLK_LOCAL_MEM_FENCE);
#define SCAN_LOCAL_MEM_FENCE_WF                 //mem_fence(CLK_LOCAL_MEM_FENCE);



/*
  Returns data = (uint4)(0, data.x, data.x+data.y, data.x+data.y+data.z) passed by reference 
  and also returns sum = data.x+data.y+data.z+data.w as a value
*/
uint prefixScanVectorEx
( 
  uint4* data
){
	uint sum = 0;
	uint tmp = data[0].x;
	data[0].x = sum;
	sum += tmp;
	tmp = data[0].y;
	data[0].y = sum;
	sum += tmp;
	tmp = data[0].z;
	data[0].z = sum;
	sum += tmp;
	tmp = data[0].w;
	data[0].w = sum;
	sum += tmp;
	return sum;
}



uint localPrefixSum
(
          uint  pData, 
          uint  wi_id, 
          uint* totalSum, 
  __local uint  cache[]
){
#if (SCAN_WG_SIZE_WF > 1)
  uint wf_id  = (wi_id/SCAN_WF_SIZE_WI);
  wi_id = (wi_id%SCAN_WF_SIZE_WI);
  uint wfOffset  = wf_id*SCAN_WF_CACHE_SIZE_WORDS;
  cache[wfOffset + wi_id] = 0;
  wfOffset += (SCAN_WF_SIZE_WI - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET);
  cache[wfOffset + wi_id] = pData;
  uint idx = wfOffset + 2*wi_id + 1;
#else
  cache[wi_id] = 0;
  cache[wi_id + SCAN_WF_SIZE_WI - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] = pData;
  uint idx = 2*wi_id + (SCAN_WF_SIZE_WI + 1) - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET;
#endif

	SCAN_LOCAL_BARRIER_WF
  
  /*Prefix sum*/
#if (SCAN_OPTIMIZATION_2LEVEL_REDUCE)
  /*
  uint u0, u1, u2;
  u0 = cache[idx-3];
  u1 = cache[idx-2];
  u2 = cache[idx-1];
  atomic_add( &cache[idx], u0+u1+u2 );
  SCAN_LOCAL_MEM_FENCE_WF

  u0 = cache[idx-12];
  u1 = cache[idx-8];
  u2 = cache[idx-4];
  atomic_add( &cache[idx], u0+u1+u2 );			
  SCAN_LOCAL_MEM_FENCE_WF

  u0 = cache[idx-48];
  u1 = cache[idx-32];
  u2 = cache[idx-16];
  atomic_add( &cache[idx], u0+u1+u2 );			
  SCAN_LOCAL_MEM_FENCE_WF
  
  cache[idx-1] += cache[idx-2];
  SCAN_LOCAL_MEM_FENCE_WF
  */
  uint u;
  u = cache[idx-3];
  u += cache[idx-2];
  u += cache[idx-1];
  atomic_add( &cache[idx], u );
  SCAN_LOCAL_MEM_FENCE_WF

  u = cache[idx-12];
  u += cache[idx-8];
  u += cache[idx-4];
  atomic_add( &cache[idx], u );
  SCAN_LOCAL_MEM_FENCE_WF

  u = cache[idx-48];
  u += cache[idx-32];
  u += cache[idx-16];
  atomic_add( &cache[idx], u );
  SCAN_LOCAL_MEM_FENCE_WF
  
  cache[idx-1] += cache[idx-2];
  SCAN_LOCAL_MEM_FENCE_WF
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

	SCAN_LOCAL_BARRIER_WF
  
#if (SCAN_WG_SIZE_WF > 1)
  *totalSum = cache[wfOffset - 1 + SCAN_WF_SIZE_WI];
	uint rank = cache[wfOffset - 1 + wi_id];

  barrier(CLK_LOCAL_MEM_FENCE);
  
  /*Offset prefix sum through reduction op*/
  if(wf_id > 0)
  {
    wfOffset = 2*SCAN_WF_SIZE_WI - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET - 1;
    uint offset = cache[(wf_id - 1)*SCAN_WF_CACHE_SIZE_WORDS + wfOffset];
    rank += offset;
    *totalSum += offset;
  }
#if (SCAN_WG_SIZE_WF > 2)
  if(wf_id > 1)
  {
    uint offset = cache[(wf_id - 2)*SCAN_WF_CACHE_SIZE_WORDS + wfOffset];
    rank += offset;
    *totalSum += offset;
  }
#if (SCAN_WG_SIZE_WF > 3)
  if(wf_id > 2)
  {
    uint offset = cache[(wf_id - 3)*SCAN_WF_CACHE_SIZE_WORDS + wfOffset];
    rank += offset;
    *totalSum += offset;
  }
#endif
#endif

#else /*(SCAN_WG_SIZE_WF == 0)*/

	*totalSum = cache[SCAN_WF_SIZE_WI*2 - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
	uint rank = cache[wi_id + SCAN_WF_SIZE_WI - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
#endif

	return rank;
}



__kernel
void
scan_histogram
( 
#if (SCAN_DEBUG_ENABLE)
  __global      uint          *gm_debug_host,
  __global      uint          *gm_debug_device,
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  __global      uint          *gm_error_code,
#endif
  __global      uint          *gm_target_neuron_histogram,
                uint          step
){
  //__local uint cache[SCAN_CACHE_SIZE_WORDS];
  __local uint cache[SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE];
  
  uint wi_id = get_local_id(0);
  uint wg_id = get_group_id(0);
  
	uint data[SCAN_HISTOGRAM_ELEMENTS_PER_WI];
  
  /*WI: load SCAN_HISTOGRAM_ELEMENTS_PER_WI elements into reg space*/
#if SCAN_HISTOGRAM_IN_TYPE == 0
  uint offset = 
    /*Time slot*/
    step*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1) + 
    /*WI work size*/
    SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id;
  
  for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
  {
    data[i] = gm_target_neuron_histogram[offset + i];
  }
#elif 0//SCAN_HISTOGRAM_IN_TYPE == 1
  /*WI: sum-accumulate corresponding elements from my chunk in a current backet*/
  for(uint i = wi_id; i < (SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE); i += SCAN_WG_SIZE_WI)
  {
    cache[i] = gm_target_neuron_histogram[i];
  }
  
  //barrier(CLK_LOCAL_MEM_FENCE);
  //mem_fence(CLK_LOCAL_MEM_FENCE);
  
  /*Iterate over bin backets and sum elements across bin backets using one-to-one mapping*/
  for(uint j=1; j<SCAN_HISTOGRAM_BIN_BACKETS; j++)
  {
    uint offset = j*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE);
    /*WI: sum-accumulate corresponding elements from my chunk in a current backet*/
    for(uint i = wi_id; i < (SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE); i += SCAN_WG_SIZE_WI)
    {
      cache[i] += gm_target_neuron_histogram[offset + i];
    }
    //barrier(CLK_LOCAL_MEM_FENCE);
    //mem_fence(CLK_LOCAL_MEM_FENCE);
  }
  
  for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
  {
    data[i] = 0;
  }
  
  /*WI: sum-accumulate corresponding elements from my chunk in a current backet*/
  for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
  {
    data[i] += cache[SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id + i];
  }
  
#elif SCAN_HISTOGRAM_IN_TYPE == 1
  for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
  {
    data[i] = 0;
  }
  
  uint offsetWi = SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id;
  
  /*Iterate over bin backets and sum elements across bin backets using one-to-one mapping*/
  for(uint j=0; j<SCAN_HISTOGRAM_BIN_BACKETS; j++)
  {
    uint offsetBacket = 
      /*Bin backet consits of radix bins (0-SCAN_HISTOGRAM_TOTAL_BINS). Each bin has 
        SCAN_HISTOGRAM_BIN_SIZE sub-bins. The scan is be done across the elements in the backet
        preserving bin/sub-bin order*/
      j*SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 
      /*WI work size corresponds to some chunk of backet elements*/
      offsetWi;
    
    /*WI: sum-accumulate corresponding elements from my chunk in a current backet*/
    for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
    {
      data[i] += gm_target_neuron_histogram[offsetBacket + i];
    }
  }
#else
  #error (Unrecognized SCAN_HISTOGRAM_IN_TYPE)
#endif

	uint4 myData = (uint4)(0,0,0,0);

  /*Each WI partitions data in 4 chunks of (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) elements and 
  computes their sums*/
	for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
	{
		myData.x += data[i];
		myData.y += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)+i];
		myData.z += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/2)+i];
		myData.w += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*3+i];
	}

	/*Store in myData its scan and return its sum*/
  uint s4 = prefixScanVectorEx( &myData );

  /*WF: compute prefix sum on my chunk*/
  uint sumPacked;
  uint rank = localPrefixSum
  (
    s4, 
    wi_id, 
    &sumPacked, 
    cache 
  );

	SCAN_LOCAL_BARRIER_WF

  /*Propagate this result to sums of smaller chunks of (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) elements*/
	uint4 scanned = myData + (uint4)( rank, rank, rank, rank );

  /*Perform scan on elements within each of 4 chunks with size of 
  (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) elements*/
  /* unroll for(uint j=0; j<4; j++) */
	{
    uint j = 0;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}
	{	uint j = 1;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}
	{	uint j = 2;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}
	{	uint j = 3;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}

  /*Propagate global scan result to the elements within chunks with size of 
  (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) items*/
	for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
	{
		data[i] += scanned.x;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)+i] += scanned.y;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/2)+i] += scanned.z;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*3+i] += scanned.w;
	}

  /*Store scan result to GM*/
  {
    uint offsetBase = 
      /*Time slot*/
      step*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1);
      
    uint offset = 
      /*Time slot*/
      offsetBase + 
      /*WI work size*/
      SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id;
      
    for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
    {
      gm_target_neuron_histogram[offset + i] = data[i];
    }
    
    /*Store total*/
    if(wi_id == SCAN_WG_SIZE_WI-1)
    {
      gm_target_neuron_histogram[offsetBase + SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE] = 
        sumPacked;
    }
  }
}



#if 0
__kernel
void
scan_histogram_old
( 
#if (SCAN_DEBUG_ENABLE)
  __global      uint          *gm_debug_host,
  __global      uint          *gm_debug_device,
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  __global      uint          *gm_error_code,
#endif
  __global      uint          *gm_target_neuron_histogram,
                uint          step
){
  uint wi_id = get_local_id(0);
  uint wg_id = get_group_id(0);
  
	__local uint cache[SCAN_WG_SIZE_WI*2 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
  
	uint data[SCAN_HISTOGRAM_ELEMENTS_PER_WI];
  
  /*WI: load SCAN_HISTOGRAM_ELEMENTS_PER_WI elements into reg space*/
#if SCAN_HISTOGRAM_IN_TYPE == 0
  uint offset = 
    /*Time slot*/
    step*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1) + 
    /*WI work size*/
    SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id;
  
  for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
  {
    data[i] = gm_target_neuron_histogram[offset + i];
  }
#elif SCAN_HISTOGRAM_IN_TYPE == 1
  for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
  {
    data[i] = 0;
  }
  
  uint offsetWi = SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id;
  
  /*Iterate over bin backets and sum elements across bin backets using one-to-one mapping*/
  for(uint j=0; j<SCAN_HISTOGRAM_BIN_BACKETS; j++)
  {
    uint offsetBacket = 
      /*Bin backet consits of radix bins (0-SCAN_HISTOGRAM_TOTAL_BINS). Each bin has 
        SCAN_HISTOGRAM_BIN_SIZE sub-bins. The scan is be done across the elements in the backet
        preserving bin/sub-bin order*/
      j*SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 
      /*WI work size corresponds to some chunk of backet elements*/
      offsetWi;
    
    /*WI: sum-accumulate corresponding elements from my chunk in a current backet*/
    for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
    {
      data[i] += gm_target_neuron_histogram[offsetBacket + i];
    }
  }
#else
  #error (Unrecognized SCAN_HISTOGRAM_IN_TYPE)
#endif

	uint4 myData = (uint4)(0,0,0,0);

  /*Each WI partitions data in 4 chunks of (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) elements and 
  computes their sums*/
	for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
	{
		myData.x += data[i];
		myData.y += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)+i];
		myData.z += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/2)+i];
		myData.w += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*3+i];
	}

	/*Store in myData its scan and return its sum*/
  uint s4 = prefixScanVectorEx( &myData );

  /*Each WF takes (2*SCAN_WF_SIZE_WI) of LM. Init it: 1/2 with 0s and 1/2 with totals of WI's 
  data chunks*/
#if ((SCAN_WG_SIZE_WF == 1) || (SCAN_WG_SIZE_WF == 2))
		cache[wi_id] = 0;
		cache[wi_id + SCAN_WG_SIZE_WI - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] = s4;
#else
	{
    uint offset = wi_id + (2*SCAN_WF_SIZE_WI)*(wi_id/((2*SCAN_WF_SIZE_WI)));
    /*Example for WG with 4 WFs: 0 for WIs 0-127, 128 for WIs 128-255*/
		cache[offset] = 0;
    /*Example for WG with 4 WFs: 128 for 0-127, 256 for 128-255*/
		cache[offset + (2*SCAN_WF_SIZE_WI) - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] = s4;
	}
#endif

  SCAN_LOCAL_BARRIER_WF

/*Perform scan on totals of data chunks*/ 
#if (SCAN_OPTIMIZATION_2LEVEL_REDUCE)
#if(SCAN_WG_SIZE_WF == 1)
  uint idx = 2*wi_id + (SCAN_WG_SIZE_WI + 1) - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET;
  
#elif(SCAN_WG_SIZE_WF == 2)
  uint idx = 2*wi_id + (SCAN_WG_SIZE_WI + 1) - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET;

  if( wi_id < SCAN_WF_SIZE_WI )
#elif(SCAN_WG_SIZE_WF == 4)
  /*idx = (SCAN_WG_SIZE_WI+1), (SCAN_WG_SIZE_WI+3), (SCAN_WG_SIZE_WI+5)... (2*SCAN_WG_SIZE_WI-1)*/
  uint idx = 2*(wi_id + SCAN_WF_SIZE_WI) + 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET;

  if( wi_id < SCAN_WF_SIZE_WI || (wi_id >= (2*SCAN_WF_SIZE_WI) && wi_id < (3*SCAN_WF_SIZE_WI)))
#endif
  {
    uint u0, u1, u2;
    
    /*Scan 3 elements to the right with pitch 1 and add them to my element*/
    u0 = cache[idx-3];
    u1 = cache[idx-2];
    u2 = cache[idx-1];
    atomic_add( &cache[idx], u0+u1+u2 ); /*TODO: are these atomics really needed? Everyone 
    gets its own*/
    mem_fence(CLK_LOCAL_MEM_FENCE);

    /*Scan 3 elements to the right with pitch 4 and add them to my element*/
    u0 = cache[idx-12];
    u1 = cache[idx-8];
    u2 = cache[idx-4];
    atomic_add( &cache[idx], u0+u1+u2 );
    mem_fence(CLK_LOCAL_MEM_FENCE);

    /*Scan 3 elements to the right with pitch 16 and add them to my element*/
    u0 = cache[idx-48];
    u1 = cache[idx-32];
    u2 = cache[idx-16];
    atomic_add( &cache[idx], u0+u1+u2 );
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
#if ( SCAN_WG_SIZE_WF > 1 )
    cache[idx] += cache[idx-SCAN_WF_SIZE_WI];
    mem_fence(CLK_LOCAL_MEM_FENCE);
#endif

    /*Add sums of my right neibour since it was not added above*/
    cache[idx-1] += cache[idx-2];
    
    mem_fence(CLK_LOCAL_MEM_FENCE);
  }
#else
  if( wi_id < SCAN_WF_SIZE_WI )
  {
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
    
#if ( SCAN_WG_SIZE_WF > 1 )
      cache[idx] += cache[idx-SCAN_WF_SIZE_WI];
      mem_fence(CLK_LOCAL_MEM_FENCE);
#endif

    cache[idx-1] += cache[idx-2];
    mem_fence(CLK_LOCAL_MEM_FENCE);
  }
#endif

	SCAN_LOCAL_BARRIER_WF
  
#if(SCAN_WG_SIZE_WF == 4)
  if(wi_id < (2*SCAN_WF_SIZE_WI))
  {
    uint sum = cache[(4*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
    uint s1 = cache[wi_id + (2*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
    uint s2 = cache[wi_id + (6*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
    cache[wi_id + (6*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] = s2 + sum;
    cache[wi_id + (4*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] = s1;
  }
  
  if(wi_id == 0)
  {
    uint sum = cache[wi_id + (6*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
    cache[(8*SCAN_WF_SIZE_WI) - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET] += sum;
  }
  
  SCAN_LOCAL_BARRIER_WF
#endif

  /*Scan result for sums of chunks (global scan)*/
	uint rank = cache[wi_id + SCAN_WG_SIZE_WI - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
  /*Propagate this result to sums of smaller chunks of (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) elements*/
	uint4 scanned = myData + (uint4)( rank, rank, rank, rank );

  /*Perform scan on elements within each of 4 chunks with size of 
  (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) elements*/
  /* unroll for(uint j=0; j<4; j++) */
	{
    uint j = 0;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}
	{	uint j = 1;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}
	{	uint j = 2;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}
	{	uint j = 3;
		uint sum = 0;
		for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
		{
			uint tmp = data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i];
			data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*j+i] = sum;
			sum += tmp;
		}
	}

  /*Propagate global scan result to the elements within chunks with size of 
  (SCAN_HISTOGRAM_ELEMENTS_PER_WI/4) items*/
	for(uint i=0; i<(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4); i++)
	{
		data[i] += scanned.x;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)+i] += scanned.y;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/2)+i] += scanned.z;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*3+i] += scanned.w;
	}

  /*Store scan result to GM*/
  {
    uint offsetBase = 
      /*Time slot*/
      step*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1);
      
    uint offset = 
      /*Time slot*/
      offsetBase + 
      /*WI work size*/
      SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id;
      
    for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
    {
      gm_target_neuron_histogram[offset + i] = data[i];
    }
    
    /*Store total*/
    if(wi_id == 0)
    {
      gm_target_neuron_histogram[offsetBase + SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE] = 
        cache[2*SCAN_WG_SIZE_WI - 1 - SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET];
    }
  }
}
#endif
