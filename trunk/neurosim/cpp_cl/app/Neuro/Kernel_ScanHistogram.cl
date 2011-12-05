
/*
  GOALS:

*/


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



/*Time slot data can't be more than WG size x elements of data per WI*/
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
  uint wi_id = get_local_id(0);
  uint wg_id = get_group_id(0);
  
	__local uint lmScanData[SCAN_WG_SIZE_WI*2];
  
  /*TODO: looks like this is in GM, but may be cached. Put it in the reg space?*/
	uint data[SCAN_HISTOGRAM_ELEMENTS_PER_WI];
  
  /*Load histogram. Each WI loads SCAN_HISTOGRAM_ELEMENTS_PER_WI contiguos elements into its private space*/
#if SCAN_HISTOGRAM_IN_TYPE == 0
	for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
	{
		data[i] = 0;
    /*TODO: reduce branching here*/
		if( (SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id + i) < SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE )
    {
      uint offset = 
      /*Time slot*/
      (step%SCAN_HISTOGRAM_TIME_SLOTS)*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1) + 
      /*WI work size*/
      SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id + i;
			data[i] = gm_target_neuron_histogram[offset];
    }
	}
#elif SCAN_HISTOGRAM_IN_TYPE == 1
	for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
	{
		data[i] = 0;
    /*TODO: reduce branching here*/
		if( (SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id + i) < SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE )
    {
      for(uint j=0; j<SCAN_HISTOGRAM_BIN_BACKETS; j++)
      {
        uint offset = 
        /*Backet*/
        j*SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 
        /*WI work size*/
        SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id + i;
        
        data[i] += gm_target_neuron_histogram[offset];
      }
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
		myData.x += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*0+i];
		myData.y += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*1+i];
		myData.z += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*2+i];
		myData.w += data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*3+i];
	}

	/*Store in myData its scan and return its sum*/
  uint s4 = prefixScanVectorEx( &myData );

  /*Each WF takes (2*SCAN_WF_SIZE_WI) of LM. Init it: 1/2 with 0s and 1/2 with totals of WI's data chunks*/
	{
    /*Example for WG with 4 WFs: 0 for WIs 0-127, 128 for WIs 128-255*/
		lmScanData[wi_id + (2*SCAN_WF_SIZE_WI)*(wi_id/((2*SCAN_WF_SIZE_WI)))] = 0;
    /*Example for WG with 4 WFs: 128 for 0-127, 256 for 128-255*/
		lmScanData[wi_id + (2*SCAN_WF_SIZE_WI)*(1 + wi_id/((2*SCAN_WF_SIZE_WI)))] = s4;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

  /*Perform scan on totals of data chunks*/
	{
		/*idx = (SCAN_WG_SIZE_WI+1), (SCAN_WG_SIZE_WI+3), (SCAN_WG_SIZE_WI+5)... (2*SCAN_WG_SIZE_WI-1)*/
    uint idx = 2*(wi_id + SCAN_WF_SIZE_WI) + 1;
    
#if (SCAN_USE_2LEVEL_REDUCE)
#if(SCAN_WG_SIZE_WF == 2)
		if( wi_id < SCAN_WF_SIZE_WI )
#elif(SCAN_WG_SIZE_WF == 4)
    if( wi_id < SCAN_WF_SIZE_WI || (wi_id >= (2*SCAN_WF_SIZE_WI) && wi_id < (3*SCAN_WF_SIZE_WI)))
#endif
		{
      uint u0, u1, u2;
      
      /*Scan 3 elements to the right with pitch 1 and add them to my element*/
			u0 = lmScanData[idx-3];
			u1 = lmScanData[idx-2];
			u2 = lmScanData[idx-1];
			atomic_add( &lmScanData[idx], u0+u1+u2 ); /*TODO: are these atomics really needed? Everyone gets its own*/
			mem_fence(CLK_LOCAL_MEM_FENCE);

      /*Scan 3 elements to the right with pitch 4 and add them to my element*/
			u0 = lmScanData[idx-12];
			u1 = lmScanData[idx-8];
			u2 = lmScanData[idx-4];
			atomic_add( &lmScanData[idx], u0+u1+u2 );
			mem_fence(CLK_LOCAL_MEM_FENCE);

      /*Scan 3 elements to the right with pitch 16 and add them to my element*/
			u0 = lmScanData[idx-48];
			u1 = lmScanData[idx-32];
			u2 = lmScanData[idx-16];
			atomic_add( &lmScanData[idx], u0+u1+u2 );
			mem_fence(CLK_LOCAL_MEM_FENCE);
      
#if ( SCAN_WG_SIZE_WI > SCAN_WF_SIZE_WI )
			{
				lmScanData[idx] += lmScanData[idx-SCAN_WF_SIZE_WI];
				mem_fence(CLK_LOCAL_MEM_FENCE);
			}
#endif

      /*Add sums of my right neibour since it was not added above*/
			lmScanData[idx-1] += lmScanData[idx-2];
      
			mem_fence(CLK_LOCAL_MEM_FENCE);
		}
#else
    if( wi_id < SCAN_WF_SIZE_WI )
		{
			lmScanData[idx] += lmScanData[idx-1];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			lmScanData[idx] += lmScanData[idx-2];			
			mem_fence(CLK_LOCAL_MEM_FENCE);
			lmScanData[idx] += lmScanData[idx-4];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			lmScanData[idx] += lmScanData[idx-8];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			lmScanData[idx] += lmScanData[idx-16];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			lmScanData[idx] += lmScanData[idx-32];
			mem_fence(CLK_LOCAL_MEM_FENCE);
			if( SCAN_WG_SIZE_WI > SCAN_WF_SIZE_WI )
			{
				lmScanData[idx] += lmScanData[idx-SCAN_WF_SIZE_WI];
				mem_fence(CLK_LOCAL_MEM_FENCE);
			}

			lmScanData[idx-1] += lmScanData[idx-2];
			mem_fence(CLK_LOCAL_MEM_FENCE);
		}
#endif
	}

	barrier(CLK_LOCAL_MEM_FENCE);
  
#if(SCAN_WG_SIZE_WF == 4)
  if(wi_id < (2*SCAN_WF_SIZE_WI))
  {
    uint sum = lmScanData[(4*SCAN_WF_SIZE_WI) - 1];
    uint s1 = lmScanData[wi_id + (2*SCAN_WF_SIZE_WI) - 1];
    uint s2 = lmScanData[wi_id + (6*SCAN_WF_SIZE_WI) - 1];
    lmScanData[wi_id + (6*SCAN_WF_SIZE_WI) - 1] = s2 + sum;
    lmScanData[wi_id + (4*SCAN_WF_SIZE_WI) - 1] = s1;
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
#endif

  /*Scan result for sums of chunks (global scan)*/
	uint rank = lmScanData[wi_id+SCAN_WG_SIZE_WI-1];
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
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*0+i] += scanned.x;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*1+i] += scanned.y;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*2+i] += scanned.z;
		data[(SCAN_HISTOGRAM_ELEMENTS_PER_WI/4)*3+i] += scanned.w;
	}

  /*Store scan result to GM*/
	for(uint i=0; i<SCAN_HISTOGRAM_ELEMENTS_PER_WI; i++)
	{
    /*TODO: reduce branching here*/
		if( (SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id+i) < SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE )
    {
      uint offset = 
      /*Time slot*/
      (step%SCAN_HISTOGRAM_TIME_SLOTS)*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1) + 
      /*WI work size*/
      SCAN_HISTOGRAM_ELEMENTS_PER_WI*wi_id + i;
      gm_target_neuron_histogram[offset] = data[i];
    }
	}
  
  /*Store total*/
  if(wi_id == 0)
  {
    uint offset = 
    /*Time slot*/
    (step%SCAN_HISTOGRAM_TIME_SLOTS)*(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE + 1) + 
    /*WI work size*/
    SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE;
    
    gm_target_neuron_histogram[offset] = lmScanData[2*SCAN_WG_SIZE_WI-1];
  }
}
