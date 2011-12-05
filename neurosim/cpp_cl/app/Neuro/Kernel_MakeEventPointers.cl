/*
  TODO:
  
    - need to verify if the following restructuring can increase perf: load, detect overlapped
      boundaries, detect private boundaries and compact private data, scan, store to GM from
      private (or store to LM and then bulk-store to to GM)
    - may need to replace a banch of WI_ID_WF_SCOPE(wi_id) and other ones used in the code 
      with variables calculated ones.
    - instead of count/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS need to operate on elements of
      MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS
*/

//TODO: map MAKE_EVENT_PTRS_.. to the ones w/o

#undef GLOBAL_WF_ID
#define GLOBAL_WF_ID(wg_id, wi_id)\
  (wg_id*MAKE_EVENT_PTRS_WG_SIZE_WF + wi_id/MAKE_EVENT_PTRS_WF_SIZE_WI)
  
#undef LOCAL_WF_ID
#define LOCAL_WF_ID(wi_id)          (wi_id/MAKE_EVENT_PTRS_WF_SIZE_WI)

#undef GLOBAL_TOTAL_WFs
#define GLOBAL_TOTAL_WFs            (MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF)

#undef ELEMENTS_PER_WF
#define ELEMENTS_PER_WF             (MAKE_EVENT_PTRS_WF_SIZE_WI*\
                                     MAKE_EVENT_PTRS_ELEMENTS_PER_WI)

#undef WI_ID_WF_SCOPE
#define WI_ID_WF_SCOPE(wi_id)       (wi_id%MAKE_EVENT_PTRS_WF_SIZE_WI)



#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
__kernel
void glue_event_pointers
(
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  __global      uint          *gm_debug_host,
  __global      uint          *gm_debug_device,
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  __global      uint          *gm_error_code,
#endif
  __global      uint          *gm_event_pointers
){

  __local uint lm_generic[GLUE_EVENT_WF_LM_SHARE_SIZE];
  
	uint wi_id = get_local_id(0);
	uint wg_id = get_group_id(0);
  uint wi_id_wf_scope = (wi_id%GLUE_EVENT_PTRS_WF_SIZE_WI);
  uint local_wf_id = (wi_id/GLUE_EVENT_PTRS_WF_SIZE_WI);

  for
  (
    uint i = wi_id_wf_scope; 
    i < MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF; 
    i += GLUE_EVENT_PTRS_WF_SIZE_WI
  ){
    uint gm_offset = MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      i*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      
    uint valSource = gm_event_pointers[gm_offset+1];
    lm_generic[i] = valSource;
  }
  
  for
  (
    uint i = wi_id_wf_scope; 
    i < MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF; 
    i += GLUE_EVENT_PTRS_WF_SIZE_WI
  ){
    uint gm_offset = MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF*
      MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      i*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      
    uint keyDest = gm_event_pointers[gm_offset];
    uint valDest = gm_event_pointers[gm_offset+1];
    if( keyDest == MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE ){continue;}
    
    for
    (
      uint j = i; 
      j < MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF; 
      j++
    ){
      uint valSource = lm_generic[j];
      if(valSource != MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE)
      {
        gm_event_pointers[keyDest*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE] = valDest;
        gm_event_pointers[keyDest*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1] = 
          (valSource-valDest)/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
        break;
      }
    }
  }
  
  for
  (
    uint i = wi_id_wf_scope; 
    i < 2*MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF; 
    i += GLUE_EVENT_PTRS_WF_SIZE_WI
  ){
    uint gm_offset = MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      i*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      
    gm_event_pointers[gm_offset]    = MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE;
    gm_event_pointers[gm_offset+1]  = MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE;
  }
}
#endif



__kernel
void make_event_pointers
(
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  __global      uint          *gm_debug_host,
  __global      uint          *gm_debug_device,
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  __global      uint          *gm_error_code,
#endif
  __global      uint          *gm_total_events, /*TODO: get rid of it somehow*/
  __global      uint          *gm_sorted_events,
  __global      uint          *gm_event_pointers
){

  __local uint lm_generic[MAKE_EVENT_PTRS_WG_SIZE_WF*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE];
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  __local uint lm_last_element[MAKE_EVENT_PTRS_WG_SIZE_WF*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE];
#endif
  
	uint wi_id = get_local_id(0);
	uint wg_id = get_group_id(0);

  /* Get totals for the events in a temporary LM location*/
  if(wi_id == 0)
  {
    lm_generic[0] = *(gm_total_events + MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET);
  }

  /*Broadcast total events to all WIs in the WG*/
  barrier(CLK_LOCAL_MEM_FENCE);
  uint total_synaptic_events = lm_generic[0];
  
  /*Compute limit address*/
  uint limit_address = total_synaptic_events*MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;

  /*Compute total chunks in the grid*/
  uint total_synaptic_event_chunks = total_synaptic_events/ELEMENTS_PER_WF;
  if(total_synaptic_event_chunks*ELEMENTS_PER_WF < total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks per WF*/
  uint chunks_per_wf = total_synaptic_event_chunks/GLOBAL_TOTAL_WFs;
  if(chunks_per_wf*GLOBAL_TOTAL_WFs  < total_synaptic_event_chunks)
  {
    chunks_per_wf++;
  }
  
  /*Compute WF start and end chunk*/
  uint chunk_start = GLOBAL_WF_ID(wg_id, wi_id)*chunks_per_wf;
  uint chunk_end = chunk_start + chunks_per_wf;
  if(chunk_start > total_synaptic_event_chunks){chunk_start = total_synaptic_event_chunks;}
  if(chunk_end > total_synaptic_event_chunks){chunk_end = total_synaptic_event_chunks;}
  
  /*May be needed for broadcast above before the next lm_generic reuse below*/
  barrier(CLK_LOCAL_MEM_FENCE);
  
  /*Initialize 1st element boundary. It must be determined as a boundary by default.*/
  if(wg_id == 0 && wi_id == 0)
  {
    lm_generic[0] = 0xffffffff;
  }
  /*WF chunk boundary overlap from preceeding segment (last element of preceeding WF)*/
  else if(WI_ID_WF_SCOPE(wi_id) == 0)
  {
    uint addr = MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS*(
      chunk_start*ELEMENTS_PER_WF - 1);
    lm_generic[LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE] = gm_sorted_events[addr];
  }

#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  /*First WI in a WF initializes inter-WF boundary elements*/
  if(WI_ID_WF_SCOPE(wi_id) == 0)
  {
    uint lm_offset = LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    lm_last_element[lm_offset] = MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE;
    lm_last_element[lm_offset + 1] = MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE;
  }
  
  /*Store limit address at the end of buffer for computing event count for the last element*/
  if(wi_id == MAKE_EVENT_PTRS_WF_SIZE_WI-1 && wg_id == MAKE_EVENT_PTRS_GRID_SIZE_WG-1)
  {
    gm_event_pointers[MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF*
      MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE-1] = limit_address;
  }
  
#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

  /*Compute offset in event buffer for WF*/
  uint globalWfOffset = MAKE_EVENT_PTRS_STRUCT_SIZE*GLOBAL_WF_ID(wg_id, wi_id)+1;

  /*A single WI stores limiting address as the first address in the last (dummy) struct 
  (to determine the size of last pointed event chunk in the case if all structs are full)*/
  if(wg_id == 0 && wi_id == 0)
  {
    gm_event_pointers[MAKE_EVENT_PTRS_STRUCT_SIZE*MAKE_EVENT_PTRS_STRUCTS + 
      MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE] = limit_address;
  }
  /*A WI in each WF stores ptr count and limiting address
  (to determine the size of last pointed event chunk for preceeding pointer)*/
  if((chunk_start == chunk_end) && (WI_ID_WF_SCOPE(wi_id) == 0))
  {
    gm_event_pointers[globalWfOffset-1] = 0;
    gm_event_pointers[globalWfOffset+1] = limit_address;
  }
#endif

  uint c = 0, j = 0, storeFirstElement = 1;
  int sign = -1;
  
  /*Iterate over chunks*/
	for(c = chunk_start, j = 0; c < chunk_end; c++, j++)
  {
    sign *= -1;
    uint addr = MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS*(
      /*Chunk for this WF*/
      c*ELEMENTS_PER_WF + 
      /*WI within chunk*/
      WI_ID_WF_SCOPE(wi_id)*MAKE_EVENT_PTRS_ELEMENTS_PER_WI);

    /*WI: load elementes in private space for detecting boundaries between elements*/
    uint dataKeys[MAKE_EVENT_PTRS_ELEMENTS_PER_WI];
		for(uint i=0; i<MAKE_EVENT_PTRS_ELEMENTS_PER_WI; i++)
    {
      uint gm_address = addr + i*MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
      /*TODO The limit address condition can be met only for the last WF in the loop.
        Its check this way by each WF is not great and could potentially be improved.
        May be with similar to 
        uint overflow = clamp(limit_address-gm_address, 0, 1); If overflow then 1, otherwise 0
        need to think more.*/
			dataKeys[i] = (gm_address < limit_address)? gm_sorted_events[gm_address] : 0xffffffff;
    }

    /*WI: detect boundaries between private elements and count them*/
    uint wiTotalBoundaries = 0;
		for(uint i=1; i<MAKE_EVENT_PTRS_ELEMENTS_PER_WI; i++)
    {
      uint boundary = dataKeys[i] - dataKeys[i-1]; /*Keys must be sorted*/
      /*TODO: try this. uint boundary = abs_diff (dataKeys[i], dataKeys[i-1]);*/
      boundary = min((uint)1, boundary); /*If keys are the same then 0, otherwise 1*/
      /*Detect limit address overflow*/
      uint overflow = 0xffffffff - dataKeys[i];
      overflow = min((uint)1, overflow); /*If overflow then 0, otherwise 1*/
      boundary *= overflow; /*Drop overflows*/
      wiTotalBoundaries += boundary;
    }

    /*WI: detect boundaries between elements of left neighbor (including inter WF) and count them*/
    uint offsetLoad = LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE + 
      MAKE_EVENT_PTRS_WF_SIZE_WI*(j%2) + sign*(WI_ID_WF_SCOPE(wi_id));
    /*Same as alternating addresses:
      even: LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE + 
        WI_ID_WF_SCOPE(wi_id)
      
      odd:  LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE + 
        MAKE_EVENT_PTRS_WF_SIZE_WI - WI_ID_WF_SCOPE(wi_id)
    */
    uint offsetStore = offsetLoad + sign;
    /*Same as alternating addresses:
      even: LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE + 
        WI_ID_WF_SCOPE(wi_id) + 1
        
      odd:  LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE + 
        MAKE_EVENT_PTRS_WF_SIZE_WI - WI_ID_WF_SCOPE(wi_id) - 1
    */

    lm_generic[offsetStore] = dataKeys[MAKE_EVENT_PTRS_ELEMENTS_PER_WI - 1];
    uint leftNeighbor = lm_generic[offsetLoad];
    uint interWiBoundary = dataKeys[0] - leftNeighbor;
    /*TODO: try this. uint interWiBoundary = abs_diff(dataKeys[0], leftNeighbor);*/
    /*Detect key change: if keys are the same then 0, otherwise 1*/
    interWiBoundary = min((uint)1, interWiBoundary); 
    /*Detect limit address overflow*/
    uint overflow = 0xffffffff - dataKeys[0];
    overflow = min((uint)1, overflow); /*If overflow then 0, otherwise 1*/
    interWiBoundary *= overflow; /*Drop overflows*/
    wiTotalBoundaries += interWiBoundary;

    /*WF: scan boundary counts for my WIs*/ 
    uint generalLmWfOffset = MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE*LOCAL_WF_ID(wi_id);
    uint generalLmWiId = WI_ID_WF_SCOPE(wi_id);

    lm_generic[generalLmWfOffset + generalLmWiId] = 0; //0-63
    lm_generic[generalLmWfOffset + generalLmWiId + MAKE_EVENT_PTRS_WF_SIZE_WI] = 
      wiTotalBoundaries; //64-127

    uint idx = 2*generalLmWiId + (MAKE_EVENT_PTRS_WF_SIZE_WI+1); //65-191

    uint u0, u1, u2;
    u0 = lm_generic[generalLmWfOffset + idx-3];
    u1 = lm_generic[generalLmWfOffset + idx-2];
    u2 = lm_generic[generalLmWfOffset + idx-1];
    atomic_add(&lm_generic[generalLmWfOffset + idx], u0+u1+u2);  //65-191	

    u0 = lm_generic[generalLmWfOffset + idx-12];
    u1 = lm_generic[generalLmWfOffset + idx-8];
    u2 = lm_generic[generalLmWfOffset + idx-4];
    atomic_add(&lm_generic[generalLmWfOffset + idx], u0+u1+u2);  //65-191

    u0 = lm_generic[generalLmWfOffset + idx-48];//17-143
    u1 = lm_generic[generalLmWfOffset + idx-32];
    u2 = lm_generic[generalLmWfOffset + idx-16];
    atomic_add(&lm_generic[generalLmWfOffset + idx], u0+u1+u2);  //65-191

    lm_generic[generalLmWfOffset + idx-1] += lm_generic[generalLmWfOffset + idx-2]; //64-190 += 63-189

    uint totalSum = lm_generic[generalLmWfOffset + MAKE_EVENT_PTRS_WF_SIZE_WI*2 - 1]; //127
    uint wiWfOffset = lm_generic[generalLmWfOffset + generalLmWiId + 
      MAKE_EVENT_PTRS_WF_SIZE_WI-1]; //63 - 126

#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    if(totalSum>0)
    {
    wiWfOffset *= MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;

#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    if(totalSum*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE >= MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE)
    {
      atomic_or(gm_error_code,MAKE_EVENT_PTRS_ERROR_CODE_1);
      totalSum = MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE/MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    }
#endif

    /*WI: store pointer to the item detected on a boundary between neighbor WIs*/
    if(interWiBoundary)
    {
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
      if(wiWfOffset >= MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE)
      {
        wiWfOffset = MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE-MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      }
#endif
      uint offset = generalLmWfOffset + wiWfOffset;
      lm_generic[offset] = dataKeys[0];
      lm_generic[offset + 1] = addr;
      wiWfOffset+=MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    }
  
    /*WI: store the rest of pointers to LM*/
		for(uint i=1; i<MAKE_EVENT_PTRS_ELEMENTS_PER_WI; i++)
    {
      uint gm_address = addr + i*MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
      if(dataKeys[i] != dataKeys[i-1] && gm_address < limit_address )
      {
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
        if(wiWfOffset >= MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE)
        {
          wiWfOffset = MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE-MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
        }
#endif
        uint offset = generalLmWfOffset + wiWfOffset;
        lm_generic[offset] = dataKeys[i];
        lm_generic[offset + 1] = gm_address;
        wiWfOffset+=MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      }
    }
    
    /*First WI in a WF restores pointer and computes count from the last element of last iteration*/
		// TODO: change these to the last WI in a WF
    if(WI_ID_WF_SCOPE(wi_id) == 0)
    {
      /*Store first element in the data chunk for inter-WF boundary handling
      TODO: this is not great. Each WI has to spend a whole register + "if" is executed every time*/
      if(storeFirstElement)
      {
        if(GLOBAL_WF_ID(wg_id, wi_id) > 0)
        {
          uint gm_offset = MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
            (GLOBAL_WF_ID(wg_id, wi_id) - 1)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
          /*uint neuronId = lm_generic[generalLmWfOffset];
          gm_event_pointers[gm_offset] = neuronId;*/
          uint address = lm_generic[generalLmWfOffset + 1];
          gm_event_pointers[gm_offset+1] = address;
        }
        storeFirstElement = 0;
      }
      else
      {
        uint offset = LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
        uint neuronId = lm_last_element[offset];
        uint address = lm_last_element[offset + 1];
        uint count = lm_generic[generalLmWfOffset + 1] - address;
        gm_event_pointers[neuronId*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE] = address;
        gm_event_pointers[neuronId*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1] = 
          count/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
      }
    }
    
    /*WF: compute counts and store pointers to GM*/
		for(uint i=WI_ID_WF_SCOPE(wi_id); i<totalSum-1; i+=MAKE_EVENT_PTRS_WF_SIZE_WI)
    {
      uint offset = generalLmWfOffset + i*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      uint neuronId = lm_generic[offset];
      uint address = lm_generic[offset + 1];
      uint count = lm_generic[offset + MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1] - address;
      gm_event_pointers[neuronId*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE] = address;
      gm_event_pointers[neuronId*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1] = 
        count/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
    }
    
    /*First WI in a WF stores last element for the next iteration*/
		if(WI_ID_WF_SCOPE(wi_id) == 0)
    {
      uint offset = generalLmWfOffset + (totalSum-1)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      uint neuronId = lm_generic[offset];
      uint address = lm_generic[offset + 1];
      offset = LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      lm_last_element[offset] = neuronId;
      lm_last_element[offset + 1] = address;
    }
    
    totalSum *= MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    }
    
    /*Last WI in a WF restores last key for the next iteration*/
    if(WI_ID_WF_SCOPE(wi_id) == MAKE_EVENT_PTRS_WF_SIZE_WI-1)
    {
      lm_generic[offsetStore] = dataKeys[MAKE_EVENT_PTRS_ELEMENTS_PER_WI - 1];  
    }
    
#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

    /*Last WI in a WF restores last key for the next iteration*/
    if(WI_ID_WF_SCOPE(wi_id) == MAKE_EVENT_PTRS_WF_SIZE_WI-1)
    {
      lm_generic[offsetStore] = dataKeys[MAKE_EVENT_PTRS_ELEMENTS_PER_WI - 1];  
    }
    
    if(totalSum>0)
    {
    wiWfOffset *= MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    totalSum *= MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;

    /*WI: store pointer to the item detected on a boundary between neighbor WIs*/
    if(interWiBoundary)
    {
      uint offset = globalWfOffset + wiWfOffset;
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
      uint offsetCheck = globalWfOffset - MAKE_EVENT_PTRS_STRUCT_SIZE*GLOBAL_WF_ID(wg_id, wi_id) + 
        wiWfOffset + 1;
      if(offsetCheck >= MAKE_EVENT_PTRS_STRUCT_SIZE)
      {
        atomic_or(gm_error_code,MAKE_EVENT_PTRS_ERROR_CODE_1);
        offset = MAKE_EVENT_PTRS_STRUCT_SIZE*GLOBAL_WF_ID(wg_id, wi_id) + 
          (MAKE_EVENT_PTRS_STRUCT_SIZE-MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE);
      }
#endif
      gm_event_pointers[offset] = dataKeys[0];
      gm_event_pointers[offset + 1] = addr;
      wiWfOffset+=MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    }
  
    /*WI: store the rest of pointers*/
    // TODO it might be better to store them to local first and then to global
		for(uint i=1; i<MAKE_EVENT_PTRS_ELEMENTS_PER_WI; i++)
    {
      uint gm_address = addr + i*MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
      if(dataKeys[i] != dataKeys[i-1] && gm_address < limit_address )
      {
        uint offset = globalWfOffset + wiWfOffset;
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
      uint offsetCheck = globalWfOffset - MAKE_EVENT_PTRS_STRUCT_SIZE*GLOBAL_WF_ID(wg_id, wi_id) + 
        wiWfOffset + 1;
      if(offsetCheck >= MAKE_EVENT_PTRS_STRUCT_SIZE)
      {
        atomic_or(gm_error_code,MAKE_EVENT_PTRS_ERROR_CODE_1);
        offset = MAKE_EVENT_PTRS_STRUCT_SIZE*GLOBAL_WF_ID(wg_id, wi_id) + 
          (MAKE_EVENT_PTRS_STRUCT_SIZE-MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE);
      }
#endif
        gm_event_pointers[offset] = dataKeys[i];
        gm_event_pointers[offset + 1] = gm_address;
        wiWfOffset+=MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      }
    }
    }
    globalWfOffset += totalSum;
#endif
  }
  
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  /*First WI in a WF stores inter-WF boundary pointers*/
  if(j && (WI_ID_WF_SCOPE(wi_id) == 0))
  {
    uint gm_offset = 
      /*Neuron-specific pointers and counters*/
      MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      /*First inter-WF elements*/
      (MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF*
      MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE) +
      /*Last inter-WF elements*/
      GLOBAL_WF_ID(wg_id, wi_id)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    uint lm_offset = LOCAL_WF_ID(wi_id)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    uint neuronId = lm_last_element[lm_offset];
    uint address = lm_last_element[lm_offset + 1];
    gm_event_pointers[gm_offset] = neuronId;
    gm_event_pointers[gm_offset+1] = address;
  }
  
#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

  /*Last WI in a WF that iterated at least once*/
  if(j && (WI_ID_WF_SCOPE(wi_id) == (MAKE_EVENT_PTRS_WG_SIZE_WF-1)))
  {
    /*Stores pointer count*/
    uint offset = MAKE_EVENT_PTRS_STRUCT_SIZE*GLOBAL_WF_ID(wg_id, wi_id);
    gm_event_pointers[offset] = (globalWfOffset - offset - 1)/MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
      uint offsetCheck = (globalWfOffset - offset - 1);
      if(offsetCheck >= MAKE_EVENT_PTRS_STRUCT_SIZE)
      {
        atomic_or(gm_error_code,MAKE_EVENT_PTRS_ERROR_CODE_1);
      }
#endif
  }
#endif
}
