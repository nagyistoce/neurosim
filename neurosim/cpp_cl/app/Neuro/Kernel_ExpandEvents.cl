
/*
  expand_events
  
  WG loads a spike packet, histograms and event counters for all time slots. Each WF picks a
  spike from the spike packet and iterates over all target synapses for the spike: loads a batch
  of synapses, computes event delivery time, computes histogram for target neurons based on
  selected bitmask/bit offset, computes pointer in global space for each synapse and stores 
  synapse at that pointer. At the end modified histograms and event counters are stored back to 
  GM.
  
  TODO:
  - Could benefit from replacing atomics with primitives like scan, sort etc
  - Consider redesigning the payload: only relocate arrival time and reference to synaptic structure.
  - Just-in-time delivery of events instead of circular buffer of time slots -> GM size reduction.
  - Bug. The following combination results in device reset:
    #define 	ENABLE_MASK	 	BIN_16(1000,0000,0000,0000)
    #define 	TOTAL_NEURON_BITS	 	17
    #define 	MAX_SYNAPSES_PER_NEURON	 	(24*64)
    #define   SYNAPSE_DEVIATION_RATIO 0.5
    #define 	EXPAND_EVENTS_TEST_MODE	 	5
    #define 	SPIKE_PACKET_SIZE	 	128
    #define 	SPIKE_PACKETS	 	512
    #define 	EVENT_DATA_BUFFER_SIZE	 	(32*1024)
    #define 	SYNAPTIC_EVENT_BUFFERS	 	64
    #define 	EXPAND_EVENTS_WG_SIZE_WF	 	2
    #define 	EXPAND_EVENTS_SPIKE_BUFFER_SIZE	 	128
    #define 	EXPAND_EVENTS_INTER_WF_COOPERATION	 	0
    #define 	PREINITIALIZE_NETWORK_STATE	 	0
    #define 	TOLERANCE_MODE	 	0

    #define 	SIMULATION_MODE	 	0
    #define 	EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT	 	(16*2+2)
    #define 	SIMULATION_TIME_STEPS	 	(16*2+1)

*/

#undef GLOBAL_TOTAL_WFs
#define GLOBAL_TOTAL_WFs                  (EXPAND_EVENTS_GRID_SIZE_WG*EXPAND_EVENTS_WG_SIZE_WF)

  /*Local memory accessors*/
#if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
  #define HISTOGRAM(i)                    cache[i]
  #define TIME_SLOT_COUNTERS(i)           cache[EXPAND_EVENTS_CACHE_OFFSET_1 + i]
  #define SPIKE_DATA(i)                   cache[EXPAND_EVENTS_CACHE_OFFSET_2 + i]
#else
  #define HISTOGRAM(i)                    cache[i]
  #define TIME_SLOT_COUNTERS(i)           cache[EXPAND_EVENTS_CACHE_OFFSET_1 + i]
  #define TOTAL_SPIKES                    cache[EXPAND_EVENTS_CACHE_OFFSET_2]
  #define SPIKE_DATA(i)                   cache[EXPAND_EVENTS_CACHE_OFFSET_3 + i]
#endif

//#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
__kernel 
void expand_events
(
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  __global      uint          *gm_debug_host,
  __global      uint          *gm_debug_device,
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  __global      uint          *gm_error_code,
#endif
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  __global      uint          *gm_target_neuron_histogram,
#endif
  __global      uint          *gm_spike_counts,
  __global      uint          *gm_spikes,
  __global      uint          *gm_event_counts,
  __global      uint          *gm_event_targets,
  __global      uint          *gm_event_delays,
  __global      uint          *gm_event_weights,
  __global      uint          *gm_synapse_targets,
  __global      DATA_TYPE     *gm_synapse_delays,
  __global      uint          *gm_synapse_weights,
  __global      uint          *gm_synapse_pointer,
                uint          step
){
#if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
  __local uint cache[EXPAND_EVENTS_CACHE_SIZE_WORDS];
  
  uint wi_id            = get_local_id(0);
  uint wg_id            = get_group_id(0);
  uint global_wf_id     = (wg_id*EXPAND_EVENTS_WG_SIZE_WF + wi_id/EXPAND_EVENTS_WF_SIZE_WI);
  uint local_wf_id      = (wi_id/EXPAND_EVENTS_WF_SIZE_WI);
  uint wi_id_wf_scope   = (wi_id%EXPAND_EVENTS_WF_SIZE_WI);

  /*WF: independently load its spike counts*/
#if EXPAND_EVENTS_SPIKE_PACKETS_PER_WF > EXPAND_EVENTS_WF_SIZE_WI || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i = wi_id_wf_scope; 
    i < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF;
    i += EXPAND_EVENTS_WF_SIZE_WI
  ){
    SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS
      + EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id + i) = 
      gm_spike_counts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + i];
  }
#else
  if(wi_id_wf_scope < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF)
  {
    SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS
      + EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id + wi_id_wf_scope) = 
      gm_spike_counts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + wi_id_wf_scope];
  }
#endif

  /*Load event counter from each time slot except the previous one, which has to be reset*/
#if EXPAND_EVENTS_TIME_SLOTS > EXPAND_EVENTS_WF_SIZE_WI || !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for(uint i = wi_id_wf_scope; i < (EXPAND_EVENTS_TIME_SLOTS); i += EXPAND_EVENTS_WF_SIZE_WI)
  {
    uint offset = EXPAND_EVENTS_TIME_SLOTS*local_wf_id;
    TIME_SLOT_COUNTERS(offset + i) = 0;
    if(i != ((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS))
    {
      TIME_SLOT_COUNTERS(offset + i) = gm_event_counts[EXPAND_EVENTS_TIME_SLOTS * global_wf_id + i];
    }
  }
#else
  if(wi_id_wf_scope < EXPAND_EVENTS_TIME_SLOTS)
  {
    uint offset = EXPAND_EVENTS_TIME_SLOTS*local_wf_id;
    TIME_SLOT_COUNTERS(offset + wi_id_wf_scope) = 0;
    if(wi_id_wf_scope != ((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS))
    {
      TIME_SLOT_COUNTERS(offset + wi_id_wf_scope) = 
        gm_event_counts[EXPAND_EVENTS_TIME_SLOTS * global_wf_id + wi_id_wf_scope];
    }
  }
#endif

  /*Load histogram totals for all time slots.*/
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  uint wfOffset1 = (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)*local_wf_id;
#if (EXPAND_EVENTS_WF_SIZE_WI < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)) ||\
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  /*TODO: need only to load EXPAND_EVENTS_TIME_SLOTS-1 bins since one bin is reset*/
	for
  (
    uint i=wi_id_wf_scope; 
    i<(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WF_SIZE_WI
  ){
    /*Offset global_wf_id, pitch GLOBAL_TOTAL_WFs: each WI brings a pc from a bin in a 
    time slot for its WF. i/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    HISTOGRAM(wfOffset1 + i) = gm_target_neuron_histogram[global_wf_id + 
      i/EXPAND_EVENTS_TIME_SLOTS + i*GLOBAL_TOTAL_WFs]; 
	}

	/*Reset previos step histogram space*/
#if (EXPAND_EVENTS_WF_SIZE_WI < EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS) || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i=wi_id_wf_scope; 
    i<(EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WF_SIZE_WI
  ){
    HISTOGRAM(((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS)*
      EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + wfOffset1 + i) = 0;
	}
#else
  if(wi_id_wf_scope < EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)
  {
    HISTOGRAM(((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS)*
      EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + wfOffset1 + wi_id_wf_scope) = 0;
	}
#endif
#else
  uint histogram_total = 0;
	if(wi_id_wf_scope < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*Init local histogram*/
    HISTOGRAM(wfOffset1 + wi_id_wf_scope) = 0;
    /*Offset global_wf_id, pitch GLOBAL_TOTAL_WFs: each WI brings a pc from a bin in a time slot 
    for its WG. wi_id_wf_scope/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    histogram_total = gm_target_neuron_histogram[global_wf_id + 
      wi_id_wf_scope/EXPAND_EVENTS_TIME_SLOTS + wi_id_wf_scope*GLOBAL_TOTAL_WFs];
	}
	/*Reset previos step histogram space*/
  if(wi_id_wf_scope < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS) && 
    ((wi_id_wf_scope/EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS) == ((EXPAND_EVENTS_TIME_SLOTS + 
    (step-1))%EXPAND_EVENTS_TIME_SLOTS)))
  {
		histogram_total = 0;
	}
#endif
#endif

  /*WF: load a batch of spike events produced in update phase*/
  /*TODO: this load algorithm has to be optimized for better workload per WI distribution.
    Besides precaching everything is bad, it waists resources*/
  uint totalSpikes = 0;

  for
  (
    uint i = 0; 
    i < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF;
    i++
  ){
    uint spikeCount = 
      SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS
      + EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id + i);
      
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    if(spikeCount + totalSpikes >= EXPAND_EVENTS_SPIKE_BUFFER_SIZE)
    {
      atomic_or(gm_error_code,EXPAND_EVENTS_ERROR_CODE_2);
      spikeCount = 0;
    }
#endif
#if EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS > EXPAND_EVENTS_WF_SIZE_WI ||\
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
    for
    (
      uint j = wi_id_wf_scope; 
      j < spikeCount*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS;
      j += EXPAND_EVENTS_WF_SIZE_WI
    ){
      SPIKE_DATA((EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS +
        EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id + 
        totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + j) = 
      gm_spikes[EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*
        global_wf_id + i*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS + j];
    }
#else
    if(wi_id_wf_scope < spikeCount*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS)
    {
      SPIKE_DATA((EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS +
        EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id + 
        totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + wi_id_wf_scope) = 
      gm_spikes[EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*
        global_wf_id + i*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS + wi_id_wf_scope];
    }
#endif
    totalSpikes += spikeCount;
  }

  /*Load synapse pointer and synapse count*/
  if(totalSpikes < EXPAND_EVENTS_WG_SIZE_WI)
  {
    if(wi_id_wf_scope < totalSpikes)
    {
      uint wfOffset = (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS +
        EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id;
        
      uint spikedNrnId = 
        SPIKE_DATA(wfOffset + wi_id_wf_scope*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS);
      SPIKE_DATA(wfOffset + wi_id_wf_scope*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS) = 
        gm_synapse_pointer[spikedNrnId];
      SPIKE_DATA(wfOffset + totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + wi_id_wf_scope) = 
        gm_synapse_pointer[spikedNrnId + 1];
    }
  }
  else
  {
    uint wfOffset = (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS +
      EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id;
        
    for(uint i = wi_id_wf_scope; i < totalSpikes; i += EXPAND_EVENTS_WG_SIZE_WI)
    {
      uint spikedNrnId = SPIKE_DATA(wfOffset + i*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS);
      SPIKE_DATA(wfOffset + i*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS) = 
        gm_synapse_pointer[spikedNrnId];
      SPIKE_DATA(wfOffset + totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + i) = 
        gm_synapse_pointer[spikedNrnId + 1];
    }
  }

  /*Each WF picks a spike*/
  /*TODO: instead of predetermined do first come first served using an atimic -> better balancing*/
  for(uint j = 0; j < totalSpikes; j++)
  {
    uint wfOffset = (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS +
      EXPAND_EVENTS_SPIKE_COUNT_PITCH)*local_wf_id;
    DATA_TYPE spikeTime = 
      as_float(SPIKE_DATA(wfOffset + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*j + 1));
    uint  synapsePointer = SPIKE_DATA(wfOffset + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*j);
    uint  synapseCount = 
      SPIKE_DATA(wfOffset + totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + j) - 
      synapsePointer;

    /*Each WF retrieves all targets of the source nrn: */
    for(uint tgt = wi_id_wf_scope; tgt < synapseCount; tgt += EXPAND_EVENTS_WF_SIZE_WI)
    {
      /*Get synapse data*/
      uint synapseOffset  = (synapsePointer + tgt);
      uint target_neuron  = gm_synapse_targets[synapseOffset];
      uint weight         = gm_synapse_weights[synapseOffset];
      DATA_TYPE delay     = gm_synapse_delays[synapseOffset];
      /*Add delay to spike time, account for time step transition by decrementing 
       (step was incremented right before this kernel)*/
      DATA_TYPE event_time = spikeTime + delay; event_time = event_time - 1.0f;
      /*Make it relative to its time slot*/
      DATA_TYPE event_time_binned = event_time - (int)event_time;
      /*Events with 0.0 time bounce back to the previous time slot*/
      /*TODO: get rid of it, possibly needs changes in update phase*/
      if(event_time_binned == 0.0f){event_time_binned = 1.0f;}
      uint bin_correction = (int)event_time_binned;
      /*Calculate time slot*/
      uint time_slot = ((step + (int)event_time - bin_correction) % EXPAND_EVENTS_TIME_SLOTS);
      
      /*Obtain offsets for storing synaptic event*/
      /*TODO: this is really bad and may be not needed since we are doing WF-based access, 
      where it can be mapped to WF ID*/
      uint event_local_ptr = atomic_inc(&TIME_SLOT_COUNTERS(EXPAND_EVENTS_TIME_SLOTS*local_wf_id + 
        time_slot));
      
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
      if(event_local_ptr >= EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
      {
        atomic_or(gm_error_code,EXPAND_EVENTS_ERROR_CODE_1);
        event_local_ptr = (EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE-1);
      }
#endif

      uint event_global_ptr = 
        /*Event data buffers*/
        global_wf_id*EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event data buffer*/
        time_slot*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event*/
        event_local_ptr;

      /*Store the event*/
      uint event_time_uint = as_uint(event_time_binned);
      gm_event_targets[event_global_ptr] = target_neuron; 
      gm_event_delays[event_global_ptr] = event_time_uint;
      gm_event_weights[event_global_ptr] = weight;
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
      /*Compute histogram key for target neuron based on MSBs*/
      uint bin = (event_time_uint>>EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT) & 
        EXPAND_EVENTS_HISTOGRAM_BIN_MASK;
      /*Increment counter for this neuron.*/
      /*TODO: verify in assembly that this results in non-returning efficient atomic*/
      uint wfHistOffset = (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)*local_wf_id;
      atomic_inc(&HISTOGRAM(wfHistOffset + time_slot*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + bin));
#endif
    }
  }

  /*Store event counter to each time slot*/
#if EXPAND_EVENTS_TIME_SLOTS > EXPAND_EVENTS_WF_SIZE_WI || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i = wi_id_wf_scope; 
    i < (EXPAND_EVENTS_TIME_SLOTS); 
    i += EXPAND_EVENTS_WF_SIZE_WI
  ){
    gm_event_counts[EXPAND_EVENTS_TIME_SLOTS*global_wf_id + i] = 
      TIME_SLOT_COUNTERS(EXPAND_EVENTS_TIME_SLOTS*local_wf_id + i);
  }
#else
  if(wi_id_wf_scope < EXPAND_EVENTS_TIME_SLOTS)
  {
    gm_event_counts[EXPAND_EVENTS_TIME_SLOTS*global_wf_id + wi_id_wf_scope] = 
      TIME_SLOT_COUNTERS(EXPAND_EVENTS_TIME_SLOTS*local_wf_id + wi_id_wf_scope);
  }
#endif

  /*Store histogram totals for all time slots.*/
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  uint wfOffset2 = (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)*local_wf_id;
#if (EXPAND_EVENTS_WF_SIZE_WI < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)) || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i=wi_id_wf_scope; 
    i<(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WF_SIZE_WI
  ){
    /*Offset global_wf_id, pitch GLOBAL_TOTAL_WFs: each WI stores a pc to a bin in a time slot.
    i/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    gm_target_neuron_histogram[global_wf_id + i/EXPAND_EVENTS_TIME_SLOTS + i*GLOBAL_TOTAL_WFs] = 
      HISTOGRAM(wfOffset2 + i); 
	}
#else
	if(wi_id_wf_scope <(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*Offset global_wf_id, pitch GLOBAL_TOTAL_WFs: each WI stores a pc to a bin in a time slot. 
    wi_id_wf_scope/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    gm_target_neuron_histogram[global_wf_id + wi_id_wf_scope/EXPAND_EVENTS_TIME_SLOTS + 
      wi_id_wf_scope*GLOBAL_TOTAL_WFs] = HISTOGRAM(wfOffset2 + wi_id_wf_scope) + 
      histogram_total;
	}
#endif
#endif

#else /*EXPAND_EVENTS_INTER_WF_COOPERATION > 0*/

  __local uint cache[EXPAND_EVENTS_CACHE_SIZE_WORDS];
  
  uint wi_id            = get_local_id(0);
  uint wg_id            = get_group_id(0);
  uint global_wf_id     = (wg_id*EXPAND_EVENTS_WG_SIZE_WF + wi_id/EXPAND_EVENTS_WF_SIZE_WI);
  uint local_wf_id      = (wi_id/EXPAND_EVENTS_WF_SIZE_WI);
  uint wi_id_wf_scope   = (wi_id%EXPAND_EVENTS_WF_SIZE_WI);

  if(wi_id == 0){TOTAL_SPIKES = 0;}
  /*Use this barrier or fence if there are none below
  barrier(CLK_LOCAL_MEM_FENCE);*/

  /*WF: independently load its spike counts*/
#if EXPAND_EVENTS_SPIKE_PACKETS_PER_WF > EXPAND_EVENTS_WF_SIZE_WI || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i = wi_id_wf_scope; 
    i < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF;
    i += EXPAND_EVENTS_WF_SIZE_WI
  ){
    SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*local_wf_id + i) = 
      gm_spike_counts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + i];
  }
#else
  if(wi_id_wf_scope < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF)
  {
    SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*local_wf_id + wi_id_wf_scope) = 
      gm_spike_counts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + wi_id_wf_scope];
  }
#endif

  /*Load event counter from each time slot except the previous one, which has to be reset*/
#if EXPAND_EVENTS_TIME_SLOTS > EXPAND_EVENTS_WG_SIZE_WI || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for(uint i = wi_id; i < (EXPAND_EVENTS_TIME_SLOTS); i += EXPAND_EVENTS_WG_SIZE_WI)
  {
    TIME_SLOT_COUNTERS(i) = 0;
    if(i != ((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS))
      TIME_SLOT_COUNTERS(i) = gm_event_counts[EXPAND_EVENTS_TIME_SLOTS * wg_id + i];
  }
#else
  if(wi_id < EXPAND_EVENTS_TIME_SLOTS)
  {
    TIME_SLOT_COUNTERS(wi_id) = 0;
    if(wi_id != ((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS))
      TIME_SLOT_COUNTERS(wi_id) = gm_event_counts[EXPAND_EVENTS_TIME_SLOTS * wg_id + wi_id];
  }
#endif

  /*Load histogram totals for all time slots.*/
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
#if (EXPAND_EVENTS_WG_SIZE_WI < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)) || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  /*TODO: need only to load EXPAND_EVENTS_TIME_SLOTS-1 bins since one bin is reset*/
	for
  (
    uint i=wi_id; 
    i<(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WG_SIZE_WI
  ){
    /*Offset wg_id, pitch EXPAND_EVENTS_GRID_SIZE_WG: each WI brings a pc from a bin in a time slot 
    for its WG. i/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    HISTOGRAM(i) = gm_target_neuron_histogram[wg_id + i/EXPAND_EVENTS_TIME_SLOTS + 
      i*EXPAND_EVENTS_GRID_SIZE_WG]; 
	}

	/*Reset previos step histogram space*/
  /*TODO: combine this barrier with some other barrier*/
  barrier(CLK_LOCAL_MEM_FENCE);
#if (EXPAND_EVENTS_WG_SIZE_WI < EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS) || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i=wi_id; 
    i<(EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WG_SIZE_WI
  ){
    HISTOGRAM(((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS)*
      EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + i) = 0;
	}
#else
  if(wi_id < EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)
  {
    HISTOGRAM(((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS)*
      EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + wi_id) = 0;
	}
#endif
#else
  uint histogram_total = 0;
	if(wi_id < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*Init local histogram*/
    HISTOGRAM(wi_id) = 0;
    /*Offset wg_id, pitch EXPAND_EVENTS_GRID_SIZE_WG: each WI brings a pc from a bin in a time slot 
    for its WG. wi_id/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    histogram_total = gm_target_neuron_histogram[wg_id + wi_id/EXPAND_EVENTS_TIME_SLOTS +   
      wi_id*EXPAND_EVENTS_GRID_SIZE_WG];
	}
  
	/*Reset previos step histogram space*/
  if(wi_id < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS) && 
  ((wi_id/EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS) == ((EXPAND_EVENTS_TIME_SLOTS + 
  (step-1))%EXPAND_EVENTS_TIME_SLOTS)))
  {
		histogram_total = 0;
	}
#endif
#endif

  /*WF: load a batch of spike events produced in update phase*/
  /*TODO: this load algorithm has to be optimized for better workload per WI distribution*/
  for
  (
    uint i = 0; 
    i < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF;
    i++
  ){
    uint spikeCount = 
      SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*local_wf_id + i);

    /*First WI in a WF computes and broadcasts ptr*/
    /*TODO: instead of doing each count atomically here a WF could compute its total and then 
      sync with other WFs using a single atomic.*/
    if(wi_id_wf_scope == 0)
    {
      SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
        EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF + local_wf_id) = 
        atomic_add(&TOTAL_SPIKES, spikeCount);
    }
    uint spikePtr = 
      SPIKE_DATA(EXPAND_EVENTS_SPIKE_BUFFER_SIZE*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
      EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF + local_wf_id);
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    if(spikeCount + spikePtr >= EXPAND_EVENTS_SPIKE_BUFFER_SIZE)
    {
      atomic_or(gm_error_code,EXPAND_EVENTS_ERROR_CODE_2);
      spikeCount = 0;
    }
#endif
#if EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS > EXPAND_EVENTS_WF_SIZE_WI || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
    for
    (
      uint j = wi_id_wf_scope; 
      j < spikeCount*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS;
      j += EXPAND_EVENTS_WF_SIZE_WI
    ){
      SPIKE_DATA(spikePtr*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + j) = 
        gm_spikes[EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*
        global_wf_id + i*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS + j];
    }
#else
    if(wi_id_wf_scope < spikeCount*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS)
    {
      SPIKE_DATA(spikePtr*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + wi_id_wf_scope) = 
        gm_spikes[EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*
        EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + i*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS + 
        wi_id_wf_scope];
    }
#endif
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  uint  totalSpikes = TOTAL_SPIKES;
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  if(totalSpikes >= EXPAND_EVENTS_SPIKE_BUFFER_SIZE)
  {
    atomic_or(gm_error_code,EXPAND_EVENTS_ERROR_CODE_2);
    totalSpikes = EXPAND_EVENTS_SPIKE_BUFFER_SIZE;
  }
#endif

  /*Load synapse pointer and synapse count*/
  if(totalSpikes < EXPAND_EVENTS_WG_SIZE_WI)
  {
    if(wi_id < totalSpikes)
    {
      uint spikedNrnId = SPIKE_DATA(EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*wi_id);
      SPIKE_DATA(EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*wi_id) = gm_synapse_pointer[spikedNrnId];
      SPIKE_DATA(totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + wi_id) = 
        gm_synapse_pointer[spikedNrnId + 1];
    }
  }
  else
  {
    for(uint i = wi_id; i < totalSpikes; i += EXPAND_EVENTS_WG_SIZE_WI)
    {
      uint spikedNrnId = SPIKE_DATA(EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*i);
      SPIKE_DATA(EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*i) = gm_synapse_pointer[spikedNrnId];
      SPIKE_DATA(totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + i) = 
        gm_synapse_pointer[spikedNrnId + 1];
    }
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);

  /*Each WF picks a spike*/
  /*TODO: instead of predetermined do first come first served using an atimic -> better balancing*/
  for(uint j = local_wf_id; j < totalSpikes; j += EXPAND_EVENTS_WG_SIZE_WF)
  {
    DATA_TYPE spikeTime = as_float(SPIKE_DATA(EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*j + 1));
    uint  synapsePointer = SPIKE_DATA(EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*j);
    uint  synapseCount = SPIKE_DATA(totalSpikes*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + j) - 
      synapsePointer;

    /*Each WF retrieves all targets of the source nrn: */
    for(uint tgt = wi_id_wf_scope; tgt < synapseCount; tgt += EXPAND_EVENTS_WF_SIZE_WI)
    {
      /*Get synapse data*/
      uint synapseOffset  = (synapsePointer + tgt);
      uint target_neuron  = gm_synapse_targets[synapseOffset];
      uint weight         = gm_synapse_weights[synapseOffset];
      DATA_TYPE delay     = gm_synapse_delays[synapseOffset];
      /*Add delay to spike time, account for time step transition by decrementing 
       (step was incremented right before this kernel)*/
      DATA_TYPE event_time = spikeTime + delay; event_time = event_time - 1.0f;
      /*Make it relative to its time slot*/
      DATA_TYPE event_time_binned = event_time - (int)event_time;
      /*Events with 0.0 time bounce back to the previous time slot*/
      /*TODO: get rid of it, possibly needs changes in update phase*/
      if(event_time_binned == 0.0f){event_time_binned = 1.0f;}
      uint bin_correction = (int)event_time_binned;
      /*Calculate time slot*/
      uint time_slot = ((step + (int)event_time - bin_correction) % EXPAND_EVENTS_TIME_SLOTS);
      
      /*Obtain offsets for storing synaptic event*/
      /*TODO: this is really bad and may be not needed since we are doing WF-based access, 
      where it can be mapped to WF ID*/
      uint event_local_ptr = atomic_inc(&TIME_SLOT_COUNTERS(time_slot));
      
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
      if(event_local_ptr >= EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
      {
        atomic_or(gm_error_code,EXPAND_EVENTS_ERROR_CODE_1);
        event_local_ptr = (EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE-1);
      }
#endif

      uint event_global_ptr = 
        /*Event data buffers*/
        wg_id*EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event data buffer*/
        time_slot*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event*/
        event_local_ptr;

      /*Store the event*/
      uint event_time_uint = as_uint(event_time_binned);
      gm_event_targets[event_global_ptr] = target_neuron; 
      gm_event_delays[event_global_ptr] = event_time_uint;
      gm_event_weights[event_global_ptr] = weight;
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
      /*Compute histogram key for target neuron based on MSBs*/
      uint bin = (event_time_uint>>EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT) & 
        EXPAND_EVENTS_HISTOGRAM_BIN_MASK;
      /*Increment counter for this neuron.*/
      /*TODO: verify in assembly that this results in non-returning efficient atomic*/
      atomic_inc(&HISTOGRAM(time_slot*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + bin));
#endif
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  /*Store event counter to each time slot*/
#if EXPAND_EVENTS_TIME_SLOTS > EXPAND_EVENTS_WG_SIZE_WI || !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
  for
  (
    uint i = wi_id; 
    i < (EXPAND_EVENTS_TIME_SLOTS); 
    i += EXPAND_EVENTS_WG_SIZE_WI
  ){
    gm_event_counts[EXPAND_EVENTS_TIME_SLOTS*wg_id + i] = TIME_SLOT_COUNTERS(i);
  }
#else
  if(wi_id < EXPAND_EVENTS_TIME_SLOTS)
  {
    gm_event_counts[EXPAND_EVENTS_TIME_SLOTS*wg_id + wi_id] = TIME_SLOT_COUNTERS(wi_id);
  }
#endif

  /*Store histogram totals for all time slots.*/
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
#if (EXPAND_EVENTS_WG_SIZE_WI < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)) || \
    !EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS
	for
  (
    uint i=wi_id; 
    i<(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WG_SIZE_WI
  ){
    /*Offset wg_id, pitch EXPAND_EVENTS_GRID_SIZE_WG: each WI stores a pc to a bin in a time slot.
    i/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    gm_target_neuron_histogram[wg_id + i/EXPAND_EVENTS_TIME_SLOTS + i*EXPAND_EVENTS_GRID_SIZE_WG] = 
      HISTOGRAM(i); 
	}
#else
	if(wi_id <(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*Offset wg_id, pitch EXPAND_EVENTS_GRID_SIZE_WG: each WI stores a pc to a bin in a time slot. 
    wi_id/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    gm_target_neuron_histogram[wg_id + wi_id/EXPAND_EVENTS_TIME_SLOTS + 
      wi_id*EXPAND_EVENTS_GRID_SIZE_WG] = 
      HISTOGRAM(wi_id) + histogram_total;
	}
#endif
#endif
#endif
}
