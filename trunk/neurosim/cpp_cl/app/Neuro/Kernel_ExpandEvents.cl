
/*
  expand_events
  
  WG loads a spike packet, histograms and event counters for all time slots. Each WF picks a
  spike from the spike packet and iterates over all target synapses for the spike: loads a batch
  of synapses, computes event delivery time, computes histogram for target neurons based on
  selected bitmask/bit offset, computes pointer in global space for each synapse and stores 
  synapse at that pointer. At the end modified histograms and event counters are stored back to 
  GM.
  
  TODO:
  - Scalability: more than 1 spike packet per WF
    - Have a separate data structure for spike counts
    - Load spike counts, then load spikes into local
    - WF-iterate over spikes, load synaptic data, calculate arrival time etc
  - Could benefit from replacing atomics with primitives like scan, sort etc
  - Redesign payload: only relocate arrival time and reference to synaptic structure.
  - Just-in-time delivery of events instead of circular buffer of time slots -> GM size reduction.
*/


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
  uint wi_id            = get_local_id(0);
  uint wg_id            = get_group_id(0);
  uint global_wf_id     = (wg_id*EXPAND_EVENTS_WG_SIZE_WF + wi_id/EXPAND_EVENTS_WF_SIZE_WI);
  uint local_wf_id      = (wi_id/EXPAND_EVENTS_WF_SIZE_WI);
  uint wi_id_wf_scope   = (wi_id%EXPAND_EVENTS_WF_SIZE_WI);

  /*TODO: combine/reuse local memory where possible to reduce its total size.*/
  __local uint lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF + 1 + 
    EXPAND_EVENTS_WG_SIZE_WF];
  /*TODO: this easts up too much memory, need to optimize (buffered on demand loads).*/
  __local uint lmSpikes[2*EXPAND_EVENTS_WG_SIZE_WF + EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*
    EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF];
  __local uint lmTimeSlotCounters[EXPAND_EVENTS_TIME_SLOTS];
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  __local uint lmTargetNeuronHistogram[EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS];
#endif

  /*WF: independently load its spike counts*/
#if EXPAND_EVENTS_SPIKE_PACKETS_PER_WF > EXPAND_EVENTS_WF_SIZE_WI
  for
  (
    uint i = wi_id_wf_scope; 
    i < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF;
    i += EXPAND_EVENTS_WF_SIZE_WI
  ){
    lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*local_wf_id + i] = 
      gm_spike_counts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + i];
  }
#else
  if(wi_id_wf_scope < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF)
  {
    lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*local_wf_id + wi_id_wf_scope] = 
      gm_spike_counts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + wi_id_wf_scope];
  }
#endif

  /*Load event counter from each time slot except the previous one, which has to be reset*/
#if EXPAND_EVENTS_TIME_SLOTS > EXPAND_EVENTS_WG_SIZE_WI
  for(uint i = wi_id; i < (EXPAND_EVENTS_TIME_SLOTS); i += EXPAND_EVENTS_WG_SIZE_WI)
  {
    lmTimeSlotCounters[i] = 0;
    if(i != ((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS))
      lmTimeSlotCounters[i] = gm_event_counts[EXPAND_EVENTS_TIME_SLOTS * wg_id + i];
  }
#else
  if(wi_id < EXPAND_EVENTS_TIME_SLOTS)
  {
    lmTimeSlotCounters[wi_id] = 0;
    if(wi_id != ((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS))
      lmTimeSlotCounters[wi_id] = gm_event_counts[EXPAND_EVENTS_TIME_SLOTS * wg_id + wi_id];
  }
#endif

#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  /*Load histogram totals for all time slots.*/
#if (EXPAND_EVENTS_WG_SIZE_WI < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
  /*TODO: need only to load EXPAND_EVENTS_TIME_SLOTS-1 bins since one bin is reset*/
	for(uint i=wi_id; i<(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); 
    i+=EXPAND_EVENTS_WG_SIZE_WI)
	{
    /*Offset wg_id, pitch EXPAND_EVENTS_GRID_SIZE_WG: each WI brings a pc from a bin in a time slot 
    for its WG. i/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    lmTargetNeuronHistogram[i] = gm_target_neuron_histogram[wg_id + i/EXPAND_EVENTS_TIME_SLOTS + 
      i*EXPAND_EVENTS_GRID_SIZE_WG]; 
	}
  barrier(CLK_LOCAL_MEM_FENCE);
	/*Reset previos step histogram space*/
  /*TODO: need to get rid of it and the barrier above somehow*/
  for(uint i=wi_id; i<(EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS); i+=EXPAND_EVENTS_WG_SIZE_WI)
	{
    lmTargetNeuronHistogram[((EXPAND_EVENTS_TIME_SLOTS + (step-1))%EXPAND_EVENTS_TIME_SLOTS)*
      EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + i] = 0;
	}
#else
  uint histogram_total = 0;
	if(wi_id < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*Init local histogram*/
    lmTargetNeuronHistogram[wi_id] = 0;
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

  if(wi_id == 0)
  {
    lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF] = 0;
  }
  /*TODO: this barrier here is only because lmSpikeCounts init to 0 above. Can it be replaced with
    mem_fence without if(wi_id == 0) above?*/
  barrier(CLK_LOCAL_MEM_FENCE);

  /*WF: independently load a batch of spike events produced in update phase*/
  /*TODO: this has to be done on a fly, integrated with the next loop. Optimize for 1 WF per WG*/
  for
  (
    uint i = 0; 
    i < EXPAND_EVENTS_SPIKE_PACKETS_PER_WF;
    i++
  ){
    uint spikeCount = lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*local_wf_id + i];

    /*First WI in a WF computes and broadcasts ptr*/
    if(wi_id_wf_scope == 0)
    {
      lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF + 1 + local_wf_id] = 
        atomic_add(&lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF], 
        spikeCount);
    }
    uint spikePtr = 
      lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF + 1 + local_wf_id];

#if EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS > EXPAND_EVENTS_WF_SIZE_WI
    for
    (
      uint j = wi_id_wf_scope; 
      j < spikeCount*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS;
      j += EXPAND_EVENTS_WF_SIZE_WI
    ){
      lmSpikes[2*EXPAND_EVENTS_WG_SIZE_WF + spikePtr*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + j] = 
        gm_spikes[EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*
        global_wf_id + i*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS + j];
    }
#else
    if(wi_id_wf_scope < spikeCount*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS)
    {
      lmSpikes[2*EXPAND_EVENTS_WG_SIZE_WF + spikePtr*EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + 
        wi_id_wf_scope] = gm_spikes[EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS*
        EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*global_wf_id + i*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS + 
        wi_id_wf_scope];
    }
#endif
  }

  barrier(CLK_LOCAL_MEM_FENCE);
  uint  totalSpikes = lmSpikeCounts[EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF];

  /*Each WF picks a spike*/
  /*TODO: instead of predetermined do first come first served using an atimic -> better balancing*/
  for(uint j = local_wf_id; j < totalSpikes; j += EXPAND_EVENTS_WG_SIZE_WF)
  {
    uint  spiked_nrn = lmSpikes[2*EXPAND_EVENTS_WG_SIZE_WF + 
    EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS*j];

    /*Get synapse pointer and synapse count*/
    if(wi_id_wf_scope < 2)
    {
      lmSpikes[2*local_wf_id + wi_id_wf_scope] = gm_synapse_pointer[spiked_nrn + wi_id_wf_scope];
    }
    uint  synapsePointer = lmSpikes[2*local_wf_id];
    uint  tgt_max_cnt = lmSpikes[2*local_wf_id + 1] - synapsePointer;

    /*Each WF retrieves all targets of the source nrn: */
    for(uint tgt = wi_id_wf_scope; tgt < tgt_max_cnt; tgt += EXPAND_EVENTS_WF_SIZE_WI)
    {
      /*Get synapse data*/
      uint synapseOffset  = (synapsePointer + tgt);
      uint target_neuron  = gm_synapse_targets[synapseOffset];
      uint weight         = gm_synapse_weights[synapseOffset];
      DATA_TYPE delay     = gm_synapse_delays[synapseOffset];

      DATA_TYPE spike_time = as_float(lmSpikes[2*EXPAND_EVENTS_WG_SIZE_WF + 
        EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * j + 1]);

      /*Add delay to spike time, account for time step transition by decrementing 
       (step was incremented right before this kernel)*/
      DATA_TYPE event_time = spike_time + delay; event_time = event_time - 1.0f;
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
      uint event_local_ptr = atomic_inc(&lmTimeSlotCounters[time_slot]);
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
      if(event_local_ptr >= EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
      {
        atomic_or(gm_error_code,EXPAND_EVENTS_ERROR_CODE_1);
        event_local_ptr = (EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE-1);
      }
#endif
      uint event_global_ptr = 
        /*Event data buffers*/
        wg_id * EXPAND_EVENTS_TIME_SLOTS * EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event data buffer*/
        time_slot * EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event*/
        event_local_ptr;

      /*Store the event*/
      /*TODO: protect from overflow, vectorize*/
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
      atomic_inc(&lmTargetNeuronHistogram[time_slot*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + bin]);
#endif
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  /*Store event counter to each time slot*/
#if EXPAND_EVENTS_TIME_SLOTS > EXPAND_EVENTS_WG_SIZE_WI
  for
  (
    uint i = wi_id; 
    i < (EXPAND_EVENTS_TIME_SLOTS); 
    i += EXPAND_EVENTS_WG_SIZE_WI
  ){
    gm_event_counts[EXPAND_EVENTS_TIME_SLOTS*wg_id + i] = lmTimeSlotCounters[i];
  }
#else
  if(wi_id < EXPAND_EVENTS_TIME_SLOTS)
  {
    gm_event_counts[EXPAND_EVENTS_TIME_SLOTS*wg_id + wi_id] = lmTimeSlotCounters[wi_id];
  }
#endif

#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
/*Store histogram totals for all time slots.*/
#if (EXPAND_EVENTS_WG_SIZE_WI < (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
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
      lmTargetNeuronHistogram[i]; 
	}
#else
	if(wi_id <(EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS))
	{
    /*Offset wg_id, pitch EXPAND_EVENTS_GRID_SIZE_WG: each WI stores a pc to a bin in a time slot. 
    wi_id/EXPAND_EVENTS_TIME_SLOTS here since each time slot has 
    (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS + 1) elements*/
    gm_target_neuron_histogram[wg_id + wi_id/EXPAND_EVENTS_TIME_SLOTS + 
      wi_id*EXPAND_EVENTS_GRID_SIZE_WG] = 
      lmTargetNeuronHistogram[wi_id] + histogram_total;
	}
#endif
#endif
}
