
/*
  TODO:
  
  - PS order limit: continue in local memory if it hits register limit.
  - Restructire spike packet: contiguous neuronID part and contiguous spike time part
  - Spike packet per group instead of per WF?
  - Spike handling kernel has to be scaled  differently than update kernel
  - Investigate "missing a spike": It is possible to miss a spike if the membrane potential only 
    has a brief superthreshold excursion and returns to subthreshold values before the end of the 
    subinterval (Hanuschkin et al., http://www.nest-initiative.org).
  - NR divergence: find more robust method
    - Adapted NR, converges at the first threshold crossing of the membrane potential. 
      See D�Haene et al. Accelerating event-driven simulation of spiking neurons with multiple 
      synaptic time constants.
    - D�Haene et al. Fast and exact simulation methods applied on a broad range of neuron models.
    - Interpolation. Morrison et al. Exact subthreshold integration with continuous spike times 
      in discrete time neural network simulations.
*/


#undef GLOBAL_WF_ID
#define GLOBAL_WF_ID(wg_id, wi_id)\
  (wg_id*UPDATE_NEURONS_WG_SIZE_WF + wi_id/UPDATE_NEURONS_WF_SIZE_WI)
  
#undef LOCAL_WF_ID
#define LOCAL_WF_ID(wi_id)          (wi_id/UPDATE_NEURONS_WF_SIZE_WI)

#undef GLOBAL_WI_ID
#define GLOBAL_WI_ID(wi_id)          (wi_id*UPDATE_NEURONS_WG_SIZE_WI)

#undef WI_ID_WF_SCOPE
#define WI_ID_WF_SCOPE(wi_id)       (wi_id%UPDATE_NEURONS_WF_SIZE_WI)

#undef ELEMENTS_PER_WF
#define ELEMENTS_PER_WF             (UPDATE_NEURONS_WF_SIZE_WI*\
                                     UPDATE_NEURONS_ELEMENTS_PER_WI)
                                     
#undef GLOBAL_TOTAL_WFs
#define GLOBAL_TOTAL_WFs            (UPDATE_NEURONS_GRID_SIZE_WG*UPDATE_NEURONS_WG_SIZE_WF)

#define UPDT_RUNTM_SZ_0                   4
#define UPDT_RUNTM_SZ_1                   (UPDATE_NEURONS_PS_ORDER_LIMIT+1)
#define UPDT_RUNTM_SZ                     (UPDT_RUNTM_SZ_0+UPDT_RUNTM_SZ_1*2)
#define V_RT                              cache[0]
#define U_RT                              cache[1]
#define GA_RT                             cache[2]
#define GG_RT                             cache[3]
#define YP1(i)                            cache[UPDT_RUNTM_SZ_0+i]
#define YP2(i)                            cache[UPDT_RUNTM_SZ_0+UPDT_RUNTM_SZ_1+i]



/*
  Parker-Sochacki integration step for IZ neuron with adaptive order error 
  tollerance test for model variables.
*/
uint gpu_iz_ps_step
(
#if (UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD)
              DATA_TYPE divergence_threshold,
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  __global    uint      *gm_error_code,
#endif
#if UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP > 0
              DATA_TYPE I,
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
              DATA_TYPE tol,
#endif
  __constant  DATA_TYPE *cm_coefficients, //__global    DATA_TYPE const* restrict cm_coefficients,
              DATA_TYPE *cache,
              DATA_TYPE k,
              DATA_TYPE E_ampa,
              DATA_TYPE E_gaba,
              DATA_TYPE E,
              DATA_TYPE a,
              DATA_TYPE b,
              DATA_TYPE dt
){
#if (UPDATE_NEURONS_TOLERANCE_MODE == 1)
  DATA_TYPE tol = UPDATE_NEURONS_PS_TOLERANCE;
#endif

  /*DATA_TYPE tau_ampa_gpu = CONST_FP(cm_coefficients,4), 
           tau_gaba_gpu = CONST_FP(cm_coefficients,5), 
           co;*/
  /*First term*/
  DATA_TYPE vchi = MUL(YP1(0),YP2(0));
  DATA_TYPE u      = U_RT;
  DATA_TYPE g_ampa = GA_RT;
  DATA_TYPE g_gaba = GG_RT;
  
#if UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP > 0
  YP1(1) = MUL( CONST_CO(cm_coefficients,0,0), (vchi - u + MUL(E_ampa, g_ampa) + 
    MUL(E_gaba, g_gaba) + I) );
#else
  YP1(1) = MUL( CONST_CO(cm_coefficients,0,0), (vchi - u + MUL(E_ampa, g_ampa) + 
    MUL(E_gaba, g_gaba)) );
#endif

  u         = MUL( CONST_CO(cm_coefficients,1,0), (MUL( b, YP1(0) ) - u) );
  g_ampa    = MUL( CONST_CO(cm_coefficients,2,0), g_ampa );
  g_gaba    = MUL( CONST_CO(cm_coefficients,3,0), g_gaba );
  YP2(1)    = MUL( k, YP1(1) ) - g_ampa - g_gaba;
  DATA_TYPE dt_pow = dt;
  
  /*Update and test error tolerance on variable value change*/
#if UPDATE_NEURONS_DT_1_0_OPTIMIZATION
  if(dt == 1.0f)
  {
    V_RT  +=  YP1(1);
    U_RT  +=  u;
    GA_RT +=  g_ampa;
    GG_RT +=  g_gaba;
  }
  else
  {
    V_RT  +=  MUL(YP1(1),dt_pow);
    U_RT  +=  MUL(u,dt_pow);
    GA_RT +=  MUL(g_ampa,dt_pow);
    GG_RT +=  MUL(g_gaba,dt_pow);
    dt_pow*=dt;
  }
#else
  V_RT  +=  MUL(YP1(1),dt_pow);
  U_RT  +=  MUL(u,dt_pow);
  GA_RT +=  MUL(g_ampa,dt_pow);
  GG_RT +=  MUL(g_gaba,dt_pow);
  dt_pow*=dt;
#endif

  /*Iterations*/
  uint p;
  for(p=1;p<(UPDATE_NEURONS_PS_ORDER_LIMIT-1);p++)
  {
    /*
    vchi = MUL(YP1(0),YP2(p)) + MUL(YP2(p),YP1(p));
    for(uint i = 1; i < p; i++){vchi+=MUL(YP1(i),YP2(p-i));}
    co = E/((DATA_TYPE)(p+1));
    YP1(p+1) = co*(vchi - u + MUL(E_ampa,g_ampa) + MUL(E_gaba,g_gaba));
    co = a/((DATA_TYPE)(p+1));
    u         = co*(MUL(b,YP1(p)) - u);
    co = -1.0/(tau_ampa_gpu*(DATA_TYPE)(p+1));
    g_ampa    = co*g_ampa;
    co = -1.0/(tau_gaba_gpu*(DATA_TYPE)(p+1));
    g_gaba    = co*g_gaba;
    YP2(p+1) = MUL(k,YP1(p+1)) - g_ampa - g_gaba;
    */

    DATA_TYPE vchi = MUL(YP1(0),YP2(p)) + MUL(YP2(0),YP1(p));
    for(uint i=1; i < p; i++){vchi += MUL(YP1(i),YP2(p-i));}
  
    YP1(p+1) = MUL( CONST_CO(cm_coefficients,0,p), (vchi - u + MUL(E_ampa, g_ampa) + 
      MUL(E_gaba, g_gaba)) );
  
    u         = MUL( CONST_CO(cm_coefficients,1,p), (MUL( b, YP1(p) ) - u) );
    g_ampa    = MUL( CONST_CO(cm_coefficients,2,p), g_ampa );
    g_gaba    = MUL( CONST_CO(cm_coefficients,3,p), g_gaba );
    YP2(p+1) = MUL( k, YP1(p+1) ) - g_ampa - g_gaba;

    /*Update and test error tolerance on variable value change*/
    DATA_TYPE v_old = V_RT;
    DATA_TYPE u_old = U_RT;
    DATA_TYPE ga_old = GA_RT;
    DATA_TYPE gg_old = GG_RT;
    
#if UPDATE_NEURONS_DT_1_0_OPTIMIZATION
    if(dt == 1.0f)
    {
      V_RT  +=  YP1(p+1);
      U_RT  +=  u;
      GA_RT +=  g_ampa;
      GG_RT +=  g_gaba;
    }
    else
    {
      V_RT  +=  MUL(YP1(p+1),dt_pow);
      U_RT  +=  MUL(u,dt_pow);
      GA_RT +=  MUL(g_ampa,dt_pow);
      GG_RT +=  MUL(g_gaba,dt_pow);
      dt_pow*=dt;
    }
#else
    V_RT  +=  MUL(YP1(p+1),dt_pow);
    U_RT  +=  MUL(u,dt_pow);
    GA_RT +=  MUL(g_ampa,dt_pow);
    GG_RT +=  MUL(g_gaba,dt_pow);
    dt_pow*=dt;
#endif

    /*Check for divergence*/
#if (UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD)
    if((fabs(YP1(p+1)) > divergence_threshold))
#else
    if((fabs(YP1(p+1)) > 10.0f))
#endif
    {
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
      atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_6);
#endif
      break;
    }
    
    /*TODO: order variables for tol check according to their statistical rate for satisfying the tol*/
#if(UPDATE_NEURONS_TOLERANCE_MODE == 0)
    DATA_TYPE tol_check = V_RT - v_old;
    if(tol_check)continue;
    tol_check = U_RT - u_old;
    if(tol_check)continue;
    tol_check = GA_RT - ga_old;
    if(tol_check)continue;
    tol_check = GG_RT - gg_old;
    if(tol_check)continue;
#else
    /*TODO: verify if using UPDATE_NEURONS_PS_TOLERANCE directly here is faster than through reg tol*/
    DATA_TYPE tol_check = fabs(V_RT - v_old);
    if(tol_check > tol)continue;
    tol_check = fabs(U_RT - u_old);
    if(tol_check > tol)continue;
    tol_check = fabs(GA_RT - ga_old);
    if(tol_check > tol)continue;
    tol_check = fabs(GG_RT - gg_old);
    if(tol_check > tol)continue;
#endif
    break;
  }
  
  p++;
 
	return p;
}



/*
  Parker-Sochacki update step for IZ neuron for updating variables according
  to the new intial conditions.
*/
/*TODO: need to optimize for 2 cases dt=1 and the rest*/
void gpu_ps_update
(
  //__global    DATA_TYPE const* restrict cm_coefficients,
  __constant  DATA_TYPE   *cm_coefficients,
  
              DATA_TYPE   *cache,
              int         var,
              int         ps_order,
              DATA_TYPE   dt,
              DATA_TYPE   *ynew,
              DATA_TYPE   a,
              DATA_TYPE   b
){
	int i;
  /*DATA_TYPE tau_ampa_gpu = CONST_FP(cm_coefficients,4), 
             tau_gaba_gpu = CONST_FP(cm_coefficients,5);
	DATA_TYPE co_g_ampa_gpu = CONST_FP(cm_coefficients,2), 
            co_g_gaba_gpu = CONST_FP(cm_coefficients,3);*/
            
	DATA_TYPE co_g_ampa_gpu = -1.0f/(UPDATE_NEURONS_TAU_AMPA/UPDATE_NEURONS_DT);
  DATA_TYPE co_g_gaba_gpu = -1.0f/(UPDATE_NEURONS_TAU_GABA/UPDATE_NEURONS_DT);

	/*Calculate first order terms PRINTF("%f vs %f\n",test,*ynew); */
  switch( var ) 
  {
    case 1 :
      YP2(1) = a*(MUL(b,YP1(0)) - YP2(0));	
      break;
    case 2 :
      YP2(1) = co_g_ampa_gpu*YP2(0);
      break;
    case 3 :
      YP2(1) = co_g_gaba_gpu*YP2(0); 
      break;
  }
  
  /*Iterations*/
  for(i=1;i<=(ps_order-1);i++)
  {
    switch( var ) 
    {
      case 1 :
        YP2(i+1) = CONST_CO(cm_coefficients,1,i)*(MUL(b,YP1(i)) - YP2(i));
        break;
      case 2 :
        YP2(i+1) = CONST_CO(cm_coefficients,2,i)*YP2(i);
        break;
      case 3 :
        YP2(i+1) = CONST_CO(cm_coefficients,3,i)*YP2(i);
        break;
    }
  }

	if(dt == 1)
  {
    for(i=1,*ynew=YP2(0);i<=ps_order;i++)
    {
      *ynew+=YP2(i);
    }
  }
	else
  {
    /*Use Horner's method*/
		*ynew = MUL(YP2(ps_order),dt) + YP2(ps_order-1);
		for(i=ps_order-2;i>=0;i--){*ynew = MUL(*ynew,dt) + YP2(i);}
	}
}



/*
  update_neurons
  
  1) Loads model variables, parameters, and synaptic events from GM to PM.
  2) Iterates through synaptic events and updates model variables.
  3) Detects spiking neurons and computes spike times
  4) Writes spike times and updated model variables back to GM
  5) Shifting from NN chunk allocation from per-WF to per-WG basis and implementing WF-based
     consummer model could result in better workload. Global consummer model could be even
     better if overhead is not big.

*/
__kernel 
void update_neurons
(
#if (UPDATE_NEURONS_DEBUG_ENABLE)
  __global      uint            *gm_debug_host,
  __global      uint            *gm_debug_device,
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  __global      uint            *gm_error_code,
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
  __constant    DATA_TYPE       *cm_ps_tol,
#endif
  //__global      DATA_TYPE const* restrict cm_coefficients,
  //constant      DATA_TYPE       *cm_coefficients __attribute__((max_constant_size (4*CONST_SIZE))),
  __constant    DATA_TYPE       *cm_coefficients,
  __global      DATA_TYPE       *gm_model_parameters,
  __global      DATA_TYPE       *gm_model_variables,
  __global      uint            *gm_spike_counts,
  __global      uint            *gm_spikes,
  __global      uint            *gm_event_ptr,
  __global      DATA_TYPE       *gm_events,
                uint            step
){

  uint wi_id = get_local_id(0);
  uint wg_id = get_group_id(0);

  /*Each WF has a spike packet*/
  __local uint lmSpikePacketCounts[UPDATE_NEURONS_WG_SIZE_WF];
  __local uint lmSpikePackets[UPDATE_NEURONS_WG_SIZE_WF*UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS];

  /*Init spike counts*/
  if(WI_ID_WF_SCOPE(wi_id) == 0)
  {
    lmSpikePacketCounts[LOCAL_WF_ID(wi_id)] = 0;
  }

  /*Compute total chunks in the grid*/
  uint totalChunks = UPDATE_NEURONS_TOTAL_NEURONS/ELEMENTS_PER_WF;
  if(totalChunks*ELEMENTS_PER_WF < UPDATE_NEURONS_TOTAL_NEURONS)
  {
    totalChunks++;
  }
  
  /*Compute total chunks per a WF*/
  uint chunksPerWf = totalChunks/GLOBAL_TOTAL_WFs;
  if(chunksPerWf*GLOBAL_TOTAL_WFs  < totalChunks)
  {
    chunksPerWf++;
  }
  
  /*Compute WF start and end chunk*/
  uint wfStartChunk = GLOBAL_WF_ID(wg_id, wi_id)*chunksPerWf;
  uint wfEndChunk = wfStartChunk + chunksPerWf;
  if(wfStartChunk > totalChunks){wfStartChunk = totalChunks;}
  if(wfEndChunk > totalChunks){wfEndChunk = totalChunks;}
  
  /*A WF works on allocated for it chunk of neurons*/
  /*TODO: instead of predetermined do first come first served using an atimic -> better balancing*/
	for(uint j = 0; j < (wfEndChunk-wfStartChunk); j++)
  {
    uint neuronId = (wfStartChunk + j)*ELEMENTS_PER_WF + WI_ID_WF_SCOPE(wi_id);
    
    /*WI: load event pointer and count, reset the count*/
    uint neuronEventPtrStart = *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE);
    int ev_count = *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1);
    ev_count++;
    
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(ev_count >= 0xFFFF)
    {
      atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_8);
      ev_count = 0;
    }
#endif

    /*Load model variables*/
    DATA_TYPE v         = gm_model_variables[neuronId];
    DATA_TYPE u         = gm_model_variables[UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE g_ampa    = gm_model_variables[2*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE g_gaba    = gm_model_variables[3*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    
    /*Load model parameters*/
    DATA_TYPE k         = gm_model_parameters[UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE l         = gm_model_parameters[2*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE b         = gm_model_parameters[3*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE v_peak    = gm_model_parameters[6*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE E_ampa    = gm_model_parameters[7*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE E_gaba    = gm_model_parameters[8*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
#if UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP > 0
    DATA_TYPE I = 0.0f;
    if(step < UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP)
    {
      I = gm_model_parameters[neuronId];
    }
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
    DATA_TYPE ps_tol = cm_ps_tol[neuronId/(UPDATE_NEURONS_TOTAL_NEURONS/UPDATE_NEURONS_TOLERANCE_CHUNKS)];
#endif
    /*
    DATA_TYPE E         = gm_model_parameters[9*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE a         = gm_model_parameters[10*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    */
    DATA_TYPE E         = (1.0f/UPDATE_NEURONS_C)*UPDATE_NEURONS_DT; 
    DATA_TYPE a         = UPDATE_NEURONS_a*UPDATE_NEURONS_DT;

    DATA_TYPE start_gpu = 0.0f, start_dt_part, spiked = -1.0f;
    DATA_TYPE cache[UPDT_RUNTM_SZ];

  while(ev_count > 0)
  {
    /*Same as: uint stp_pass = (ev_count > 1);*/
    uint stp_pass = clamp(ev_count-1, 0, 1);
    
    /*Read event time from valid ptr else make it eq 1: */
    DATA_TYPE in_t_now = 1.0f;
    if(stp_pass)
    {
      in_t_now = gm_events[neuronEventPtrStart + UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE];
    }
    
    /*Compute time difference btwn events: */
    if(spiked < 0)
    {
      start_dt_part = in_t_now - start_gpu;
    }

    if(start_dt_part > 0)
    {
    /*Load parameters into work memory space: */
    YP1(0)  = v;
    V_RT    = v;
    U_RT    = u;
    GA_RT   = g_ampa;
    GG_RT   = g_gaba;
    YP2(0)  = MUL(k,v) - g_ampa - g_gaba + l;

    /*numerical integration step function*/
    uint ps_order = gpu_iz_ps_step
    (
#if (UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD)
      v_peak,
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
      gm_error_code,
#endif
#if UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP > 0
      I,
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
      ps_tol,
#endif
      cm_coefficients, 
      cache,
      k,
      E_ampa,
      E_gaba,
      E,
      a,
      b,
      start_dt_part
    );
    
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(ps_order >= UPDATE_NEURONS_PS_ORDER_LIMIT)
    {
      ps_order = UPDATE_NEURONS_PS_ORDER_LIMIT-1;
      atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_1);
    }
#endif

    DATA_TYPE vnew = V_RT;

    spiked = vnew - v_peak;

    /*if a spike occures: record neuron data and break*/
    if(spiked >= 0)
    {
      *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE) = neuronEventPtrStart;
      *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1) = ev_count;

      /*Record spike and schedule events:*/
      /*TODO: find if storing directly to GM is faster (totals are still accumulated in LM)*/
      uint spikePacketOffset = UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*LOCAL_WF_ID(wi_id);
      uint index = spikePacketOffset + atomic_inc(lmSpikePacketCounts+LOCAL_WF_ID(wi_id))*
        UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS;
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
      if(index > spikePacketOffset + UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS - 
        UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS)
      {
        atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_5);
        index = spikePacketOffset + UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS - 
          UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS;
      }
#endif
      lmSpikePackets[index] = neuronId;
      lmSpikePackets[index+1] = as_uint(start_gpu);
      break;
    }
    else
    {
      v = V_RT;
      u = U_RT;
      g_ampa = GA_RT;
      g_gaba = GG_RT;
      start_gpu = in_t_now;
    }
    }

    if(spiked < 0)
    {
      if(stp_pass)
      {
        DATA_TYPE w = gm_events[neuronEventPtrStart + 2*UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE];
        uint g_select = (w > 0);
        g_ampa += g_select*w;
        g_gaba -= (!g_select)*w;
      }
      
      neuronEventPtrStart++;
      ev_count--;
    }
  }
  
    /*Store model variables*/
    gm_model_variables[neuronId]                                = v;
    gm_model_variables[UPDATE_NEURONS_TOTAL_NEURONS+neuronId]   = u;
    gm_model_variables[2*UPDATE_NEURONS_TOTAL_NEURONS+neuronId] = g_ampa;
    gm_model_variables[3*UPDATE_NEURONS_TOTAL_NEURONS+neuronId] = g_gaba;
    if(spiked < 0)
    {
      *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1) = 0;
    }
  }

  uint spikePacketSizeWords = (*(lmSpikePacketCounts + LOCAL_WF_ID(wi_id)))* 
    UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS;

  /*Store spike packet*/
#if UPDATE_NEURONS_WF_SIZE_WI >= UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS
  {
    uint i = WI_ID_WF_SCOPE(wi_id);
    if(i < spikePacketSizeWords)
    {
      gm_spikes[UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*GLOBAL_WF_ID(wg_id, wi_id) + i] = 
        lmSpikePackets[UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*LOCAL_WF_ID(wi_id) + i];
    }
  }
#else
  for
  (
    uint i = WI_ID_WF_SCOPE(wi_id); 
    i < spikePacketSizeWords;
    i += UPDATE_NEURONS_WF_SIZE_WI
  ){
    gm_spikes[UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*GLOBAL_WF_ID(wg_id, wi_id) + i] = 
      lmSpikePackets[UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*LOCAL_WF_ID(wi_id) + i];
  }
#endif

  /*Store spike counts*/
  if(WI_ID_WF_SCOPE(wi_id) == 0)
  {
    gm_spike_counts[GLOBAL_WF_ID(wg_id, wi_id)] = lmSpikePacketCounts[LOCAL_WF_ID(wi_id)];
  }
}



/*
  update_spiked_neurons
  
  1) Loads model variables, parameters, and synaptic events from GM to PM.
  2) Iterates through synaptic events and updates model variables.
  3) Detects spiking neurons and computes spike times
  4) Writes spike times and updated model variables back to GM
  
*/
/*TODO: this can be further divided into 2 kernels. 
The last kernel (same as above) can be executed concurently with next propagation kernel.*/
__kernel 
void update_spiked_neurons
(
#if (UPDATE_NEURONS_DEBUG_ENABLE)
  __global      uint            *gm_debug_host,
  __global      uint            *gm_debug_device,
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  __global      uint            *gm_error_code,
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
  __constant    DATA_TYPE       *cm_ps_tol,
#endif
  //__global      DATA_TYPE const* restrict cm_coefficients,
  //constant      DATA_TYPE       *cm_coefficients __attribute__((max_constant_size (4*CONST_SIZE))),
  __constant    DATA_TYPE       *cm_coefficients,
  __global      DATA_TYPE       *gm_model_parameters,
  __global      DATA_TYPE       *gm_model_variables,
  __global      uint            *gm_spike_counts,
  __global      uint            *gm_spikes,
  __global      uint            *gm_event_ptr,
  __global      DATA_TYPE       *gm_events,
                uint            step
){

  uint wi_id = get_local_id(0);
  uint wg_id = get_group_id(0);

  /*Each WF has a spike packet
  __local uint lmSpikePackets[UPDATE_NEURONS_WG_SIZE_WF];
  
  if(WI_ID_WF_SCOPE(wi_id) == 0)
  {
    lmSpikePackets[LOCAL_WF_ID(wi_id)] = 
      gm_spikes[UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*GLOBAL_WF_ID(wg_id, wi_id)];
  }

  barrier(CLK_LOCAL_MEM_FENCE);
  uint  totalSpikes =  lmSpikePackets[LOCAL_WF_ID(wi_id)];*/
  
  /*TODO: verify if this short version is faster than the one above*/
  uint  totalSpikes =  gm_spike_counts[GLOBAL_WF_ID(wg_id, wi_id)];
  
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  if(totalSpikes >= 0xFFFF)
  {
    atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_7);
    totalSpikes = 0;
  }
#endif

  /*A WF works on allocated for it spike packet*/
  for(uint j = WI_ID_WF_SCOPE(wi_id); j < totalSpikes; j += UPDATE_NEURONS_WF_SIZE_WI)
  {
    uint spikeOffset = UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS*GLOBAL_WF_ID(wg_id, wi_id) +
      j*UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS;
    uint  neuronId = gm_spikes[spikeOffset];
    DATA_TYPE start_gpu = as_float(gm_spikes[spikeOffset + 1]);

    /*WI: load event pointer and count, reset the count*/
    uint neuronEventPtrStart = *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE);
    int ev_count = *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1);

#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(ev_count >= 0xFFFF)
    {
      atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_8);
      ev_count = 0;
    }
#endif
    
    /*Load model variables*/
    DATA_TYPE v         = gm_model_variables[neuronId];
    DATA_TYPE u         = gm_model_variables[UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE g_ampa    = gm_model_variables[2*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE g_gaba    = gm_model_variables[3*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    
    /*Load model parameters*/
    DATA_TYPE k         = gm_model_parameters[UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE l         = gm_model_parameters[2*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE b         = gm_model_parameters[3*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE v_reset   = gm_model_parameters[4*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE u_step    = gm_model_parameters[5*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE v_peak    = gm_model_parameters[6*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE E_ampa    = gm_model_parameters[7*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE E_gaba    = gm_model_parameters[8*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
#if UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP > 0
    DATA_TYPE I = 0.0f;
    if(step < UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP)
    {
      I = gm_model_parameters[neuronId];
    }
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
    DATA_TYPE ps_tol = cm_ps_tol[neuronId/(UPDATE_NEURONS_TOTAL_NEURONS/UPDATE_NEURONS_TOLERANCE_CHUNKS)];
#endif
    /*
    DATA_TYPE E         = gm_model_parameters[9*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    DATA_TYPE a         = gm_model_parameters[10*UPDATE_NEURONS_TOTAL_NEURONS+neuronId];
    */
    DATA_TYPE E         = (1.0f/UPDATE_NEURONS_C)*UPDATE_NEURONS_DT; 
    DATA_TYPE a         = UPDATE_NEURONS_a*UPDATE_NEURONS_DT;

    DATA_TYPE start_dt_part, spiked = -1.0f;
    DATA_TYPE cache[UPDT_RUNTM_SZ];

#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  uint spikeCount = 0;
#endif
  
  while(ev_count > 0)
  {
    /*Same as: uint stp_pass = (ev_count > 1);*/
    uint stp_pass = clamp(ev_count-1, 0, 1);
    
    /*Read event time from valid ptr else make it eq 1: */
    DATA_TYPE in_t_now = 1.0f;
    if(stp_pass)
    {
      in_t_now = gm_events[neuronEventPtrStart + UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE];
    }
    
    /*Compute time difference btwn events: */
    if(spiked < 0)
    {
      start_dt_part = in_t_now - start_gpu;
    }

    if(start_dt_part > 0)
    {
    /*Load parameters into work memory space: */
    YP1(0)  = v;
    V_RT    = v;
    U_RT    = u;
    GA_RT   = g_ampa;
    GG_RT   = g_gaba;
    YP2(0)  = MUL(k,v) - g_ampa - g_gaba + l;

    /*numerical integration step function*/
    uint ps_order = gpu_iz_ps_step
    (
#if (UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD)
      v_peak,
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
      gm_error_code,
#endif
#if UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP > 0
      I,
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
      ps_tol,
#endif
      cm_coefficients, 
      cache,
      k,
      E_ampa,
      E_gaba,
      E,
      a,
      b,
      start_dt_part
    );
    
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(ps_order >= UPDATE_NEURONS_PS_ORDER_LIMIT)
    {
      ps_order = UPDATE_NEURONS_PS_ORDER_LIMIT-1;
      atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_1);
    }
#endif

    DATA_TYPE vnew = V_RT;

    spiked = vnew - v_peak;

    /*if a spike occured (rare)*/
    if(spiked >= 0)
    {
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
      spikeCount++;
      if(spikeCount > 1)
      {
        atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_2);
      }
#endif
      
      DATA_TYPE u_tmp,g_ampa_tmp,g_gaba_tmp;
      DATA_TYPE dv,dx,dx_old,dt_part;
      
      /*v-v_peak shifted for root finding*/
      YP1(0) = v - v_peak;
      
      /*First NR step*/
      dt_part = -(DIV(YP1(0),YP1(1)));
      dx_old = 100.0f;
      
      /*NR: */
      int nr_order;
#if (UPDATE_NEURONS_TOLERANCE_MODE == 0)
      DATA_TYPE nr_tol = (DATA_TYPE)UPDATE_NEURONS_NR_ZERO_TOLERANCE;
#elif (UPDATE_NEURONS_TOLERANCE_MODE == 1)
      DATA_TYPE nr_tol = (DATA_TYPE)UPDATE_NEURONS_NR_TOLERANCE;
#elif (UPDATE_NEURONS_TOLERANCE_MODE > 1)
      DATA_TYPE nr_tol = (DATA_TYPE)UPDATE_NEURONS_NR_ZERO_TOLERANCE;
      if(ps_tol != 0){nr_tol = ps_tol;}
#endif
      
      for (nr_order=0; nr_order<UPDATE_NEURONS_NR_ORDER_LIMIT; nr_order++)
      {
        vnew = MUL(YP1(ps_order),dt_part) + YP1(ps_order-1);
        dv = YP1(ps_order);

        for(int j=ps_order-2;j>=0; j--)
        {
          dv = vnew + MUL(dv,dt_part);
          vnew = YP1(j) + MUL(vnew,dt_part);
        }

        /*TODO: divergence may come from division: */
        dx = DIV(vnew,dv);
        dt_part -= dx; 

        if(fabs(dx)<nr_tol)break;
        /*For oscillations*/
        if(fabs(dx+dx_old)<nr_tol)break; 
        dx_old=dx;
      }
      
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
      if(nr_order == UPDATE_NEURONS_NR_ORDER_LIMIT)
      {
        atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_3);
      }
#endif

      /*Handle divergence error: */
      if(dt_part>start_dt_part || dt_part<0)
      {
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
        atomic_or(gm_error_code,UPDATE_NEURONS_ERROR_CODE_4);
#endif
        dt_part = start_dt_part/2;
      }

      /*Record spike and schedule events:*/
      gm_spikes[spikeOffset+1] = as_uint(start_gpu+dt_part);

      /*Evaluate u, g_ampa, g_gaba at corrected spike time: */
      YP1(0) = v;
      YP2(0) = u;
      gpu_ps_update(cm_coefficients, cache,1,ps_order,dt_part,&u_tmp,a,b);
      YP2(0) = g_ampa;
      gpu_ps_update(cm_coefficients, cache,2,ps_order,dt_part,&g_ampa_tmp,a,b);
      YP2(0) = g_gaba;
      gpu_ps_update(cm_coefficients, cache,3,ps_order,dt_part,&g_gaba_tmp,a,b); 
      
      start_dt_part = start_dt_part-dt_part;
      u_tmp += u_step;

      v = v_reset;
      u = u_tmp;
      g_ampa = g_ampa_tmp;
      g_gaba = g_gaba_tmp;
    }
    else
    {
      v = V_RT;
      u = U_RT;
      g_ampa = GA_RT;
      g_gaba = GG_RT;
      start_gpu = in_t_now;
    }
    }

    if(spiked < 0)
    {
      if(stp_pass)
      {
        DATA_TYPE w = gm_events[neuronEventPtrStart + 2*UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE];
        uint g_select = (w > 0);
        g_ampa += g_select*w;
        g_gaba -= (!g_select)*w;
      }
      
      neuronEventPtrStart++;
      ev_count--;
    }
  }
  
    /*Store model variables*/
    gm_model_variables[neuronId]                                = v;
    gm_model_variables[UPDATE_NEURONS_TOTAL_NEURONS+neuronId]   = u;
    gm_model_variables[2*UPDATE_NEURONS_TOTAL_NEURONS+neuronId] = g_ampa;
    gm_model_variables[3*UPDATE_NEURONS_TOTAL_NEURONS+neuronId] = g_gaba;
    *(gm_event_ptr + neuronId*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1) = 0;
  }
}
