
/* ===============================================================================================



  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Operator_Expand.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



#if ENABLE_OPERATOR_EXPAND



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



#if EXPAND_EVENTS_DEBUG_ENABLE
void 
Operator_Expand::invalidateDebug
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_EXPAND_VALID_DEBUG ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ERROR_TRACK_ENABLE
void 
Operator_Expand::invalidateError
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_EXPAND_VALID_ERROR ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ENABLE
void 
Operator_Expand::expand
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt,
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  Data_SpikeEvents                    &spikeEvents,
  Data_SynapticEvents                 &synapticEvents,
  Data_Connectome                     &connectome
)
/**************************************************************************************************/
{
  cl_uint relativeTimeStep = currentTimeStep % timeSlotCount;

#if (EXPAND_EVENTS_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_EXPAND_VALID_DEBUG);
#endif

  if(this->setKernelArguments)
  {
    this->setKernelArguments = false;
    
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelExpandEvents, this->dataExpandEventsDebugHostBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, this->dataExpandEventsDebugDeviceBuffer, 
      this->argNumExpandEvents++);
#endif

#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelExpandEvents, this->dataExpandEventsErrorBuffer, 
      this->argNumExpandEvents++);
#endif

#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
    SET_KERNEL_ARG(this->kernelExpandEvents, synapticEvents.dataHistogramBuffer, 
      this->argNumExpandEvents++);
#endif

    SET_KERNEL_ARG(this->kernelExpandEvents, spikeEvents.dataSpikePacketCountsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, spikeEvents.dataSpikePacketsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, synapticEvents.dataUnsortedEventCountsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, synapticEvents.dataUnsortedEventTargetsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, synapticEvents.dataUnsortedEventDelaysBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, synapticEvents.dataUnsortedEventWeightsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, connectome.dataSynapseTargetsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, connectome.dataSynapseDelaysBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, connectome.dataSynapseWeightsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG(this->kernelExpandEvents, connectome.dataSynapsePointerBuffer, 
      this->argNumExpandEvents++);
  }
  
  SET_KERNEL_ARG(this->kernelExpandEvents, relativeTimeStep, this->argNumExpandEvents);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelExpandEvents, 
    *(this->globalThreadsExpandEvents), *(this->localThreadsExpandEvents), NULL, ndrEvt, 
    "kernelExpandEvents");

#if DEVICE_HOST_DATA_COHERENCE
  synapticEvents.invalidateEvents();
#endif

#if (EXPAND_EVENTS_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_EXPAND_VALID_DEBUG);
#endif

#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Expand::expand:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_EXPAND_VALID_ERROR);}, 
    dataExpandEventsError
  );
#endif
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ENABLE
void 
Operator_Expand::expand
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt,
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  cl_uint                             overWriteSpikesUntilStep,
  double                              spikeBufferMinPercentFill,
  double                              spikeBufferMaxPercentFill,
  double                              spikeNeuronsPermil,
  Data_SpikeEvents                    &spikeEvents,
  Data_SynapticEvents                 &synapticEvents,
  Data_Connectome                     &connectome
)
/**************************************************************************************************/
{
  if(currentTimeStep < overWriteSpikesUntilStep)
  {
    spikeEvents.setEvents
    (
      queue,
      CL_TRUE,
      spikeBufferMinPercentFill,
      spikeBufferMaxPercentFill,
      spikeNeuronsPermil
    );

    if(currentTimeStep == overWriteSpikesUntilStep-1)
    {
      std::cout << "\nCompleted injecting spikes at step " << currentTimeStep << std::endl;
    }
  }

  this->expand
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
    queue,
    ndrEvt,
    currentTimeStep,
    timeSlotCount,
    spikeEvents,
    synapticEvents,
    connectome
  );
}
/**************************************************************************************************/
#endif



#if ENABLE_UNIT_TEST_EXPAND_EVENTS
void 
Operator_Expand::expandUnitTest
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt,
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  cl_uint                             resetTimeSlot,
  cl_int                              testMode,
  double                              spikeBufferMinPercentFill,
  double                              spikeBufferMaxPercentFill,
  double                              spikeNeuronsPermil,
  double                              gabaPercent,
  double                              minDelay,
  double                              maxDelay,
  Data_SpikeEvents                    &spikeEvents,
  Data_SynapticEvents                 &synapticEvents,
  Data_Connectome                     &connectome
)
/**************************************************************************************************/
{
  /*Reset the data with new values every resetAtSlot steps for unit test for better represenation*/
  cl_uint relativeTimeStep = currentTimeStep % timeSlotCount;
  relativeTimeStep = relativeTimeStep % resetTimeSlot;
  bool expandEventsReset = (!relativeTimeStep);
  
  if(expandEventsReset)
  {
    connectome.resetConnections
    (
      queue,
      CL_TRUE,
      gabaPercent,
      minDelay,
      maxDelay
    );
  
    spikeEvents.setEvents
    (
      queue,
      CL_TRUE,
      spikeBufferMinPercentFill,
      spikeBufferMaxPercentFill,
      spikeNeuronsPermil
    );

    synapticEvents.clearUnsortedEvents
    (
      queue, 
      CL_FALSE,
      0x3
    );
  }
  else
  {
    if(testMode < 0)
    {
      if(relativeTimeStep == 1)
      {
        spikeEvents.setEvents
        (
          queue,
          CL_TRUE,
          0.0, 
          0.0, 
          NULL
        );
      }
      else if(relativeTimeStep == resetTimeSlot - 1)
      {
        spikeEvents.setEvents
        (
          queue,
          CL_TRUE,
          100.0, 
          100.0, 
          NULL
        );
      }
      else
      {
        spikeEvents.setEvents
        (
          queue,
          CL_TRUE,
          0.0, 
          (double)(100*relativeTimeStep/resetTimeSlot), 
          NULL
        );
      }
    }
    else if((testMode >= 0) && (testMode <= 1000))
    {
      spikeEvents.setEvents
      (
        queue,
        CL_TRUE,
        0.0, 
        100.0, 
        (double)testMode
      );
    }
#if CLASS_VALIDATION_ENABLE
    else
    {
      THROW_SIMEX("Operator_Expand::expand: invalid test mode: " << testMode);
    }
#endif
  }

  this->expand
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
    queue,
    ndrEvt,
    currentTimeStep,
    timeSlotCount,
    spikeEvents,
    synapticEvents,
    connectome
  );
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_VERIFY_ENABLE
void 
Operator_Expand::verifyExpand
(
  cl::CommandQueue        &queue,
  bool                    reset,
  cl_uint                 timeStep,
  cl_uint                 spikePacketCountWG,
  cl_uint                 histogramBitShift,
  cl_uint                 histogramBinMask,
  double                  maxDelay,
  double                  minDelay,
  Data_SpikeEvents        &spikeEvents,
  Data_SynapticEvents     &synapticEvents,
  Data_Connectome         &connectome
)
/**************************************************************************************************/
{
  bool result = true;
  std::stringstream ss;
  cl_uint error_event_totals = 0;
  cl_uint* dataUnsortedEventCountsVerify = NULL;
  cl_uint* dataUnsortedEventsTargetsVerify = NULL;
  cl_uint* dataUnsortedEventsDelaysVerify = NULL;
  cl_uint* dataUnsortedEventsWeightsVerify = NULL;
  cl_uint* dataHistogramVerify = NULL;
  
  try
  {
  /* get data required synapticEvents parameters*/
  cl_uint eventBufferCount = synapticEvents.getEventBufferCount();
  cl_uint eventTimeSlots = synapticEvents.getEventTimeSlots();
  cl_uint eventBufferSize = synapticEvents.getEventBufferSize();
  cl_uint eventHistogramBinCount = synapticEvents.getEventHistogramBinCount();
  
  /* allocate memory for synaptic events for verification*/
  size_t size = eventBufferCount * eventTimeSlots;
  dataUnsortedEventCountsVerify = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  /* allocate memory for synaptic events for verification*/
  size = eventBufferCount * eventTimeSlots * eventBufferSize;
  dataUnsortedEventsTargetsVerify = (cl_uint *)calloc(size, sizeof(cl_uint));
  dataUnsortedEventsDelaysVerify = (cl_uint *)calloc(size, sizeof(cl_uint));
  dataUnsortedEventsWeightsVerify = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  /* allocate memory for histogram*/
  size = eventTimeSlots*(eventHistogramBinCount*eventBufferCount + 1);
  dataHistogramVerify = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  if(!reset)
  {
    cl_uint targetNeuron = 0, eventTime = 0, weight = 0, eventCount = 0;
    for(cl_uint s = 0; s < eventTimeSlots; s++)
    {
      for(cl_uint b = 0; b < eventBufferCount; b++)
      {
        eventCount = synapticEvents.getEventCount(queue, b, s, Data_SynapticEvents::PREVIOUS);
        dataUnsortedEventCountsVerify[eventTimeSlots*b + s] = eventCount;

        for(cl_uint j = 0; j < eventCount; j++)
        {
          synapticEvents.getEvent(queue, b, j, s, Data_SynapticEvents::PREVIOUS, 
            targetNeuron, eventTime, weight);
            
          cl_uint ptr = 
            /*Event data buffers*/
            b * eventTimeSlots * eventBufferSize +
            /*Current event data buffer*/
            s * eventBufferSize +
            /*Current event*/
            j;
            
          dataUnsortedEventsTargetsVerify[ptr] = targetNeuron;
          dataUnsortedEventsDelaysVerify[ptr] = eventTime;
          dataUnsortedEventsWeightsVerify[ptr] = weight; 
        }
      }
    }
    
    memcpy(dataHistogramVerify, synapticEvents.dataHistogram, 
      synapticEvents.dataHistogramSizeBytes);
  }

  /*Init verification data*/
  cl_uint maxSynapticBuffer = 0;
  cl_uint timeSlot = timeStep%eventTimeSlots;

  /*Reset histogram data in previous step*/
  for
  (
    cl_uint i = 0; 
    i < (eventHistogramBinCount*eventBufferCount+1); 
    i++
  ){
    cl_uint offset = 
    /*time slot*/
    (((eventTimeSlots + (timeSlot-1))%eventTimeSlots)*
    (eventHistogramBinCount*eventBufferCount+1) + 
    /*data*/
    i);
    dataHistogramVerify[offset] = 0;
  }
  
  /*Reset event counter in previous step*/
  for(cl_uint b = 0; b < eventBufferCount; b++)
  {
    cl_uint offset = 
      /*Buffer*/
      eventTimeSlots * b + 
      /*Time slot*/
      ((eventTimeSlots + (timeSlot-1))%eventTimeSlots);
      
    dataUnsortedEventCountsVerify[offset] = 0;
  }
  
  /*Iterate through spike packets and expand packets into events*/
  for(cl_uint packet = 0; packet < spikeEvents.getSpikePacketCount(); packet++)
  {
    cl_uint total_spikes = spikeEvents.getSpikeCount(queue, packet);

    /*Iterate through spikes in a current packet*/
    for(cl_uint i = 0; i < total_spikes; i++)
    {
      cl_uint spiked_neuron = 0;
      cl_float spike_time = 0;
      spikeEvents.getSpike(queue, packet, i, spiked_neuron, spike_time);

      cl_uint synapse_count = connectome.getSynapseCount(spiked_neuron);
      
      /*Iterate through synapses of spiked neuron*/
      for(cl_uint j = 0; j < synapse_count; j++)
      {
        cl_uint buffer = (packet/spikePacketCountWG);
        cl_uint target_neuron = 0; cl_float weight = 0; cl_float delay = 0;
        
        connectome.getSynapse(spiked_neuron, j, target_neuron, delay, weight);
        
        /*Add delay to spike time, decrement*/
        cl_float event_time = spike_time + delay; event_time = event_time - 1.0f;
        /*Make it relative to its time slot*/
        cl_float event_time_binned = event_time - (int)event_time;
        /*Events with 0.0 time bounce back to the previous time slot*/
        if(event_time_binned == 0.0f){event_time_binned = 1.0f;}
        int bin_correction = (int)event_time_binned;
        /*Calculate time slot*/
        cl_uint time_slot = 
          ((timeSlot + (int)event_time - bin_correction)%eventTimeSlots);
        
        /*Obtain offsets for storing synaptic event*/
        cl_uint event_local_ptr = dataUnsortedEventCountsVerify[eventTimeSlots * buffer + 
          time_slot];
        
        /*Catch overflows and record max overflow*/
        if(event_local_ptr >= eventBufferSize)
        {
          if(maxSynapticBuffer < event_local_ptr)
          {
            maxSynapticBuffer = event_local_ptr;
          }
          event_local_ptr = eventBufferSize-1;
        }

        /*Increment event counter*/
        dataUnsortedEventCountsVerify[eventTimeSlots * buffer + time_slot]++;
        
        cl_uint event_global_ptr = 
          /*Event data buffers*/
          buffer * eventTimeSlots * eventBufferSize +
          /*Current event data buffer*/
          time_slot * eventBufferSize +
          /*Current event*/
          event_local_ptr;

        /*Store the event*/
        dataUnsortedEventsTargetsVerify[event_global_ptr] = target_neuron; 
        *((cl_float *)(&dataUnsortedEventsDelaysVerify[event_global_ptr])) = event_time_binned; 
        *((cl_float *)(&dataUnsortedEventsWeightsVerify[event_global_ptr])) = weight;
        
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
        /*Compute histogram key for target neuron based on MSBs*/
        cl_uint event_time_uint = *((cl_uint *)(&event_time_binned));
        cl_uint bin = (event_time_uint >> histogramBitShift) & histogramBinMask;
          
        /*Offset is based on time slot, bin, WG*/
        cl_uint offset = 
          /*WG offset*/
          buffer + 
          /*time slot + bin with eventBufferCount as a pitch*/
          (time_slot*(eventHistogramBinCount*eventBufferCount+1) + 
          bin*eventBufferCount);
          
        /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
        dataHistogramVerify[offset]++;
#endif
      }
    }
  }
  
  /*Fail if overflow*/
  if(maxSynapticBuffer != 0)
  {
    /*Correct time bin counters in case if simulation still proceeds*/
    for(cl_uint packet = 0; packet < spikeEvents.getSpikePacketCount(); packet++)
    {
      for(cl_uint time_slot = 0; time_slot < eventTimeSlots; time_slot++)
      {
        cl_uint buffer = (packet/spikePacketCountWG);

        if(dataUnsortedEventCountsVerify[eventTimeSlots * buffer + time_slot] >= 
          eventBufferSize)
        {
          dataUnsortedEventCountsVerify[eventTimeSlots * buffer + time_slot] = 
            eventBufferSize-1;
        }
      }
    }
    
    ss << "Operator_Expand::verifyExpand: Event buffer overflow. " << 
      "Need to increase eventBufferSize above " 
      << maxSynapticBuffer << std::endl;
    result = false;
  }

  if(result)
  {
    for(cl_uint s = 0; s < eventTimeSlots; s++)
    {
      /*Iterate through synaptic counters*/
      for(cl_uint i = 0; i < eventBufferCount; i++)
      {
        cl_uint eventCount = synapticEvents.getEventCount(queue, i, s, Data_SynapticEvents::RECENT);
        
        if(dataUnsortedEventCountsVerify[eventTimeSlots * i + s] != eventCount)
        {
          error_event_totals++;
          /**/
          ss << i <<"->(" << dataUnsortedEventCountsVerify[eventTimeSlots * i + 
            s] << "," << eventCount << ");";
        }
      }
      
      if(error_event_totals)
      {
        ss << "\n" 
          << "Operator_Expand::verifyExpand: Failed to match synaptic event time slot counters " 
          << error_event_totals << " times." << std::endl;
        result = false; break;
      }
    }
  }

  if(result)
  {
    for(cl_uint s = 0; s < eventTimeSlots; s++)
    {
      /*Iterate through events*/
      for(cl_uint p = 0; p < eventBufferCount; p++)
      {
        cl_uint checksum_event_target_neuron_host = 0, checksum_event_target_neuron_device = 0;
        cl_uint checksum_event_delay_host = 0, checksum_event_delay_device = 0;
        cl_uint checksum_event_weight_host = 0, checksum_event_weight_device = 0;
        cl_uint event_total = dataUnsortedEventCountsVerify[eventTimeSlots*p + s];
        
        /*Iterate through synapses of spiked neuron*/
        for(cl_uint j = 0; j < event_total; j++)
        {
          cl_uint event_global_ptr = 
            /*Event data buffers*/
            p * eventTimeSlots * 
            eventBufferSize +
            /*Current event data buffer*/
            s * eventBufferSize +
            /*Current event*/
            j;
            
          cl_uint targetNeuron = 0, eventTime = 0, weight = 0;
          synapticEvents.getEvent(queue, p, j, s, Data_SynapticEvents::RECENT, 
            targetNeuron, eventTime, weight);
            
          /*Compute checksums*/
          //checksum_event_target_neuron_host ^= dataUnsortedEventsTargetsVerify[event_global_ptr];
          CHECKSUM01(checksum_event_target_neuron_host, 
            dataUnsortedEventsTargetsVerify[event_global_ptr]);
          //checksum_event_target_neuron_device ^= targetNeuron;
          CHECKSUM01(checksum_event_target_neuron_device, targetNeuron);

          //checksum_event_delay_host ^= dataUnsortedEventsDelaysVerify[event_global_ptr];
          CHECKSUM01(checksum_event_delay_host, dataUnsortedEventsDelaysVerify[event_global_ptr]);
          //checksum_event_delay_device ^= eventTime;
          CHECKSUM01(checksum_event_delay_device, eventTime);
          cl_float timeA = *((cl_float *)(&eventTime));
          cl_float timeV = *((cl_float *)(&dataUnsortedEventsDelaysVerify[event_global_ptr]));
          if(timeA < 0.0f || timeA > (float)(maxDelay-minDelay))
          {
            ss << "Operator_Expand::verifyExpand: Found an event in actual data with time outside "
              << "of valid range: value " << timeA << ", range " << 0.0f
              << " - " << (float)(maxDelay-minDelay) << std::endl;
            result = false; break;
          }
          if(timeV< 0.0f || timeV > (float)(maxDelay-minDelay))
          {
            ss << "Operator_Expand::verifyExpand: Found an event in verification data with delay "
              << "outside of valid range: value " << timeV << ", range " << 0.0f
              << " - " << (float)(maxDelay-minDelay) << std::endl;
            result = false; break;
          }
          
          //checksum_event_weight_host ^= dataUnsortedEventsWeightsVerify[event_global_ptr];
          CHECKSUM01(checksum_event_weight_host, dataUnsortedEventsWeightsVerify[event_global_ptr]);
          //checksum_event_weight_device ^= weight;
          CHECKSUM01(checksum_event_weight_device, weight);
        }
        
        if(result == false){break;}
        
        if(
          (checksum_event_target_neuron_host != checksum_event_target_neuron_device) ||
          (checksum_event_delay_host != checksum_event_delay_device) ||
          (checksum_event_weight_host != checksum_event_weight_device)
        ){
          ss << "Operator_Expand::verifyExpand: Failed to verify sunaptic data checksum(s), time slot " 
          << s << ": " << std::endl;
          ss << "\tTarget neuron: host = " << checksum_event_target_neuron_host 
          << " vs device = " << 
            checksum_event_target_neuron_device << std::endl;
          ss << "\tDelay: host = " << checksum_event_delay_host << " vs device = " << 
            checksum_event_delay_device << std::endl;
          ss << "\tWeight host = " << checksum_event_weight_host << " vs device = " << 
            checksum_event_weight_device << std::endl;
            
          result = false; break;
        }
      }
      if(result == false){break;}
    }
  }

#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  if(result)
  {
    for(cl_uint s = 0; s < eventTimeSlots; s++)
    {
      /* init scan and verify data */
      for(cl_uint j = 0; j < eventHistogramBinCount; j++)
      {
        for(cl_uint k = 0; k < eventBufferCount; k++)
        {
          /*Offset is based on time slot, bin, WG*/
          cl_uint offset = 
          /*Time slot space = (time slot #) x (bins per WG) x WGs*/
          s*(eventHistogramBinCount*eventBufferCount+1) + 
          /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
          j*eventBufferCount + k;
          /*Incr counter for a target neuron failing into slot/bin/WG specified by the offset.*/
          error_event_totals += 
            (synapticEvents.getHistogramItem(queue, s, j, k, Data_SynapticEvents::RECENT) != 
            dataHistogramVerify[offset]);
        }
      }
      
      if(error_event_totals)
      {
        ss << "Operator_Expand::verifyExpand: Failed to match histogram " << error_event_totals 
          << " times." << std::endl;
        result = false; break;
      }
    }
  }
#endif
  }
  CATCH(ss, Operator_Expand::verifyExpand, result = false;)
  
  if(dataUnsortedEventCountsVerify)
    free(dataUnsortedEventCountsVerify);
  if(dataUnsortedEventsTargetsVerify)
    free(dataUnsortedEventsTargetsVerify);
  if(dataUnsortedEventsDelaysVerify)
    free(dataUnsortedEventsDelaysVerify);
  if(dataUnsortedEventsWeightsVerify)
    free(dataUnsortedEventsWeightsVerify);
  if(dataHistogramVerify)
    free(dataHistogramVerify);

  if(!result)
  {
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/
#endif



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Operator_Expand::initialize
(
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  size_t                              debugBufferSizeWords,
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  size_t                              errorBufferSizeWords,
#endif
#if (EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  cl::CommandQueue                    &queue,
  cl_bool                             block,
#endif
  cl::Context                         &context,
  cl::Device                          &device,
  struct kernelStatistics             &kernelStats,
  size_t                              cacheSizeWords
)
/**************************************************************************************************/
{
#ifdef WIN32
  QueryPerformanceFrequency((LARGE_INTEGER *)&(this->performanceFrequency));
#endif

#if EXPAND_EVENTS_ENABLE
  /* register device local memory buffer for stats */
  size_t lmExpandEvents = sizeof(cl_uint)*cacheSizeWords; 
  size_t lmExpandEventsSizeBytes = lmExpandEvents;
  REGISTER_MEMORY_O(device, EXPAND_EVENTS_KERNEL_NAME, MEM_LOCAL, lmExpandEvents, kernelStats);
  
#if EXPAND_EVENTS_DEBUG_ENABLE
  /* allocate memory for debug host buffer */
  CALLOC(dataExpandEventsDebugHost, cl_uint, debugBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataExpandEventsDebugHost, kernelStats);
  /* allocate memory for debug device buffer */
  CALLOC(dataExpandEventsDebugDevice, cl_uint, debugBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataExpandEventsDebugDevice, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataExpandEventsDebugHostBuffer, 
    this->dataExpandEventsDebugHostSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataExpandEventsDebugDeviceBuffer, 
    this->dataExpandEventsDebugDeviceSizeBytes);
    
  this->storeData(queue, block, OPERATOR_EXPAND_VALID_DEBUG);
#endif

#if EXPAND_EVENTS_ERROR_TRACK_ENABLE
  /* allocate memory for debug host buffer */
  CALLOC(dataExpandEventsError, cl_uint, errorBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataExpandEventsError, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataExpandEventsErrorBuffer, 
    this->dataExpandEventsErrorSizeBytes);
    
  this->storeData(queue, block, OPERATOR_EXPAND_VALID_ERROR);
#endif

  createKernel
  (
#if LOG_SIMULATION
    this->dataToSimulationLogFile,
#endif
    context,
    device,
    this->kernelExpandEvents,
    EXPAND_EVENTS_KERNEL_FILE_NAME,
    EXPAND_EVENTS_KERNEL_NAME,
    "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF),
    this->blockSizeX_KernelExpandEvents,
    this->blockSizeY_KernelExpandEvents
  );
#endif
}
/**************************************************************************************************/



#if EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE
void
Operator_Expand::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
#if EXPAND_EVENTS_DEBUG_ENABLE
  IF_HIT_READ(selectBitMask, OPERATOR_EXPAND_VALID_DEBUG, this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataExpandEventsDebugHostBuffer, 
      this->dataExpandEventsDebugHostSizeBytes, this->dataExpandEventsDebugHost);
      
    ENQUEUE_READ_BUFFER(block, queue, this->dataExpandEventsDebugDeviceBuffer, 
      this->dataExpandEventsDebugDeviceSizeBytes, this->dataExpandEventsDebugDevice);
      
    this->dataValid |= OPERATOR_EXPAND_VALID_DEBUG;
  }
#endif

#if EXPAND_EVENTS_ERROR_TRACK_ENABLE
  IF_HIT_READ(selectBitMask, OPERATOR_EXPAND_VALID_ERROR, this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataExpandEventsErrorBuffer, 
      this->dataExpandEventsErrorSizeBytes, this->dataExpandEventsError);
      
    this->dataValid |= OPERATOR_EXPAND_VALID_ERROR;
  }
#endif
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE
void
Operator_Expand::storeData
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           selectBitMask
)
/**************************************************************************************************/
{
#if EXPAND_EVENTS_DEBUG_ENABLE
  if(selectBitMask & OPERATOR_EXPAND_VALID_DEBUG)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataExpandEventsDebugHostBuffer, 
      this->dataExpandEventsDebugHostSizeBytes, this->dataExpandEventsDebugHost);
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataExpandEventsDebugDeviceBuffer, 
      this->dataExpandEventsDebugDeviceSizeBytes, this->dataExpandEventsDebugDevice);
    
    this->dataValid |= OPERATOR_EXPAND_VALID_DEBUG;
  }
#endif

#if EXPAND_EVENTS_ERROR_TRACK_ENABLE
  if(selectBitMask & OPERATOR_EXPAND_VALID_ERROR)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataExpandEventsErrorBuffer, 
      this->dataExpandEventsErrorSizeBytes, this->dataExpandEventsError);
    
    this->dataValid |= OPERATOR_EXPAND_VALID_ERROR;
  }
#endif
}
/**************************************************************************************************/
#endif



#endif /*ENABLE_OPERATOR_EXPAND*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
