
/* ===============================================================================================



  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Operator_Group.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



#if ENABLE_OPERATOR_GROUP



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



#if GROUP_EVENTS_DEBUG_ENABLE
void 
Operator_Group::invalidateDebug
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_GROUP_VALID_DEBUG ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ERROR_TRACK_ENABLE
void 
Operator_Group::invalidateError
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_GROUP_VALID_ERROR ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V00
void 
Operator_Group::group_v00
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   currentTimeStep,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
  cl_uint relativeTimeStep = currentTimeStep % (this->timeSlots);

#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v00)
  {
    this->setKernelArguments_v00 = false;
    
#if (GROUP_EVENTS_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV00, this->dataGroupEventsDebugHostBuffer, 
      this->argNumGroupEventsV00++);
    SET_KERNEL_ARG(this->kernelGroupEventsV00, this->dataGroupEventsDebugDeviceBuffer, 
      this->argNumGroupEventsV00++);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV00, this->dataGroupEventsErrorBuffer, 
      this->argNumGroupEventsV00++);
#endif

    SET_KERNEL_ARG(this->kernelGroupEventsV00, synapticEvents.dataUnsortedEventCountsBuffer, 
      this->argNumGroupEventsV00++);
  }
  
  SET_KERNEL_ARG(this->kernelGroupEventsV00, synapticEvents.dataSortedEventsHistogramWriteBuffer, 
      this->argNumGroupEventsV00);
  SET_KERNEL_ARG(this->kernelGroupEventsV00, synapticEvents.dataSortedEventsWriteBuffer, 
    this->argNumGroupEventsV00 + 1);
  SET_KERNEL_ARG(this->kernelGroupEventsV00, synapticEvents.dataUnsortedEventsHistogramBuffer, 
    this->argNumGroupEventsV00 + 2);
  SET_KERNEL_ARG(this->kernelGroupEventsV00, synapticEvents.dataUnsortedEventDelaysBuffer, 
      this->argNumGroupEventsV00 + 3);
  SET_KERNEL_ARG(this->kernelGroupEventsV00, relativeTimeStep, this->argNumGroupEventsV00 + 4);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelGroupEventsV00, 
    *(this->globalThreadsGroupEventsV00), *(this->localThreadsGroupEventsV00), NULL, ndrEvt, 
    "kernelGroupEventsV00");

  swap2(synapticEvents.dataSortedEventsReadBuffer, 
    synapticEvents.dataSortedEventsWriteBuffer);
  swap2(synapticEvents.dataSortedEventsHistogramReadBuffer, 
    synapticEvents.dataSortedEventsHistogramWriteBuffer);

#if DEVICE_HOST_DATA_COHERENCE
    synapticEvents.invalidateSortedEvents();
    synapticEvents.invalidateSortedHistogram();
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Group::group_v00:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_GROUP_VALID_ERROR);}, 
    this->dataGroupEventsError
  );
#endif
}
/**************************************************************************************************/
#endif



#if ENABLE_UNIT_TEST_GROUP_EVENTS_V00
void 
Operator_Group::groupUnitTest_v00
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   currentTimeStep,
  int                       testMode,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
  cl_uint reset = !(currentTimeStep % 5);
  
  if(reset)
  {
    /*Set histogram with non-zero since kernel should zero it out*/
    synapticEvents.setSortedHistogram(queue, CL_TRUE, 0x3, 0xF);
    
    synapticEvents.clearUnsortedEvents(queue, CL_TRUE, 0x2);
    
    double perecentInh = 0, deltaDev = 0;
    int detla = 0;
    
    if(testMode < 0)
    {
      perecentInh = 5.0;
      deltaDev = 50.0; 
      detla = int((((this->eventDataSrcBufMaxSize)*deltaDev)/100)/this->eventDataSrcBufCount);
    }
    else if((testMode >= 0) && (testMode <= 100))
    {
      perecentInh = 5.0;
      deltaDev = (double)testMode; 
      detla = -1;
    }
#if CLASS_VALIDATION_ENABLE
    else
    {
      THROW_SIMEX("Operator_Group::groupUnitTest_v00: invalid test mode: " << testMode);
    }
#endif

    synapticEvents.setUnsortedEvents
    (
#if CLASS_VALIDATION_ENABLE
      (this->eventDataDstBufMaxSize),
#endif
      queue,
      CL_TRUE,
      true,
      detla,
      this->histogramBinShift_v00,
      this->histogramBinMask,
      this->neuronCount,
      deltaDev,
      perecentInh,
      this->minimumDelay,
      0.0
    );
  }
  
  this->group_v00
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
    queue,
    ndrEvt,
    currentTimeStep,
    synapticEvents
  );
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V00
void 
Operator_Group::verifyGroup_v00
(
  cl::CommandQueue        &queue,
  bool                    print,
  cl_uint                 currentTimeSlot,
  Data_SynapticEvents     &synapticEvents
)
/**************************************************************************************************/
{
  bool result = 1;
  std::stringstream *ss = new std::stringstream("");
  std::stringstream *ss_log = new std::stringstream("");

  /*Allocate data for verification*/
  cl_uint *eventPointersCopy = NULL;
  
  synapticEvents.getHistogram
  (
    queue, 
    &eventPointersCopy,
#if CLASS_VALIDATION_ENABLE
    (this->timeSlots)*((this->histogramBinCount)*(this->histogramBinSize) + 1),
#endif
    Data_SynapticEvents::RECENT
  );

  size_t outputHistogramReferenceSize = (this->gridSizeWg)*(this->histogramOutGridSize)*
    (this->histogramOutBinCount);
  cl_uint *outputHistogramReference = (cl_uint *)calloc(outputHistogramReferenceSize, 
    sizeof(cl_uint));
    
  size_t outputEventsReferenceSize = (this->eventDataDstBufMaxSize)*(this->eventDataUnitSize);
  cl_uint *outputEventsReference = (cl_uint *)calloc(outputEventsReferenceSize, sizeof(cl_uint));

  cl_uint total_synaptic_events = synapticEvents.getEventCount(queue, currentTimeSlot, 
    Data_SynapticEvents::RECENT);
  
  print ? *ss_log << "Total synaptic events: " << total_synaptic_events << std::endl,true : false;
  
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = total_synaptic_events/((this->wgSize)*(this->wiDataCount));
  
  if(total_synaptic_event_chunks*((this->wgSize)*(this->wiDataCount)) < total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = total_synaptic_event_chunks/(this->histogramOutGridSize);
  
  if(wg_chunk_size*(this->histogramOutGridSize)  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  
  wg_chunk_size *= ((this->wgSize)*(this->wiDataCount));

  /*Group events for verification*/
  for(cl_uint b = 0; b < this->eventDataSrcBufCount; b++)
  {
    cl_uint event_total = synapticEvents.getEventCount(queue, b, currentTimeSlot, 
      Data_SynapticEvents::RECENT);

    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Access event*/
      cl_uint targetNeuron = 0, key = 0, weight = 0;
      synapticEvents.getEvent(queue, b, e, currentTimeSlot, Data_SynapticEvents::RECENT, 
        targetNeuron, key, weight);
      
      /*Compute offset key for target neuron*/
      cl_uint bin = (key >> (this->histogramBinShift_v00))&
        (this->histogramBinMask);
        
      /*Offset is based on time slot, bin, WG*/
      cl_uint bin_offset = 
        /*Time slot*/
        currentTimeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
        /*WG offset*/
        b +
        /*time slot + bin with histogramBinSize as a pitch*/
        bin*(this->histogramBinSize);

      /*Check for offset overlap*/
      if(eventPointersCopy[bin_offset] >= 
        synapticEvents.getCurrentUnsortedHistogramItem(queue, currentTimeSlot, bin, b+1))
      {
        *ss << "Operator_Group::verifyGroup_v00: Destination event bin pointer overlaps "
          << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot 
          << ", buffer " << b << std::endl;

        result = 0; 
        break;
      }
      
      /*Calculate offset in the grouped data space*/
      cl_uint dest_offset = eventPointersCopy[bin_offset];
        
      /*Compute histogram key for target neuron*/
      size_t hist_out_ptr = 
        /*WG offset*/
        (b/(this->wgEventBufCount)) * 
        ((this->histogramOutGridSize)*(this->histogramOutBinCount)) + 
        /*WG offset for histogram out*/
        (dest_offset/wg_chunk_size) + 
        /*bin*/
        (this->histogramOutGridSize)*
        ((key>>(this->histogramOutBinShift_v00))&
        (this->histogramOutBinMask));
        
      /*Verify*/
      if(hist_out_ptr > outputHistogramReferenceSize)
      {
        *ss << "Operator_Group::verifyGroup_v00: Pointer to an element in output histogram "
          << "is outside of its range" << std::endl;
        result = 0;
        break;
      }
      /*Increment counter for this bin.*/
      outputHistogramReference[hist_out_ptr]++;
      
      /*Store event at its group location (grouped by bins)*/
      outputEventsReference[dest_offset] = key;

      if(this->valuesMode_v00)
      {
        /*Compute pointer to event data*/
        cl_uint ptr = 
          /*Event data buffers*/
          b * (this->timeSlots) * this->eventDataSrcBufMaxSize +
          /*Current event data buffer*/
          currentTimeSlot * this->eventDataSrcBufMaxSize +
          /*Current event*/
          e;

        outputEventsReference[(this->eventDataDstBufMaxSize)+dest_offset] = ptr;
      }

      /*Increment ptr for next data item*/
      eventPointersCopy[bin_offset]++;
    }
    if(!result){break;}
  }
  
  free(eventPointersCopy);
  
  print ? *ss_log << "Time slot " << currentTimeSlot << ": " << std::endl, true : false;

  /*Verify event data*/
  if(result)
  {
    result = this->verifyGroup
    (
      print,
      false,
      (this->sourceEventsDataStructureType_v00 == 0),
      (this->valuesMode_v00 > 0),
      (this->relocateValues_v00),
      true,
      (this->enableStepShift_v00),
      currentTimeSlot,
      0,
      (this->eventDataDstBufMaxSize),
      (this->histogramBinSize),
      (this->histogramBinCount),
      (this->histogramBinMask),
      (this->histogramBinShift_v00),
      (this->histogramOutBinCount),
      (this->histogramOutGridSize),
      (this->gridSizeWg),
      outputHistogramReference,
      outputEventsReference,
      ss,
      ss_log,
      queue,
      synapticEvents
    );
  }

  free(outputHistogramReference);
  free(outputEventsReference);
  
  if(print){LOG_SIM((*ss_log).str());}
  
  if(!result)
  {
    throw SimException((*ss).str());
  }
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V01
void 
Operator_Group::group_v01
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   sortStep,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v01)
  {
    this->setKernelArguments_v01 = false;
    
#if (GROUP_EVENTS_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV01, this->dataGroupEventsDebugHostBuffer, 
      this->argNumGroupEventsV01++);
    SET_KERNEL_ARG(this->kernelGroupEventsV01, this->dataGroupEventsDebugDeviceBuffer, 
      this->argNumGroupEventsV01++);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV01, this->dataGroupEventsErrorBuffer, 
      this->argNumGroupEventsV01++);
#endif
  }
  
  SET_KERNEL_ARG(this->kernelGroupEventsV01, synapticEvents.dataSortedEventsHistogramWriteBuffer, 
    this->argNumGroupEventsV01);
  SET_KERNEL_ARG(this->kernelGroupEventsV01, synapticEvents.dataSortedEventsWriteBuffer, 
    this->argNumGroupEventsV01 + 1);
  SET_KERNEL_ARG(this->kernelGroupEventsV01, synapticEvents.dataSortedEventsHistogramReadBuffer, 
    this->argNumGroupEventsV01 + 2);
  SET_KERNEL_ARG(this->kernelGroupEventsV01, synapticEvents.dataSortedEventsReadBuffer, 
    this->argNumGroupEventsV01 + 3);
  SET_KERNEL_ARG(this->kernelGroupEventsV01, sortStep, 
    this->argNumGroupEventsV01 + 4);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelGroupEventsV01, 
    *(this->globalThreadsGroupEventsV01), *(this->localThreadsGroupEventsV01), NULL, ndrEvt, 
    "kernelGroupEventsV01");

  swap2(synapticEvents.dataSortedEventsReadBuffer, 
    synapticEvents.dataSortedEventsWriteBuffer);
  swap2(synapticEvents.dataSortedEventsHistogramReadBuffer, 
    synapticEvents.dataSortedEventsHistogramWriteBuffer);

#if DEVICE_HOST_DATA_COHERENCE
    synapticEvents.invalidateSortedEvents();
    synapticEvents.invalidateSortedHistogram();
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Group::group_v01:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_GROUP_VALID_ERROR);}, 
    this->dataGroupEventsError
  );
#endif
}
/**************************************************************************************************/
#endif



#if ENABLE_UNIT_TEST_GROUP_EVENTS_V01
void 
Operator_Group::groupUnitTest_v01
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   sortStep,
  double                    testMode,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
  int step = ((int)sortStep)-1;
  cl_uint keyOffset = 1, shiftFirstStage = 0, shiftNextStage = 0;

  if(testMode < 0)
  {
    testMode = 0;
  }
  else if(testMode > 100)
  {
    throw SimException("Operator_Group::groupUnitTest_v01: Incorrect test mode");
  }

  if(step >= 0)
  {
    shiftFirstStage = (this->histogramBinShift_v01) * (cl_uint)step;
    shiftNextStage = (this->histogramOutBinShift_v01) * ((cl_uint)step + 1);
  }

  synapticEvents.initializeGrouppedEvents
  (
    false,
    true,
    (this->neuronCount),
    keyOffset,
    (this->histogramBinMask),
    (this->histogramOutBinMask),
    shiftFirstStage,
    shiftNextStage,
    (this->wgSize),
    (this->wiDataCount),
    10.0,
    (this->minimumDelay),
    testMode,
    queue
  );

  this->group_v01
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    currentTimeStep,
#endif
    queue,
    ndrEvt,
    sortStep,
    synapticEvents
  );
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V01
void 
Operator_Group::verifyGroup_v01
(
  cl::CommandQueue        &queue,
  cl_uint                 currentTimeSlot,
  cl_uint                 sortStep,
  Data_SynapticEvents     &synapticEvents
)
/**************************************************************************************************/
{
  this->setReferenceAndVerifyGroup
  (
    false,
    true,
    (this->relocateValues_v01),
    false,
    (this->sourceEventsDataStructureType_v01 == 0),
    true,
    this->enableStepShift_v01,
    -1,
    currentTimeSlot,
    sortStep,
    this->valuesMode_v01,
    this->histogramBinShift_v01,
    this->histogramOutBinShift_v01,
    synapticEvents,
    queue
  );
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V02
void 
Operator_Group::group_v02
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   sortStep,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v02)
  {
    this->setKernelArguments_v02 = false;
    
#if (GROUP_EVENTS_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV02, this->dataGroupEventsDebugHostBuffer, 
      this->argNumGroupEventsV02++);
    SET_KERNEL_ARG(this->kernelGroupEventsV02, this->dataGroupEventsDebugDeviceBuffer, 
      this->argNumGroupEventsV02++);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV02, this->dataGroupEventsErrorBuffer, 
      this->argNumGroupEventsV02++);
#endif

    SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataUnsortedEventTargetsBuffer, 
      this->argNumGroupEventsV02++);
    SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataUnsortedEventDelaysBuffer, 
      this->argNumGroupEventsV02++);
    SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataUnsortedEventWeightsBuffer, 
      this->argNumGroupEventsV02++);
  }
  
  SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataSortedEventsHistogramWriteBuffer, 
    this->argNumGroupEventsV02);
  SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataSortedEventsWriteBuffer, 
    this->argNumGroupEventsV02 + 1);
  SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataSortedEventsHistogramReadBuffer, 
    this->argNumGroupEventsV02 + 2);
  SET_KERNEL_ARG(this->kernelGroupEventsV02, synapticEvents.dataSortedEventsReadBuffer, 
    this->argNumGroupEventsV02 + 3);
  SET_KERNEL_ARG(this->kernelGroupEventsV02, sortStep, 
    this->argNumGroupEventsV02 + 4);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelGroupEventsV02, 
    *(this->globalThreadsGroupEventsV02), *(this->localThreadsGroupEventsV02), NULL, ndrEvt, 
    "kernelGroupEventsV02");

  swap2(synapticEvents.dataSortedEventsReadBuffer, 
    synapticEvents.dataSortedEventsWriteBuffer);
  swap2(synapticEvents.dataSortedEventsHistogramReadBuffer, 
    synapticEvents.dataSortedEventsHistogramWriteBuffer);

#if DEVICE_HOST_DATA_COHERENCE
    synapticEvents.invalidateSortedEvents();
    synapticEvents.invalidateSortedHistogram();
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Group::group_v02:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_GROUP_VALID_ERROR);}, 
    this->dataGroupEventsError
  );
#endif
}
/**************************************************************************************************/
#endif



#if ENABLE_UNIT_TEST_GROUP_EVENTS_V02
void 
Operator_Group::groupUnitTest_v02
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   sortStep,
  double                    testMode,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
  int step = ((int)sortStep)-1;
  cl_uint keyOffset = 1, shiftFirstStage = 0, shiftNextStage = 0;

  if(testMode < 0)
  {
    testMode = 0;
  }
  else if(testMode > 100)
  {
    throw SimException("Operator_Group::groupUnitTest_v02: Incorrect test mode");
  }

  if(step >= 0)
  {
    shiftFirstStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS * (cl_uint)step;
    shiftNextStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT * ((cl_uint)step + 1);
  }

  synapticEvents.initializeGrouppedEvents
  (
    false,
    false,
    (this->neuronCount),
    keyOffset,
    (this->histogramBinMask),
    (this->histogramOutBinMask),
    shiftFirstStage,
    shiftNextStage,
    (this->wgSize),
    (this->wiDataCount),
    10.0,
    (this->minimumDelay),
    testMode,
    queue
  );

  this->group_v02
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    currentTimeStep,
#endif
    queue,
    ndrEvt,
    sortStep,
    synapticEvents
  );
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V02
void 
Operator_Group::verifyGroup_v02
(
  cl::CommandQueue        &queue,
  cl_uint                 currentTimeSlot,
  cl_uint                 sortStep,
  Data_SynapticEvents     &synapticEvents
)
/**************************************************************************************************/
{
  this->setReferenceAndVerifyGroup
  (
    false,
    true,
    (this->relocateValues_v02),
    false,
    (this->sourceEventsDataStructureType_v02 == 0),
    false,
    this->enableStepShift_v02,
    0,
    currentTimeSlot,
    sortStep,
    this->valuesMode_v02,
    this->histogramBinShift_v02,
    this->histogramOutBinShift_v02,
    synapticEvents,
    queue
  );
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V03
void 
Operator_Group::group_v03
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   sortStep,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v03)
  {
    this->setKernelArguments_v03 = false;
    
#if (GROUP_EVENTS_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV03, this->dataGroupEventsDebugHostBuffer, 
      this->argNumGroupEventsV03++);
    SET_KERNEL_ARG(this->kernelGroupEventsV03, this->dataGroupEventsDebugDeviceBuffer, 
      this->argNumGroupEventsV03++);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelGroupEventsV03, this->dataGroupEventsErrorBuffer, 
      this->argNumGroupEventsV03++);
#endif

    SET_KERNEL_ARG(this->kernelGroupEventsV03, synapticEvents.dataUnsortedEventTargetsBuffer, 
      this->argNumGroupEventsV03++);
    SET_KERNEL_ARG(this->kernelGroupEventsV03, synapticEvents.dataUnsortedEventDelaysBuffer, 
      this->argNumGroupEventsV03++);
    SET_KERNEL_ARG(this->kernelGroupEventsV03, synapticEvents.dataUnsortedEventWeightsBuffer, 
      this->argNumGroupEventsV03++);
  }

  SET_KERNEL_ARG(this->kernelGroupEventsV03, synapticEvents.dataSortedEventsWriteBuffer, 
    this->argNumGroupEventsV03);
  SET_KERNEL_ARG(this->kernelGroupEventsV03, synapticEvents.dataSortedEventsHistogramReadBuffer, 
    this->argNumGroupEventsV03 + 1);
  SET_KERNEL_ARG(this->kernelGroupEventsV03, synapticEvents.dataSortedEventsReadBuffer, 
    this->argNumGroupEventsV03 + 2);
  SET_KERNEL_ARG(this->kernelGroupEventsV03, sortStep, 
    this->argNumGroupEventsV03 + 3);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelGroupEventsV03, 
    *(this->globalThreadsGroupEventsV03), *(this->localThreadsGroupEventsV03), NULL, ndrEvt, 
    "kernelGroupEventsV03");

  swap2(synapticEvents.dataSortedEventsReadBuffer, 
    synapticEvents.dataSortedEventsWriteBuffer);

#if DEVICE_HOST_DATA_COHERENCE
  synapticEvents.invalidateSortedEvents();
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_GROUP_EVENTS_VALID_DEBUG);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Group::group_v03:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_GROUP_VALID_ERROR);}, 
    this->dataGroupEventsError
  );
#endif
}
/**************************************************************************************************/
#endif



#if ENABLE_UNIT_TEST_GROUP_EVENTS_V03
void 
Operator_Group::groupUnitTest_v03
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt,
  cl_uint                   sortStep,
  double                    testMode,
  Data_SynapticEvents       &synapticEvents
)
/**************************************************************************************************/
{
  int step = ((int)sortStep)-1;
  cl_uint keyOffset = 0, shiftFirstStage = 0, shiftNextStage = 0;

  if(testMode < 0)
  {
    testMode = 0;
  }
  else if(testMode > 100)
  {
    throw SimException("Operator_Group::groupUnitTest_v02: Incorrect test mode");
  }

  if(step >= 0)
  {
    shiftFirstStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS * (cl_uint)step;
    shiftNextStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT * ((cl_uint)step + 1);
  }

  synapticEvents.initializeGrouppedEvents
  (
    false,
    false,
    (this->neuronCount),
    keyOffset,
    (this->histogramBinMask),
    (this->histogramOutBinMask),
    shiftFirstStage,
    shiftNextStage,
    (this->wgSize),
    (this->wiDataCount),
    10.0,
    (this->minimumDelay),
    testMode,
    queue
  );

  this->group_v03
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    currentTimeStep,
#endif
    queue,
    ndrEvt,
    sortStep,
    synapticEvents
  );
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V03
void 
Operator_Group::verifyGroup_v03
(
  cl::CommandQueue        &queue,
  cl_uint                 currentTimeSlot,
  cl_uint                 sortStep,
  Data_SynapticEvents     &synapticEvents
)
/**************************************************************************************************/
{
  this->setReferenceAndVerifyGroup
  (
    false,
    false,
    (this->relocateValues_v03),
    true,
    (this->sourceEventsDataStructureType_v03 == 0),
    true,
    this->enableStepShift_v03,
    0,
    currentTimeSlot,
    sortStep,
    this->valuesMode_v03,
    this->histogramBinShift_v03,
    this->histogramOutBinShift_v03,
    synapticEvents,
    queue
  );
}
/**************************************************************************************************/
#endif



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Operator_Group::initialize
(
#if (GROUP_EVENTS_DEBUG_ENABLE)
  size_t                              debugBufferSizeWords,
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  size_t                              errorBufferSizeWords,
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
  cl::CommandQueue                    &queue,
  cl_bool                             block,
#endif
#if GROUP_EVENTS_ENABLE_V00
  size_t                              cacheSizeWordsV00,
#endif
#if GROUP_EVENTS_ENABLE_V01
  size_t                              cacheSizeWordsV01,
#endif
#if GROUP_EVENTS_ENABLE_V02
  size_t                              cacheSizeWordsV02,
#endif
#if GROUP_EVENTS_ENABLE_V03
  size_t                              cacheSizeWordsV03,
#endif
  cl::Context                         &context,
  cl::Device                          &device,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
  QueryPerformanceFrequency((LARGE_INTEGER *)&(this->performanceFrequency));
#endif
  
#if GROUP_EVENTS_DEBUG_ENABLE
  /* allocate memory for debug host buffer */
  CALLOC(dataGroupEventsDebugHost, cl_uint, debugBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataGroupEventsDebugHost, kernelStats);
  /* allocate memory for debug device buffer */
  CALLOC(dataGroupEventsDebugDevice, cl_uint, debugBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataGroupEventsDebugDevice, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataGroupEventsDebugHostBuffer, 
    this->dataGroupEventsDebugHostSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataGroupEventsDebugDeviceBuffer, 
    this->dataGroupEventsDebugDeviceSizeBytes);
    
  this->storeData(queue, block, OPERATOR_GROUP_VALID_DEBUG);
#endif

#if GROUP_EVENTS_ERROR_TRACK_ENABLE
  /* allocate memory for debug host buffer */
  CALLOC(dataGroupEventsError, cl_uint, errorBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataGroupEventsError, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataGroupEventsErrorBuffer, 
    this->dataGroupEventsErrorSizeBytes);
    
  this->storeData(queue, block, OPERATOR_GROUP_VALID_ERROR);
#endif

#if GROUP_EVENTS_ENABLE_V00
  {
    /* register device local memory buffer for stats */
    std::string kernelTag = std::string(GROUP_EVENTS_KERNEL_NAME "_V00");
    size_t lmGroupEvents = sizeof(cl_uint)*cacheSizeWordsV00; 
    size_t lmGroupEventsSizeBytes = lmGroupEvents;
    REGISTER_MEMORY_O(device, kernelTag, MEM_LOCAL, lmGroupEvents, kernelStats);

    createKernel
    (
#if LOG_SIMULATION
      this->dataToSimulationLogFile,
#endif
      context,
      device,
      this->kernelGroupEventsV00,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V00=1",
      this->blockSizeX_kernelGroupEventsV00,
      this->blockSizeY_kernelGroupEventsV00
    );
  }
#endif

#if GROUP_EVENTS_ENABLE_V01
  {
    /* register device local memory buffer for stats */
    std::string kernelTag = std::string(GROUP_EVENTS_KERNEL_NAME "_V01");
    size_t lmGroupEvents = sizeof(cl_uint)*cacheSizeWordsV01; 
    size_t lmGroupEventsSizeBytes = lmGroupEvents;
    REGISTER_MEMORY_O(device, kernelTag, MEM_LOCAL, lmGroupEvents, kernelStats);

    createKernel
    (
#if LOG_SIMULATION
      this->dataToSimulationLogFile,
#endif
      context,
      device,
      this->kernelGroupEventsV01,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V01=1",
      this->blockSizeX_kernelGroupEventsV01,
      this->blockSizeY_kernelGroupEventsV01
    );
  }
#endif

#if GROUP_EVENTS_ENABLE_V02
  {
    /* register device local memory buffer for stats */
    std::string kernelTag = std::string(GROUP_EVENTS_KERNEL_NAME "_V02");
    size_t lmGroupEvents = sizeof(cl_uint)*cacheSizeWordsV02; 
    size_t lmGroupEventsSizeBytes = lmGroupEvents;
    REGISTER_MEMORY_O(device, kernelTag, MEM_LOCAL, lmGroupEvents, kernelStats);

    createKernel
    (
#if LOG_SIMULATION
      this->dataToSimulationLogFile,
#endif
      context,
      device,
      this->kernelGroupEventsV02,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V02=1",
      this->blockSizeX_kernelGroupEventsV02,
      this->blockSizeY_kernelGroupEventsV02
    );
  }
#endif

#if GROUP_EVENTS_ENABLE_V03
  {
    /* register device local memory buffer for stats */
    std::string kernelTag = std::string(GROUP_EVENTS_KERNEL_NAME "_V03");
    size_t lmGroupEvents = sizeof(cl_uint)*cacheSizeWordsV03; 
    size_t lmGroupEventsSizeBytes = lmGroupEvents;
    REGISTER_MEMORY_O(device, kernelTag, MEM_LOCAL, lmGroupEvents, kernelStats);

    createKernel
    (
#if LOG_SIMULATION
      this->dataToSimulationLogFile,
#endif
      context,
      device,
      this->kernelGroupEventsV03,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V03=1",
      this->blockSizeX_kernelGroupEventsV03,
      this->blockSizeY_kernelGroupEventsV03
    );
  }
#endif
}
/**************************************************************************************************/



#if GROUP_EVENTS_VERIFY_ENABLE
void 
Operator_Group::setReferenceAndVerifyGroup
(
  bool                    print,
  bool                    enableOutputHistogram,
  bool                    relocateValues,
  bool                    useCurrentHistogram,
  bool                    useUnsortedHistogram,
  bool                    verifyKeyBins,
  bool                    enableStepShift,
  int                     replacementKeyOffset,
  cl_uint                 currentTimeSlot,
  cl_uint                 sortStep,
  cl_uint                 valuesMode,
  cl_uint                 histogramBinShift,
  cl_uint                 histogramOutBinShift,
  Data_SynapticEvents     &synapticEvents,
  cl::CommandQueue        &queue
)
/**************************************************************************************************/
{
  bool result = 1;
  std::stringstream *ss = new std::stringstream("");
  std::stringstream *ss_log = new std::stringstream("");

  /* allocate memory for offset copies used in verification data generation*/
  cl_uint *eventPointersCopy = NULL;
  
  if(useCurrentHistogram)
  {
    synapticEvents.getCurrentSortedHistogram
    (
#if CLASS_VALIDATION_ENABLE
      (this->gridSizeWg)*(this->histogramOutGridSize)*(this->histogramOutBinCount),
#endif
      queue, 
      &eventPointersCopy
    );
  }
  else
  {
    synapticEvents.getPreviousSortedHistogram
    (
#if CLASS_VALIDATION_ENABLE
      (this->gridSizeWg)*(this->histogramOutGridSize)*(this->histogramOutBinCount),
#endif
      queue, 
      &eventPointersCopy
    );
  }

  /*Allocate data for verification*/
  cl_uint *outputHistogramReference = NULL;
  
  size_t outputHistogramReferenceSize = (this->gridSizeWg)*(this->histogramOutGridSize)*
    (this->histogramOutBinCount);

  if(enableOutputHistogram)
  {
    outputHistogramReference = (cl_uint *)calloc(outputHistogramReferenceSize, 
      sizeof(cl_uint));
  }
  
  size_t outputEventsReferenceSize = (this->eventDataDstBufMaxSize)*(this->eventDataUnitSize);
  cl_uint *outputEventsReference = (cl_uint *)calloc(outputEventsReferenceSize, sizeof(cl_uint));
  
  cl_uint event_total = 0;
  
  if(useCurrentHistogram)
  {
    event_total = synapticEvents.getCurrentSortedHistogramItem(queue, 
      (this->histogramBinCount)*(this->histogramBinSize));
  }
  else
  {
    event_total = synapticEvents.getPreviousSortedHistogramItem(queue, 
      (this->histogramBinCount)*(this->histogramBinSize));
  }
  
  print ? *ss_log << "Total synaptic events: " << event_total << std::endl,true : false;
    
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = event_total/((this->wgSize)*(this->wiDataCount));
    
  if(total_synaptic_event_chunks*((this->wgSize)*(this->wiDataCount)) < event_total)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks for WG*/
  cl_uint total_wg_chunks = total_synaptic_event_chunks/(this->gridSizeWg);
  
  if(total_wg_chunks*(this->gridSizeWg)  < total_synaptic_event_chunks)
  {
    total_wg_chunks++;
  }
  
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = 1;
  
  if(enableOutputHistogram)
  {
    wg_chunk_size = total_synaptic_event_chunks/(this->histogramOutGridSize);
    
    if(wg_chunk_size*(this->histogramOutGridSize)  < total_synaptic_event_chunks)
    {
      wg_chunk_size++;
    }
    
    wg_chunk_size *= ((this->wgSize)*(this->wiDataCount));
  }
  
  bool verifyValue0 = false, verifyValue1 = false;

  /*Group events for verification*/
  for(cl_uint e = 0; e < event_total; e++)
  {
    /*Access event*/
    cl_uint key = synapticEvents.getEvent(queue, e, Data_SynapticEvents::PREVIOUS);
    cl_uint value = synapticEvents.getEvent(queue, (this->eventDataDstBufMaxSize) + e, 
      Data_SynapticEvents::PREVIOUS);
    
    cl_uint new_key = key, time = 0, weight = 0;
    
    if(replacementKeyOffset >= 0)
    {
      synapticEvents.getEvent(queue, value+replacementKeyOffset, Data_SynapticEvents::RECENT, 
        new_key, time, weight);
    }

    /*Compute WG which is working on this element*/
    cl_uint wg_id = e/(total_wg_chunks*((this->wgSize)*(this->wiDataCount)));

    /*Compute bin for neuron*/
    cl_uint bin = (key>>histogramBinShift)&(this->histogramBinMask);
    
    if(enableStepShift)
    {
      bin = (key>>(histogramBinShift*sortStep))&(this->histogramBinMask);
    }

    /*Offset is based on time slot, bin, WG*/
    cl_uint bin_offset = 
    /*WG offset*/
    wg_id +
    /*time slot + bin with (this->histogramBinSize) as a pitch*/
    bin*(this->histogramBinSize);

    /*Check for offset overlap*/
    if(useCurrentHistogram)
    {
      if
      (
        eventPointersCopy[bin_offset] >= 
        synapticEvents.getCurrentSortedHistogramItem(queue, bin_offset+1)
      )
      {
        *ss << "Operator_Group::setReferenceAndVerifyGroup: Destination event bin pointer "
          "overlaps with next pointer for bin " << bin << ", time slot " << currentTimeSlot 
          << ", WG " << wg_id << ", pointer " << e << ", sort step " << sortStep << std::endl;
        result = 0; 
        break;
      }

    }
    else
    {
      if
      (
        eventPointersCopy[bin_offset] >= 
        synapticEvents.getPreviousSortedHistogramItem(queue, bin_offset+1)
      )
      {
        *ss << "Operator_Group::setReferenceAndVerifyGroup: Destination event bin pointer "
          "overlaps with next pointer for bin " << bin << ", time slot " << currentTimeSlot 
          << ", WG " << wg_id << ", pointer " << e << ", sort step " << sortStep << std::endl;
        result = 0; 
        break;
      }
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = eventPointersCopy[bin_offset];

    if(enableOutputHistogram)
    {
      /*Compute histogram key for target neuron*/
      size_t hist_out_ptr = 
        /*WG offset*/
        (wg_id/(this->wgEventBufCount)) * 
        ((this->histogramOutGridSize)*(this->histogramOutBinCount)) + 
        /*WG offset for histogram out*/
        dest_offset/wg_chunk_size;
        
      /*bin*/
      if(enableStepShift)
      {
        hist_out_ptr += (this->histogramOutGridSize)*
        ((new_key>>(histogramOutBinShift*(sortStep+1)))&(this->histogramOutBinMask));
      }
      else
      {
        hist_out_ptr += (this->histogramOutGridSize)*
        ((new_key>>histogramOutBinShift)&(this->histogramOutBinMask));
      }
        
      /*Verify overflow*/
      if(hist_out_ptr > outputHistogramReferenceSize)
      {
        *ss << "Operator_Group::setReferenceAndVerifyGroup: Pointer to an element in output "
          "histogram is outside of its range" << std::endl;
        result = 0;
        break;
      }
      /*Increment counter for this bin.*/
      outputHistogramReference[hist_out_ptr]++;
    }
    
    /*Store event at its group location (grouped by bins)*/
    outputEventsReference[dest_offset] = new_key;

    if(valuesMode == 1)
    {
      verifyValue0 = true;
      outputEventsReference[(this->eventDataDstBufMaxSize) + dest_offset] = e;
    }
    else if(valuesMode == 2)
    {
      verifyValue0 = true;
      
      if(relocateValues)
      {
        verifyValue1 = true;
        outputEventsReference[(this->eventDataDstBufMaxSize) + dest_offset] = time;
        outputEventsReference[2*(this->eventDataDstBufMaxSize) + dest_offset] = weight;
      }
      else
      {
        outputEventsReference[(this->eventDataDstBufMaxSize) + dest_offset] = value;
      }
    }
    
    /*Increment ptr for next data item*/
    eventPointersCopy[bin_offset]++;
  }
  
  free(eventPointersCopy);
  
  print ? *ss_log << "Time slot " << currentTimeSlot << ": " << std::endl, true : false;

  if(result)
  {
    result = this->verifyGroup
    (
      print,
      useCurrentHistogram,
      useUnsortedHistogram,
      verifyValue0,
      verifyValue1,
      verifyKeyBins,
      enableStepShift,
      currentTimeSlot,
      sortStep,
      (this->eventDataDstBufMaxSize),
      (this->histogramBinSize),
      (this->histogramBinCount),
      (this->histogramBinMask),
      histogramBinShift,
      (this->histogramOutBinCount),
      (this->histogramOutGridSize),
      (this->gridSizeWg),
      outputHistogramReference,
      outputEventsReference,
      ss,
      ss_log,
      queue,
      synapticEvents
    );
  }

  if(outputHistogramReference)
    free(outputHistogramReference);

  free(outputEventsReference);
  
  if(print){LOG_SIM((*ss_log).str());}
  
  if(!result)
  {
    throw SimException((*ss).str());
  }
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_VERIFY_ENABLE
bool 
Operator_Group::verifyGroup
(
  bool  print,
  bool  useCurrentHistogram,
  bool  useUnsortedHistogram,
  bool  verifyValue0,
  bool  verifyValue1,
  bool verifyKeyBins,
  bool enableStepShift,
  cl_uint timeSlot,
  cl_uint step,
  cl_uint destinationBufferSize,
  cl_uint histogramBinSize,
  cl_uint histogramTotalBins,
  cl_uint histogramBitMask,
  cl_uint histogramBitShift,
  cl_uint histogramOutTotalBins,
  cl_uint histogramOutTotalGroups,
  cl_uint groupSize,
  cl_uint *outputHistogramReference,
  cl_uint *outputEventsReference,
  std::stringstream     *ss,
  std::stringstream     *ss_log,
  cl::CommandQueue      &queue,
  Data_SynapticEvents   &synapticEvents
)
/**************************************************************************************************/
{
  bool result = 1;
  
  /*Verify event data*/
  for(cl_uint j = 0; j < histogramTotalBins; j++)
  {
    /*Error counters*/
    unsigned long long verify_checksum_target_neuron = 0;
    unsigned long long actual_checksum_target_neuron = 0;
    unsigned long long verify_checksum_value_01 = 0;
    unsigned long long actual_checksum_value_01 = 0;
    unsigned long long verify_checksum_value_02 = 0;
    unsigned long long actual_checksum_value_02 = 0;
    unsigned long long verify_error_count_target_neuron = 0;
    unsigned long long actual_error_count_target_neuron = 0;
    
    /*Start and end bin ptr*/
    cl_uint start = 0, end = 0;

    if(useUnsortedHistogram)
    {
      start = synapticEvents.getCurrentUnsortedHistogramItem(queue, timeSlot, j, 0);
      end = synapticEvents.getCurrentUnsortedHistogramItem(queue, timeSlot, j+1, 0);
    }
    else
    {
      if(useCurrentHistogram)
      {
        start = synapticEvents.getCurrentSortedHistogramItem(queue, j*histogramBinSize);
        end = synapticEvents.getCurrentSortedHistogramItem(queue, (j+1)*histogramBinSize);
      }
      else
      {
        start = synapticEvents.getPreviousSortedHistogramItem(queue, j*histogramBinSize);
        end = synapticEvents.getPreviousSortedHistogramItem(queue, (j+1)*histogramBinSize);
      }
    }

    print ? *ss_log << "Start-End: " << start << "-" << end << std::endl, true : false;
    
    /*Verify correct bin and checksum in that bin for current time slot*/
    for(cl_uint p = start; p < end; p++)
    {
      cl_uint v_key = 0, a_key = 0;
      
      if(enableStepShift)
      {
        v_key = (outputEventsReference[p] >> (histogramBitShift*step)) & histogramBitMask;
        a_key = ((synapticEvents.getEvent(queue, p, Data_SynapticEvents::RECENT)) >> 
          (histogramBitShift*step)) & histogramBitMask;
      }
      else
      {
        v_key = (outputEventsReference[p] >> histogramBitShift) & histogramBitMask;
        a_key = ((synapticEvents.getEvent(queue, p, Data_SynapticEvents::RECENT)) >> 
          histogramBitShift) & histogramBitMask;
      }
      
      if(verifyKeyBins)
      {
        verify_error_count_target_neuron += (v_key != j);
        actual_error_count_target_neuron += (a_key != j);
        
        (print && (v_key != j)) ? *ss_log << "V," << p << ": (" << v_key << " != " << j << ")" 
          << "; ", true : false;
        (print && (a_key != j)) ? *ss_log << "A," << p << ": (" << a_key << " != " << j << ")" 
          << "; ", true : false;
      }
      
      CHECKSUM01(verify_checksum_target_neuron, outputEventsReference[p]);
      /*verify_checksum_target_neuron += outputEventsReference[p]; */
      CHECKSUM01(actual_checksum_target_neuron, 
        (synapticEvents.getEvent(queue, p, Data_SynapticEvents::RECENT)));
      /*actual_checksum_target_neuron += 
        (synapticEvents.getEvent(queue, p, Data_SynapticEvents::RECENT)); */
      
      if(verifyValue0)
      {
        CHECKSUM01(verify_checksum_value_01, outputEventsReference[destinationBufferSize + p]);
        /*verify_checksum_value_01 += outputEventsReference[destinationBufferSize + p]; */
        CHECKSUM01(actual_checksum_value_01, (synapticEvents.getEvent(queue, 
          destinationBufferSize + p, Data_SynapticEvents::RECENT)));
        /*actual_checksum_value_01 += (synapticEvents.getEvent(queue, destinationBufferSize + p, 
          Data_SynapticEvents::RECENT)); */
      }
      
      if(verifyValue1)
      {
        CHECKSUM01(verify_checksum_value_02, outputEventsReference[2*destinationBufferSize + p]);
        /*verify_checksum_value_02 += outputEventsReference[2*destinationBufferSize + p]; */
        CHECKSUM01(actual_checksum_value_02, (synapticEvents.getEvent(queue, 
          2*destinationBufferSize + p, Data_SynapticEvents::RECENT)));
        /*actual_checksum_value_02 += (synapticEvents.getEvent(queue, 
          2*destinationBufferSize + p, Data_SynapticEvents::RECENT)); */
      }
    }
    
    print ? *ss_log << std::endl, true : false;
    
    if(verifyKeyBins)
    {
      if(verify_error_count_target_neuron)
      {
        *ss << "Operator_Group::verifyGroup: Failed to validate correct bin for " 
          << verify_error_count_target_neuron 
          << " keys out of " << (end-start) << " in verification data in bin " << j 
          << ", time slot " << timeSlot << std::endl;
        result = 0;
        break;
      }
      
      if(actual_error_count_target_neuron)
      {
        *ss << "Operator_Group::verifyGroup: Failed to validate correct bin for " 
          << actual_error_count_target_neuron 
          << " keys out of " << (end-start) << " in actual data in bin " << j 
          << ", time slot " << timeSlot << std::endl;
        result = 0;
        break;
      }
    }
    
    if(verify_checksum_target_neuron != actual_checksum_target_neuron)
    {
      *ss << "Operator_Group::verifyGroup: Failed to match neuron checksum in bin " 
        << j << ", time slot " << timeSlot << std::endl;
      result = 0;
      break;
    }
    
    if(verifyValue0 && (verify_checksum_value_01 != actual_checksum_value_01))
    {
      *ss << "Operator_Group::verifyGroup: Failed to match value 0 checksum in bin " 
        << j << ", time slot " 
        << timeSlot << std::endl;
      result = 0;
      break;
    }
    
    if(verifyValue1 && (verify_checksum_value_02 != actual_checksum_value_02))
    {
      *ss << "Operator_Group::verifyGroup: Failed to match value 1 checksum in bin " 
        << j << ", time slot " 
        << timeSlot << std::endl;
      result = 0;
      break;
    }
    
    if(print && useUnsortedHistogram && (!result))
    {
      *ss_log << "WG Pointers for bin " << j << "\nWG,Ptr Start,Ptr End,Index,Key Verified,"
        "Key Actual,F/P" << std::endl;
        
      for(cl_uint w = 0; w < histogramBinSize; w++)
      {
        cl_uint start = synapticEvents.getCurrentUnsortedHistogramItem(queue, timeSlot, j, w);
        cl_uint end = synapticEvents.getCurrentUnsortedHistogramItem(queue, timeSlot, j, w+1) - 1;
          
        for(cl_uint p = start; p < end; p++)
        {
          cl_uint v_key = outputEventsReference[p];
          cl_uint a_key = synapticEvents.getEvent(queue, p, Data_SynapticEvents::RECENT);
          
          *ss_log << w << "," << start << "," << end << "," << p << "," << v_key << "," << a_key 
            << ","; 
          if(v_key != a_key){*ss_log << "F";}
          *ss_log << std::endl;
        }
      }
    }
  }

  if(result)
  {
    print ? *ss_log << std::endl, true : false;
    
    /*Verify event data*/
    if(outputHistogramReference != NULL)
    {
      for(cl_uint w = 0; w < (groupSize); w++)
      {
        for(cl_uint b = 0; b < (histogramOutTotalBins); b++)
        {
          unsigned long long verify_checksum_histogram_out = 0;
          unsigned long long actual_checksum_histogram_out = 0;
          
          for(cl_uint j = 0; j < (histogramOutTotalGroups); j++)
          {
            cl_uint p = 
              /*WG offset*/
              w*histogramOutTotalGroups*histogramOutTotalBins + 
              /*WG offset for histogram out*/
              j +
              /*bin*/
              b*histogramOutTotalGroups;
              
            CHECKSUM01(verify_checksum_histogram_out, outputHistogramReference[p]);
            /*verify_checksum_histogram_out += outputHistogramReference[p]; */
            CHECKSUM01(actual_checksum_histogram_out, 
              synapticEvents.getCurrentSortedHistogramItem(queue, w, b, j));
            /*actual_checksum_histogram_out += 
              synapticEvents.getCurrentSortedHistogramItem(queue, w, b, j); */
          }
          
          /*Verify*/
          if(verify_checksum_histogram_out != actual_checksum_histogram_out)
          {
            *ss << "Operator_Group::verifyGroup: Failed to match output histogram "
              << "checksum in bin " << b << ", WG " << w << ", time slot " << timeSlot 
              << std::endl;
            result = 0;
            
            print ? *ss_log << "Output histogram:" << std::endl, true : false;
            
            for(cl_uint j = 0; j < histogramOutTotalGroups; j++)
            {
              cl_uint p = 
                /*WG offset*/
                w*histogramOutTotalGroups*histogramOutTotalBins + 
                /*WG offset for histogram out*/
                j +
                /*bin*/
                b*histogramOutTotalGroups;
                
              print ? *ss_log << p << "," << outputHistogramReference[p] << ","
                << synapticEvents.getCurrentSortedHistogramItem(queue, w, b, j) 
                << std::endl, true : false;
            }
            
            break;
          }
        }
        if(!result){break;}
      }
    }
  }

  return result;
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
void
Operator_Group::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
#if GROUP_EVENTS_DEBUG_ENABLE
  IF_HIT_READ(selectBitMask, OPERATOR_GROUP_VALID_DEBUG, this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataGroupEventsDebugHostBuffer, 
      this->dataGroupEventsDebugHostSizeBytes, this->dataGroupEventsDebugHost);
      
    ENQUEUE_READ_BUFFER(block, queue, this->dataGroupEventsDebugDeviceBuffer, 
      this->dataGroupEventsDebugDeviceSizeBytes, this->dataGroupEventsDebugDevice);
      
    this->dataValid |= OPERATOR_GROUP_VALID_DEBUG;
  }
#endif

#if GROUP_EVENTS_ERROR_TRACK_ENABLE
  IF_HIT_READ(selectBitMask, OPERATOR_GROUP_VALID_ERROR, this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataGroupEventsErrorBuffer, 
      this->dataGroupEventsErrorSizeBytes, this->dataGroupEventsError);
      
    this->dataValid |= OPERATOR_GROUP_VALID_ERROR;
  }
#endif
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
void
Operator_Group::storeData
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           selectBitMask
)
/**************************************************************************************************/
{
#if GROUP_EVENTS_DEBUG_ENABLE
  if(selectBitMask & OPERATOR_GROUP_VALID_DEBUG)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataGroupEventsDebugHostBuffer, 
      this->dataGroupEventsDebugHostSizeBytes, this->dataGroupEventsDebugHost);
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataGroupEventsDebugDeviceBuffer, 
      this->dataGroupEventsDebugDeviceSizeBytes, this->dataGroupEventsDebugDevice);
    
    this->dataValid |= OPERATOR_GROUP_VALID_DEBUG;
  }
#endif

#if GROUP_EVENTS_ERROR_TRACK_ENABLE
  if(selectBitMask & OPERATOR_GROUP_VALID_ERROR)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataGroupEventsErrorBuffer, 
      this->dataGroupEventsErrorSizeBytes, this->dataGroupEventsError);
    
    this->dataValid |= OPERATOR_GROUP_VALID_ERROR;
  }
#endif
}
/**************************************************************************************************/
#endif



#endif /*ENABLE_OPERATOR_GROUP*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
