
/* ===============================================================================================



  =============================================================================================== */


  
/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Data_SynapticEvents.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



void 
Data_SynapticEvents::clearSortedEvents
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           clearBitMask
)
/**************************************************************************************************/
{
  if(clearBitMask & 0x1)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataSortedEvents, this->dataPastSortedEvents);
    swap2(this->dataSortedEventsHistogram, this->dataPastSortedEventsHistogram);
#endif

    memset(this->dataSortedEvents, 0, this->dataSortedEventsSizeBytes);
    memset(this->dataSortedEventsHistogram, 0, this->dataSortedEventsHistogramSizeBytes);
    
    this->storeData(queue, block, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS | 
      SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM);
  }

#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  if(clearBitMask & 0x2)
  {
    memset(this->dataPastSortedEvents, 0, this->dataSortedEventsSizeBytes);
    memset(this->dataPastSortedEventsHistogram, 0, this->dataSortedEventsHistogramSizeBytes);
  }
#endif
}
/**************************************************************************************************/



void 
Data_SynapticEvents::clearUnsortedEvents
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           clearBitMask
)
/**************************************************************************************************/
{
  if(clearBitMask & 0x1)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataUnsortedEventCounts, this->dataPastUnsortedEventCounts);
    swap2(this->dataUnsortedEventTargets, this->dataPastUnsortedEventTargets);
    swap2(this->dataUnsortedEventDelays, this->dataPastUnsortedEventDelays);
    swap2(this->dataUnsortedEventWeights, this->dataPastUnsortedEventWeights);
    swap2(this->dataUnsortedEventsHistogram, this->dataPastUnsortedEventsHistogram);
#endif

    memset(this->dataUnsortedEventCounts, 0, this->dataUnsortedEventCountsSizeBytes);
    memset(this->dataUnsortedEventTargets, 0, this->dataUnsortedEventTargetsSizeBytes);
    memset(this->dataUnsortedEventDelays, 0, this->dataUnsortedEventDelaysSizeBytes);
    memset(this->dataUnsortedEventWeights, 0, this->dataUnsortedEventWeightsSizeBytes);
    memset(this->dataUnsortedEventsHistogram, 0, this->dataUnsortedEventsHistogramSizeBytes);
    
    this->storeData(queue, block, (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS | 
      SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS | SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM));
  }
  
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  if(clearBitMask & 0x2)
  {
    memset(this->dataPastUnsortedEventCounts, 0, this->dataUnsortedEventCountsSizeBytes);
    memset(this->dataPastUnsortedEventTargets, 0, this->dataUnsortedEventTargetsSizeBytes);
    memset(this->dataPastUnsortedEventDelays, 0, this->dataUnsortedEventDelaysSizeBytes);
    memset(this->dataPastUnsortedEventWeights, 0, this->dataUnsortedEventWeightsSizeBytes);
    memset(this->dataPastUnsortedEventsHistogram, 0, this->dataUnsortedEventsHistogramSizeBytes);
  }
#endif
}
/**************************************************************************************************/



void
Data_SynapticEvents::initializeGrouppedEvents
/**************************************************************************************************/
(
  bool      print,
  bool      useTempUnsortedEvents,
  cl_uint   totalNeurons,
  cl_uint   keyOffset,
  cl_uint   histogramBinMask,
  cl_uint   histogramOutBinMask,
  cl_uint   histogramBitShift,
  cl_uint   histogramOutBitShift,
  cl_uint   wgSize,
  cl_uint   wiDataCount,
  cl_float  percentInhibitory,
  double    minDelay,
  double    percentBufferSizeDeviation,
  cl::CommandQueue  &queue
){
  bool result = 1;
  cl_uint max_offset = 0, runningSum = 0, runningSize = 0;
  std::stringstream ss;

  cl_uint *dataUnsortedEventCounts = NULL;
  cl_uint *dataUnsortedEventTargets = NULL;
  cl_uint *dataUnsortedEventDelays = NULL;
  cl_uint *dataUnsortedEventWeights = NULL;
  cl_uint *dataUnsortedEventsHistogram = NULL;

  if(useTempUnsortedEvents)
  {
    cl_uint size = (this->eventBufferCount);
    dataUnsortedEventCounts = (cl_uint *)calloc(size, sizeof(cl_uint));
    size = (this->eventBufferCount) * (this->eventBufferSize);
    dataUnsortedEventTargets = (cl_uint *)calloc(size, sizeof(cl_uint));
    dataUnsortedEventDelays = (cl_uint *)calloc(size, sizeof(cl_uint));
    dataUnsortedEventWeights = (cl_uint *)calloc(size, sizeof(cl_uint));
    size = ((this->histogramBinCount)*(this->histogramBinSize) + 1);
    dataUnsortedEventsHistogram = (cl_uint *)calloc(size, sizeof(cl_uint));
  }
  else
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataUnsortedEventCounts, this->dataPastUnsortedEventCounts);
    swap2(this->dataUnsortedEventTargets, this->dataPastUnsortedEventTargets);
    swap2(this->dataUnsortedEventDelays, this->dataPastUnsortedEventDelays);
    swap2(this->dataUnsortedEventWeights, this->dataPastUnsortedEventWeights);
    swap2(this->dataUnsortedEventsHistogram, this->dataPastUnsortedEventsHistogram);
#endif

    memset(this->dataUnsortedEventCounts, 0, this->dataUnsortedEventCountsSizeBytes);
    memset(this->dataUnsortedEventTargets, 0, this->dataUnsortedEventTargetsSizeBytes);
    memset(this->dataUnsortedEventDelays, 0, this->dataUnsortedEventDelaysSizeBytes);
    memset(this->dataUnsortedEventWeights, 0, this->dataUnsortedEventWeightsSizeBytes);
    memset(this->dataUnsortedEventsHistogram, 0, this->dataUnsortedEventsHistogramSizeBytes);

    dataUnsortedEventCounts = this->dataUnsortedEventCounts;
    dataUnsortedEventTargets = this->dataUnsortedEventTargets;
    dataUnsortedEventDelays = this->dataUnsortedEventDelays;
    dataUnsortedEventWeights = this->dataUnsortedEventWeights;
    dataUnsortedEventsHistogram = this->dataUnsortedEventsHistogram;
  }
  
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < this->eventBufferCount; b++)
  {
    cl_uint eventCount = cl_uint((this->eventBufferSize)*(1 - percentBufferSizeDeviation/100.0) +
        abs((this->eventBufferSize)*(percentBufferSizeDeviation/100.0)*(double)rand()/
        ((double)RAND_MAX)));

    if(eventCount >= (this->eventBufferSize)){eventCount = (this->eventBufferSize);}

    dataUnsortedEventCounts[b] = eventCount;
    
    cl_uint inhibitorySynapseCount = cl_uint(eventCount*(percentInhibitory/100.0));

    this->initializeEventBuffer
    (
      keyOffset,
      eventCount,
      0,
      b,
      totalNeurons,
      1,
      inhibitorySynapseCount,
      histogramBitShift,
      histogramBinMask,
      minDelay,
      1.0,
      dataUnsortedEventTargets,
      dataUnsortedEventDelays,
      dataUnsortedEventWeights,
      dataUnsortedEventsHistogram
    );
  }

  /*Compute offsets based on histogram.*/
  print ? std::cout 
    << "Data_SynapticEvents::initializeGrouppedEvents: Number of synaptic events in bins: " 
    << std::endl, true : false;

  for(cl_uint j = 0; j < (this->histogramBinCount)*(this->histogramBinSize); j++)
  {
    cl_uint temp = dataUnsortedEventsHistogram[j];
    dataUnsortedEventsHistogram[j] = runningSum;
    runningSum += temp;

    if(j%(this->histogramBinSize) == 0 && j != 0)
    {
      print ? std::cout << runningSum-runningSize << ", ", true : false;
      runningSize = runningSum;
    }
  }
  
  dataUnsortedEventsHistogram[(this->histogramBinCount)*(this->histogramBinSize)] = runningSum;
  
  print ? std::cout << runningSum-runningSize << "] = " << runningSum << 
    std::endl, true : false;
  
  if(max_offset < runningSum){max_offset = runningSum;}

  if(max_offset > (this->sortedEventBufferSize))
  {
    ss << "Data_SynapticEvents::initializeGrouppedEvents: Destination event buffer overflow. " 
      << "Need to increase sorted events buffer size, which is currently " << 
      (this->sortedEventBufferSize) << " above " << max_offset << std::endl;
    result = 0;
  }

  if(result)
  {
    try
    {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
      swap2(this->dataSortedEventsHistogram, this->dataPastSortedEventsHistogram);
#endif

      memset(this->dataSortedEventsHistogram, 0, this->dataSortedEventsHistogramSizeBytes);
    
      this->groupEvents
      (
        queue,
        CL_TRUE,
        1,
        this->dataSortedEventsHistogramSize,
        keyOffset,
        (this->eventBufferCount),
        (this->eventBufferSize),
        histogramBinMask,
        histogramBitShift,
        histogramOutBitShift,
        (this->sortedEventsHistogramBinSize),
        histogramOutBinMask,
        (this->sortedEventsHistogramBinSize),
        wgSize,
        wiDataCount,
        dataUnsortedEventCounts,
        dataUnsortedEventTargets,
        dataUnsortedEventDelays,
        dataUnsortedEventWeights,
        dataUnsortedEventsHistogram,
        this->dataSortedEventsHistogram
      );
    }
    CATCH(ss, Data_SynapticEvents::initializeGrouppedEvents, result = 0;)
  }
  
  if(useTempUnsortedEvents)
  {
    free(dataUnsortedEventsHistogram);
    free(dataUnsortedEventTargets);
    free(dataUnsortedEventDelays);
    free(dataUnsortedEventWeights);
    free(dataUnsortedEventCounts);
  }

  if(result)
  {
    runningSum = this->dataSortedEventsHistogram[0];
    this->dataSortedEventsHistogram[0] = 0;
    
    /*Compute offsets*/
    for(cl_uint j = 1; j < ((this->sortedEventsHistogramBinCount)*
      (this->sortedEventsHistogramBinSize) + 1); j++)
    {
      cl_uint d = this->dataSortedEventsHistogram[j];
      this->dataSortedEventsHistogram[j] = runningSum;
      runningSum += d;
    }
  }

  this->storeData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM);
  
  if(!useTempUnsortedEvents)
  {
    this->storeData(queue, CL_TRUE, (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS | 
      SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS | SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM));
  }

  if(!result)
  {
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/



void 
Data_SynapticEvents::setUnsortedEvents
(
#if CLASS_VALIDATION_ENABLE
  cl_uint             eventDestinationBufferSize,
#endif
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_bool             initHistogram,
  int                 unsortedEventTimeSlotDelta,
  cl_uint             histogramBitShift,
  cl_uint             histogramBinMask,
  cl_uint             neuronCount,
  double              percentTimeSlotDeltaDeviation,
  double              percentInh,
  double              minDelay,
  double              percentDelayBoarder
)
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  swap2(this->dataUnsortedEventCounts, this->dataPastUnsortedEventCounts);
  swap2(this->dataUnsortedEventTargets, this->dataPastUnsortedEventTargets);
  swap2(this->dataUnsortedEventDelays, this->dataPastUnsortedEventDelays);
  swap2(this->dataUnsortedEventWeights, this->dataPastUnsortedEventWeights);
  swap2(this->dataUnsortedEventsHistogram, this->dataPastUnsortedEventsHistogram);
#endif
  
  memset(this->dataUnsortedEventCounts, 0, this->dataUnsortedEventCountsSizeBytes);
  memset(this->dataUnsortedEventTargets, 0, this->dataUnsortedEventTargetsSizeBytes);
  memset(this->dataUnsortedEventDelays, 0, this->dataUnsortedEventDelaysSizeBytes);
  memset(this->dataUnsortedEventWeights, 0, this->dataUnsortedEventWeightsSizeBytes);
  memset(this->dataUnsortedEventsHistogram, 0, this->dataUnsortedEventsHistogramSizeBytes);

  this->initializeEventBuffers
  (
    unsortedEventTimeSlotDelta,
    histogramBitShift,
    histogramBinMask,
    neuronCount,
    minDelay,
    percentDelayBoarder,
    percentTimeSlotDeltaDeviation,
    percentInh
  );

  if(initHistogram)
  {
    this->initializeHistogram
    (
#if CLASS_VALIDATION_ENABLE
      eventDestinationBufferSize
#endif
    );
  }

  this->storeData(queue, block, (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS | 
    SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS | SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM));
}
/**************************************************************************************************/



void 
Data_SynapticEvents::setSortedHistogram
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           clearBitMask,
  int               value
)
/**************************************************************************************************/
{
  if(clearBitMask & 0x1)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataSortedEventsHistogram, this->dataPastSortedEventsHistogram);
#endif

    memset(this->dataSortedEventsHistogram, value, this->dataSortedEventsHistogramSizeBytes);
    
    this->storeData(queue, block, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM);
  }
  
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  if(clearBitMask & 0x2)
  {
    memset(this->dataPastSortedEventsHistogram, value, this->dataSortedEventsHistogramSizeBytes);
  }
#endif
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getEventCount
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(timeSlot >= this->timeSlots)
  {
    throw SimException("Data_SynapticEvents::getEventCount: time slot exceeds total time slots");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM);

  /*Pointer is based on time slot, bin, WG*/
  cl_uint ptr = 
  /*Time slot space = (time slot #) x (bins per WG) x WGs*/
  timeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
  /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
  (this->histogramBinCount)*(this->histogramBinSize);
  
  if(type == Data_SynapticEvents::RECENT)
  {
    return this->dataUnsortedEventsHistogram[ptr];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    return this->dataPastUnsortedEventsHistogram[ptr];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getEventCount: incorrect item type");
  }
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getEventCount
(
  cl::CommandQueue  &queue,
  cl_uint           bufferID,
  cl_uint           timeSlot,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(bufferID >= this->eventBufferCount)
  {
    throw SimException("Data_SynapticEvents::getEventCount: buffer ID exceeds buffer count");
  }
  
  if(timeSlot >= this->timeSlots)
  {
    throw SimException("Data_SynapticEvents::getEventCount: time slot exceeds total time slots");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS);
  
  if(type == Data_SynapticEvents::RECENT)
  {
    return this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    return this->dataPastUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getEventCount: incorrect event count type");
  }
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getEvent
(
  cl::CommandQueue  &queue,
  cl_uint           eventID,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(eventID >= (this->sortedEventBufferCount)*(this->sortedEventBufferSize))
  {
    throw SimException("Data_SynapticEvents::getEvent: event ID exceeds total event count");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS);

  if(type == Data_SynapticEvents::RECENT)
  {
    return this->dataSortedEvents[eventID];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    return this->dataPastSortedEvents[eventID];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getEvent: incorrect item type");
  }
}
/**************************************************************************************************/



void
Data_SynapticEvents::getEvent
(
  cl::CommandQueue  &queue,
  cl_uint           bufferID,
  cl_uint           eventID,
  cl_uint           timeSlot,
  const int         type,
  cl_uint           &targetNeuron,
  cl_uint           &eventTime,
  cl_uint           &weight
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(bufferID >= this->eventBufferCount)
  {
    throw SimException("Data_SynapticEvents::getEvent: buffer ID exceeds buffer count");
  }
  
  if(timeSlot >= this->timeSlots)
  {
    throw SimException("Data_SynapticEvents::getEvent: time slot exceeds total time slots");
  }

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS | 
    SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS);
  
  if(type == Data_SynapticEvents::RECENT)
  {
    if
    (
      (eventID >= this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot]) &&
      (eventID != 0) &&
      (this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot] != 0)
    )
    {
      THROW_SIMEX("Data_SynapticEvents::getEvent: event ID " << eventID 
        << " exceeds RECENT event count " 
        << this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot]);
    }
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    if
    (
      (eventID >= this->dataPastUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot]) &&
      (eventID != 0) &&
      (this->dataPastUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot] != 0)
    )
    {
      THROW_SIMEX("Data_SynapticEvents::getEvent: event ID " << eventID 
        << " exceeds PREVIOUS event count " 
        << this->dataPastUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot]);
    }
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getEvent: incorrect event type");
  }
#endif

  cl_uint ptr = 
    /*Event data buffers*/
    bufferID * (this->timeSlots) * (this->eventBufferSize) +
    /*Current event data buffer*/
    timeSlot * (this->eventBufferSize) +
    /*Current event*/
    eventID;
  
  this->getEvent(queue, ptr, type, targetNeuron, eventTime, weight);
}
/**************************************************************************************************/



void
Data_SynapticEvents::getEvent
(
  cl::CommandQueue  &queue,
  cl_uint           ptr,
  const int         type,
  cl_uint           &targetNeuron,
  cl_uint           &eventTime,
  cl_uint           &weight
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if
  (
    (ptr >= this->dataUnsortedEventTargetsSize) || 
    (ptr >= this->dataUnsortedEventDelaysSize) || 
    (ptr >= this->dataUnsortedEventWeightsSize)
  )
  {
    throw SimException("Data_SynapticEvents::getEvent: pointer exceeds data size");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS);
  
  if(type == Data_SynapticEvents::RECENT)
  {
    targetNeuron = this->dataUnsortedEventTargets[ptr];
    eventTime = this->dataUnsortedEventDelays[ptr];
    weight = this->dataUnsortedEventWeights[ptr];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    targetNeuron = this->dataPastUnsortedEventTargets[ptr];
    eventTime = this->dataPastUnsortedEventDelays[ptr];
    weight = this->dataPastUnsortedEventWeights[ptr];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getEvent: incorrect event type");
  }
}
/**************************************************************************************************/



void
Data_SynapticEvents::getHistogram
(
  cl::CommandQueue  &queue,
  cl_uint           **dataHistogramCopy,
#if CLASS_VALIDATION_ENABLE
  size_t            size,
#endif
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(size > this->dataUnsortedEventsHistogramSize)
  {
    throw SimException("Data_SynapticEvents::getHistogram: insufficient buffer size");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM);

  if(type == Data_SynapticEvents::RECENT)
  {
    if(*dataHistogramCopy)
      free(*dataHistogramCopy);
    
    *dataHistogramCopy = (cl_uint *)calloc(this->dataUnsortedEventsHistogramSize, sizeof(cl_uint));
    memcpy(*dataHistogramCopy, this->dataUnsortedEventsHistogram, 
      this->dataUnsortedEventsHistogramSizeBytes);
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    if(*dataHistogramCopy)
      free(*dataHistogramCopy);
    
    *dataHistogramCopy = (cl_uint *)calloc(this->dataUnsortedEventsHistogramSize, sizeof(cl_uint));
    memcpy(*dataHistogramCopy, this->dataPastUnsortedEventsHistogram, 
      this->dataUnsortedEventsHistogramSizeBytes);
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getHistogram: incorrect item type");
  }
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getPreviousUnsortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Data_SynapticEvents* mySelf = (Data_SynapticEvents*)objectToCallThisFunction;
  return mySelf->getUnsortedHistogramItem(queue, timeSlot, binID, bufferID, 
    Data_SynapticEvents::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getPreviousUnsortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID
)
/**************************************************************************************************/
{
  return this->getUnsortedHistogramItem(queue, timeSlot, binID, bufferID, 
    Data_SynapticEvents::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getCurrentUnsortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Data_SynapticEvents* mySelf = (Data_SynapticEvents*)objectToCallThisFunction;
  return mySelf->getUnsortedHistogramItem(queue, timeSlot, binID, bufferID, 
    Data_SynapticEvents::RECENT);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getCurrentUnsortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID
)
/**************************************************************************************************/
{
  return this->getUnsortedHistogramItem(queue, timeSlot, binID, bufferID, 
    Data_SynapticEvents::RECENT);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getPreviousSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Data_SynapticEvents* mySelf = (Data_SynapticEvents*)objectToCallThisFunction;
  return mySelf->getSortedHistogramItem(queue, backetID, binID, itemID, 
    Data_SynapticEvents::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getPreviousSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID
)
/**************************************************************************************************/
{
  return this->getSortedHistogramItem(queue, backetID, binID, itemID, 
    Data_SynapticEvents::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getPreviousSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           pointer
)
/**************************************************************************************************/
{
  return this->getSortedHistogramItem(queue, pointer, Data_SynapticEvents::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getCurrentSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Data_SynapticEvents* mySelf = (Data_SynapticEvents*)objectToCallThisFunction;
  return mySelf->getSortedHistogramItem(queue, backetID, binID, itemID, 
    Data_SynapticEvents::RECENT);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getCurrentSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID
)
/**************************************************************************************************/
{
  return this->getSortedHistogramItem(queue, backetID, binID, itemID, 
    Data_SynapticEvents::RECENT);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getCurrentSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           pointer
)
/**************************************************************************************************/
{
  return this->getSortedHistogramItem(queue, pointer, Data_SynapticEvents::RECENT);
}
/**************************************************************************************************/



void
Data_SynapticEvents::getPreviousSortedHistogram
(
#if CLASS_VALIDATION_ENABLE
  size_t            size,
#endif
  cl::CommandQueue  &queue,
  cl_uint           **dataHistogramCopy
)
/**************************************************************************************************/
{
  this->getSortedHistogram
  (
    queue, 
    dataHistogramCopy,
#if CLASS_VALIDATION_ENABLE
    size,
#endif
    Data_SynapticEvents::PREVIOUS
  );
}
/**************************************************************************************************/



void
Data_SynapticEvents::getCurrentSortedHistogram
(
#if CLASS_VALIDATION_ENABLE
  size_t            size,
#endif
  cl::CommandQueue  &queue,
  cl_uint           **dataHistogramCopy
)
/**************************************************************************************************/
{
  this->getSortedHistogram
  (
    queue, 
    dataHistogramCopy,
#if CLASS_VALIDATION_ENABLE
    size,
#endif
    Data_SynapticEvents::RECENT
  );
}
/**************************************************************************************************/



void 
Data_SynapticEvents::refreshSortedEventsHistogram
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  this->getData
  (
    queue, 
    block, 
    SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM
  );
}
/**************************************************************************************************/



void 
Data_SynapticEvents::invalidateSortedEvents
()
/**************************************************************************************************/
{
  this->dataValid &= (SYNAPTIC_EVENTS_VALID_SORTED_EVENTS ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
Data_SynapticEvents::invalidateUnsortedEvents
()
/**************************************************************************************************/
{
  this->dataValid &= 
    (
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS  ^ 0xFFFFFFFF) & 
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS        ^ 0xFFFFFFFF) & 
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM         ^ 0xFFFFFFFF)
    );
}
/**************************************************************************************************/



void 
Data_SynapticEvents::invalidateSortedHistogram
()
/**************************************************************************************************/
{
  this->dataValid &= (SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
Data_SynapticEvents::invalidateUnsortedHistogram
()
/**************************************************************************************************/
{
  this->dataValid &= (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
Data_SynapticEvents::refreshUnsortedHistogram
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  this->getData
  (
    queue, 
    block, 
    SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM
  );
}
/**************************************************************************************************/



void 
Data_SynapticEvents::refreshSortedHistogram
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  this->getData
  (
    queue, 
    block, 
    SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM
  );
}
/**************************************************************************************************/



void 
Data_SynapticEvents::refreshSortedEvents
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  this->getData
  (
    queue, 
    block, 
    SYNAPTIC_EVENTS_VALID_SORTED_EVENTS
  );
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventBufferCount
()
/**************************************************************************************************/
{
  return this->eventBufferCount;
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventTimeSlots
()
/**************************************************************************************************/
{
  return this->timeSlots;
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventBufferSize
()
/**************************************************************************************************/
{
  return this->eventBufferSize;
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventHistogramBinCount
()
/**************************************************************************************************/
{
  return this->histogramBinCount;
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Data_SynapticEvents::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
  size_t size = (this->eventBufferCount)*(this->timeSlots);
  
  CALLOC(dataUnsortedEventCounts, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventCounts, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventCounts = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  size = (this->eventBufferCount)*(this->timeSlots)*(this->eventBufferSize);
  
  CALLOC(dataUnsortedEventTargets, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventTargets, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventTargets = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif

  CALLOC(dataUnsortedEventDelays, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventDelays, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventDelays = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  CALLOC(dataUnsortedEventWeights, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventWeights, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventWeights = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif

  size = (this->timeSlots)*((this->histogramBinCount)*(this->histogramBinSize) + 1);
  CALLOC(dataUnsortedEventsHistogram, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventsHistogram, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventsHistogram = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif

  size = (this->sortedEventsHistogramBacketCount)*(this->sortedEventsHistogramBinSize)*
    (this->sortedEventsHistogramBinCount);
  CALLOC(dataSortedEventsHistogram, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSortedEventsHistogram, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastSortedEventsHistogram = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  size = (this->sortedEventBufferCount)*(this->sortedEventBufferSize);
  CALLOC(dataSortedEvents, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSortedEvents, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastSortedEvents = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventCountsBuffer, 
    this->dataUnsortedEventCountsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventTargetsBuffer, 
    this->dataUnsortedEventTargetsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventDelaysBuffer, 
    this->dataUnsortedEventDelaysSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventWeightsBuffer, 
    this->dataUnsortedEventWeightsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventsHistogramBuffer, 
    this->dataUnsortedEventsHistogramSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSortedEventsReadBuffer, 
    this->dataSortedEventsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSortedEventsWriteBuffer, 
    this->dataSortedEventsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSortedEventsHistogramWriteBuffer, 
    this->dataSortedEventsHistogramSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSortedEventsHistogramReadBuffer, 
    this->dataSortedEventsHistogramSizeBytes);
    
  this->storeData(queue, block, SYNAPTIC_EVENTS_VALID_ALL);
}
/**************************************************************************************************/



void 
Data_SynapticEvents::initializeEventBuffers
(
  int       unsortedEventTimeSlotDelta,
  cl_uint   histogramBitShift,
  cl_uint   histogramBinMask,
  cl_uint   neuronCount,
  double    minDelay,
  double    percentDelayBoarder,
  double    percentTimeSlotDeltaDeviation,
  double    percentInhibitory
)
/**************************************************************************************************/
{
  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG_SIM("Data_SynapticEvents::initializeEventBuffers: set srand seed to " << this->srandSeed);
  
  std::stringstream ss;

  /*Initialize syn events and neuron counters*/
  for(cl_uint s = 0; s < this->timeSlots; s++)
  {
    cl_uint timeSlotTotal = 0;
    
    for(cl_uint b = 0; b < this->eventBufferCount; b++)
    {
      cl_uint eventCount = this->eventBufferSize;
      
      if(unsortedEventTimeSlotDelta >= 0)
      {
        eventCount = cl_uint(unsortedEventTimeSlotDelta*
          (this->timeSlots - s - 1)*(1 - percentTimeSlotDeltaDeviation/100.0) +
          abs(unsortedEventTimeSlotDelta*(this->timeSlots - s - 1)*
          (percentTimeSlotDeltaDeviation/100.0)*(double)rand()/((double)RAND_MAX)));
      }
      else
      {
        eventCount = cl_uint(this->eventBufferSize*(1 - percentTimeSlotDeltaDeviation/100.0) +
          abs(this->eventBufferSize*(percentTimeSlotDeltaDeviation/100.0)*
          (double)rand()/((double)RAND_MAX)));
      }
      
      if(eventCount >= this->eventBufferSize){eventCount = this->eventBufferSize;}

      this->dataUnsortedEventCounts[this->timeSlots*b + s] = eventCount;
      timeSlotTotal += eventCount;
      
      cl_uint inhibitorySynapseCount = cl_uint(eventCount*(percentInhibitory/100.0));
      
      
      this->initializeEventBuffer
      (
        1,
        eventCount,
        s,
        b,
        neuronCount,
        this->timeSlots,
        inhibitorySynapseCount,
        histogramBitShift,
        histogramBinMask,
        minDelay,
        percentDelayBoarder,
        this->dataUnsortedEventTargets,
        this->dataUnsortedEventDelays,
        this->dataUnsortedEventWeights,
        this->dataUnsortedEventsHistogram
      );
    }
    
    ss << s << "(" << timeSlotTotal << "), ";
  }
  
  LOG_SIM("Data_SynapticEvents::initializeEventBuffers: allocation of events to time slots " 
    << ss.str());
}
/**************************************************************************************************/



void 
Data_SynapticEvents::initializeEventBuffer
(
  cl_uint   keyOffset,
  cl_uint   eventCount,
  cl_uint   timeSlot,
  cl_uint   bufferId,
  cl_uint   neuronCount,
  cl_uint   timeSlotCount,
  cl_uint   inhibitorySynapseCount,
  cl_uint   histogramBitShift,
  cl_uint   histogramBinMask,
  double    minDelay,
  double    percentDelayBoarder,
  cl_uint   *dataUnsortedEventTargets,
  cl_uint   *dataUnsortedEventDelays,
  cl_uint   *dataUnsortedEventWeights,
  cl_uint   *dataUnsortedEventsHistogram
)
/**************************************************************************************************/
{
  for(cl_uint e = 0; e < eventCount; e++)
  {
    /*Compute pointer to event data*/
    cl_uint ptr = 
      /*Event data buffers*/
      bufferId * timeSlotCount * (this->eventBufferSize) +
      /*Current event data buffer*/
      timeSlot * (this->eventBufferSize) +
      /*Current event*/
      e;

    /*Compute event data*/
    /*target neuron*/
    cl_uint target_neuron = 
      cl_uint(abs(((neuronCount)-1)*((double)rand()/((double)RAND_MAX))));
    /*weight*/
    cl_float weight = 6.0f/1.4f;
    if (inhibitorySynapseCount != 0)
    {
      inhibitorySynapseCount--;
      weight = -67.0f/1.4f;
    }
    cl_float ev_weight = cl_float(weight*((double)rand()/((double)RAND_MAX)));
    /*delay: avoid close proximity to the boarders because of decrement on the host*/
    cl_float ev_delay = cl_float(minDelay*(percentDelayBoarder/100.0) + 
        abs(minDelay*(1.0-2*percentDelayBoarder/100.0)*((double)rand()/((double)RAND_MAX))));

    /*Store event*/
    dataUnsortedEventTargets[ptr] = target_neuron;
    *((cl_float *)(&(dataUnsortedEventWeights[ptr]))) = ev_weight;
    /*This reduction in accuracy is needed to match decrement operation on the host*/
    ev_delay += (timeSlot+1);
    *((cl_float *)(&(dataUnsortedEventDelays[ptr]))) = ev_delay;
    *((cl_float *)(&(dataUnsortedEventDelays[ptr]))) -= (timeSlot+1);

    /*Compute histogram key for target neuron based on MSBs*/
    cl_uint bin = 0;
    if(keyOffset == 0)
    {
      bin = (dataUnsortedEventTargets[ptr]>>histogramBitShift) & histogramBinMask;
    }
    else if(keyOffset == 1)
    {
      bin = (dataUnsortedEventDelays[ptr]>>histogramBitShift) & histogramBinMask;
    }
    else if(keyOffset == 2)
    {
      bin = (dataUnsortedEventWeights[ptr]>>histogramBitShift) & histogramBinMask;
    }
      
    /*Offset is based on time slot, bin, WG*/
    cl_uint offset = 
      /*Time slot*/
      timeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
      /*WG offset*/
      bufferId +
      /*time slot + bin with histogramBinSize as a pitch*/
      bin*(this->histogramBinSize);

    /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
    dataUnsortedEventsHistogram[offset]++;
  }
}
/**************************************************************************************************/



void
Data_SynapticEvents::groupEvents
/**************************************************************************************************/
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  bool    enableValues,
  size_t  histogramOutSize,
  cl_uint keyOffset,
  cl_uint totalBuffers,
  cl_uint bufferSize,
  cl_uint histogramBitMask,
  cl_uint histogramBitShift,
  cl_uint histogramOutBitShift,
  cl_uint histogramOutTotalGroups,
  cl_uint histogramOutBitMask,
  cl_uint histogramOutGridSize,
  cl_uint wgSize,
  cl_uint wiDataCount,
  cl_uint *dataUnsortedEventCounts,
  cl_uint *dataUnsortedEventTargets,
  cl_uint *dataUnsortedEventDelays,
  cl_uint *dataUnsortedEventWeights,
  cl_uint *dataHistogram,
  cl_uint *dataHistogramOut
){
  bool result = 1;
  std::stringstream ss;

#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  swap2(this->dataSortedEvents, this->dataPastSortedEvents);
#endif

  memset(this->dataSortedEvents, 0, this->dataSortedEventsSizeBytes);
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = (this->histogramBinCount)*(this->histogramBinSize) + 1;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogram, size*sizeof(cl_uint));
  
  /*Init data for verification*/
  cl_uint total_synaptic_events = 
    dataHistogram[(this->histogramBinCount)*(this->histogramBinSize)];
  
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = total_synaptic_events/(wgSize*wiDataCount);
  if(total_synaptic_event_chunks*(wgSize*wiDataCount) < total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = total_synaptic_event_chunks/histogramOutGridSize;
  if(wg_chunk_size*histogramOutGridSize < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (wgSize*wiDataCount);

  /*Group events for verification*/
  for(cl_uint b = 0; b < totalBuffers; b++)
  {
    cl_uint event_total = dataUnsortedEventCounts[b];
    
    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Compute pointer to event data*/
      cl_uint ptr = 
        /*Event data buffers*/
        b *
        (bufferSize) +
        /*Current event*/
        e;

      /*Access event*/
      cl_uint key = 0;
      if(keyOffset == 0)
      {
        key = dataUnsortedEventTargets[ptr];
      }
      else if(keyOffset == 1)
      {
        key = dataUnsortedEventDelays[ptr];
      }
      else if(keyOffset == 2)
      {
        key = dataUnsortedEventWeights[ptr];
      }
      
      /*Compute offset key for target neuron*/
      cl_uint bin = (key>>histogramBitShift)&histogramBitMask;
        
      /*Offset is based on time slot, bin, WG*/
      cl_uint bin_offset = 
      /*WG offset*/
      b +
      /*time slot + bin with histogramBinSize as a pitch*/
      bin*(this->histogramBinSize);

      /*Check for offset overlap*/
      if(dataOffsetGroupEventsCopy[bin_offset] >= dataHistogram[bin_offset+1])
      {
        ss << "Data_SynapticEvents::groupEvents: Destination event bin pointer "
          "overlaps with next pointer for bin " << bin << ", buffer " << b << std::endl;
        result = 0; 
        break;
      }
      
      /*Calculate offset in the grouped data space*/
      cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
        
      /*Compute histogram key for target neuron*/
      cl_uint hist_out_ptr = 
        /*WG offset for histogram out*/
        (dest_offset/wg_chunk_size) + 
        /*bin*/
        histogramOutTotalGroups*
        ((key>>histogramOutBitShift)&histogramOutBitMask);
      /*Verify*/
      if(hist_out_ptr > histogramOutSize)
      {
        ss << "Data_SynapticEvents::groupEvents: Pointer to an element in output "
          "histogram is outside of its range" << std::endl;
        result = 0;
        break;
      }
      /*Increment counter for this bin.*/
      dataHistogramOut[hist_out_ptr]++;
      
      /*Store event at its group location (grouped by bins)*/
      this->dataSortedEvents[dest_offset] = key;
      if(enableValues)
      {
        this->dataSortedEvents[(this->sortedEventBufferSize) + dest_offset] = ptr;
      }
      /*Increment ptr for next data item*/
      dataOffsetGroupEventsCopy[bin_offset]++;
    }
    if(!result){break;}
  }
  
  free(dataOffsetGroupEventsCopy);
  
  this->storeData(queue, block, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS);

  if(!result)
  {
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/



void 
Data_SynapticEvents::initializeHistogram
(
#if CLASS_VALIDATION_ENABLE
  cl_uint   eventDestinationBufferSize
#endif
)
/**************************************************************************************************/
{
  cl_uint print_bins = 0;
  
  /*Compute offsets based on counters.*/
  cl_uint max_offset = 0;
  
  print_bins ? std::cout << "Data_SynapticEvents::initializeHistogram: " << 
    "Number of synaptic events in bins: " << std::endl, true : false;
  
  for(cl_uint i = 0; i < this->timeSlots; i++)
  {
    cl_uint offset = i*((this->histogramBinCount)*(this->histogramBinSize) + 1);
    cl_uint runningSum = 0;
    cl_uint runningSize = 0;
    
    print_bins ? std::cout << "Data_SynapticEvents::initializeHistogram: Time slot " << i 
      << ": [", true : false;
    
    for(cl_uint j = 0; j < (this->histogramBinCount)*(this->histogramBinSize); j++)
    {
      cl_uint temp = this->dataUnsortedEventsHistogram[offset + j];
      this->dataUnsortedEventsHistogram[offset + j] = runningSum;
      runningSum += temp;

      if(j%(this->histogramBinSize) == 0 && j != 0)
      {
        print_bins ? std::cout << runningSum-runningSize << ", ", true : false;
        runningSize = runningSum;
      }
    }
    
    this->dataUnsortedEventsHistogram[offset + (this->histogramBinCount)*(this->histogramBinSize)] = 
      runningSum;
    
    print_bins ? std::cout << runningSum-runningSize << "] = " << runningSum << 
      std::endl, true : false;
    
    if(max_offset < runningSum){max_offset = runningSum;}
  }
#if CLASS_VALIDATION_ENABLE
  if(max_offset > eventDestinationBufferSize)
  {
    THROW_SIMEX("Data_SynapticEvents::initializeHistogram: Destination event buffer overflow. " 
      << "Need to increase eventDestinationBufferSize, which is currently " 
      << eventDestinationBufferSize << " above " << max_offset);
  }
#endif
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getUnsortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(timeSlot >= this->timeSlots)
  {
    throw SimException("Data_SynapticEvents::getUnsortedHistogramItem: time slot exceeds total time slots");
  }
  
  if(binID*(this->histogramBinSize) + bufferID > (this->histogramBinCount)*(this->histogramBinSize))
  {
    throw SimException("Data_SynapticEvents::getUnsortedHistogramItem: buffer ID and bin ID combination "
      "exceeds buffer count");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM);

  /*Pointer is based on time slot, bin, WG*/
  cl_uint ptr = 
  /*Time slot space = (time slot #) x (bins per WG) x WGs*/
  timeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
  /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
  binID*(this->histogramBinSize) + bufferID;
  
  if(type == Data_SynapticEvents::RECENT)
  {
    return this->dataUnsortedEventsHistogram[ptr];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    return this->dataPastUnsortedEventsHistogram[ptr];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getUnsortedHistogramItem: incorrect item type");
  }
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(backetID >= this->sortedEventsHistogramBacketCount)
  {
    throw SimException("Data_SynapticEvents::getSortedHistogramItem: backetID exceeds "
      "total backet count");
  }
  
  if(binID >= this->sortedEventsHistogramBinCount)
  {
    throw SimException("Data_SynapticEvents::getSortedHistogramItem: binID exceeds bin count");
  }
  
  if(itemID >= this->sortedEventsHistogramBinSize)
  {
    throw SimException("Data_SynapticEvents::getSortedHistogramItem: itemID exceeds bin size");
  }
#endif

  cl_uint ptr = 
    /*WG offset*/
    backetID*((this->sortedEventsHistogramBinSize)*(this->sortedEventsHistogramBinCount)) + 
    /*WG offset for histogram out*/
    itemID +
    /*bin*/
    binID*(this->sortedEventsHistogramBinSize);
  
  return this->getSortedHistogramItem(queue, ptr, type);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getSortedHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           pointer,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(pointer >= this->dataSortedEventsHistogramSize)
  {
    throw SimException("Data_SynapticEvents::getSortedHistogramItem: pointer exceeds "
      "total histogram size");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM);

  if(type == Data_SynapticEvents::RECENT)
  {
    return this->dataSortedEventsHistogram[pointer];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    return this->dataPastSortedEventsHistogram[pointer];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getSortedHistogramItem: incorrect item type");
  }
}
/**************************************************************************************************/



void
Data_SynapticEvents::getSortedHistogram
(
  cl::CommandQueue  &queue,
  cl_uint           **dataHistogramCopy,
#if CLASS_VALIDATION_ENABLE
  size_t            size,
#endif
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(size > this->dataSortedEventsHistogramSize)
  {
    throw SimException("Data_SynapticEvents::getSortedHistogram: insufficient buffer size");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM);

  if(type == Data_SynapticEvents::RECENT)
  {
    if(*dataHistogramCopy)
      free(*dataHistogramCopy);
    
    *dataHistogramCopy = (cl_uint *)calloc(this->dataSortedEventsHistogramSize, sizeof(cl_uint));
    memcpy(*dataHistogramCopy, this->dataSortedEventsHistogram, 
      this->dataSortedEventsHistogramSizeBytes);
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    if(*dataHistogramCopy)
      free(*dataHistogramCopy);
    
    *dataHistogramCopy = (cl_uint *)calloc(this->dataSortedEventsHistogramSize, sizeof(cl_uint));
    memcpy(*dataHistogramCopy, this->dataPastSortedEventsHistogram, 
      this->dataSortedEventsHistogramSizeBytes);
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getSortedHistogram: incorrect item type");
  }
}
/**************************************************************************************************/



void
Data_SynapticEvents::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataUnsortedEventCounts, this->dataPastUnsortedEventCounts);
#endif
    
    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventCountsBuffer, 
      this->dataUnsortedEventCountsSizeBytes, this->dataUnsortedEventCounts);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS;
  }
  
  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataUnsortedEventTargets, this->dataPastUnsortedEventTargets);
    swap2(this->dataUnsortedEventDelays, this->dataPastUnsortedEventDelays);
    swap2(this->dataUnsortedEventWeights, this->dataPastUnsortedEventWeights);
#endif

    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventTargetsBuffer, 
      this->dataUnsortedEventTargetsSizeBytes, this->dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventDelaysBuffer, 
      this->dataUnsortedEventDelaysSizeBytes, this->dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventWeightsBuffer, 
      this->dataUnsortedEventWeightsSizeBytes, this->dataUnsortedEventWeights);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS;
  }

  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataUnsortedEventsHistogram, this->dataPastUnsortedEventsHistogram);
#endif

    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventsHistogramBuffer, 
      this->dataUnsortedEventsHistogramSizeBytes, this->dataUnsortedEventsHistogram);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM;
  }

  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataSortedEvents, this->dataPastSortedEvents);
#endif

    ENQUEUE_READ_BUFFER(block, queue, this->dataSortedEventsReadBuffer, 
      this->dataSortedEventsSizeBytes, this->dataSortedEvents);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_SORTED_EVENTS;
  }

  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(this->dataSortedEventsHistogram, this->dataPastSortedEventsHistogram);
#endif

    ENQUEUE_READ_BUFFER(block, queue, this->dataSortedEventsHistogramReadBuffer, 
      this->dataSortedEventsHistogramSizeBytes, this->dataSortedEventsHistogram);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM;
  }
}
/**************************************************************************************************/



void
Data_SynapticEvents::storeData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataUnsortedEventCountsBuffer, 
      this->dataUnsortedEventCountsSizeBytes, this->dataUnsortedEventCounts);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS;
  }
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataUnsortedEventTargetsBuffer, 
      this->dataUnsortedEventTargetsSizeBytes, this->dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataUnsortedEventDelaysBuffer, 
      this->dataUnsortedEventDelaysSizeBytes, this->dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataUnsortedEventWeightsBuffer, 
      this->dataUnsortedEventWeightsSizeBytes, this->dataUnsortedEventWeights);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS;
  }
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataUnsortedEventsHistogramBuffer, 
      this->dataUnsortedEventsHistogramSizeBytes, this->dataUnsortedEventsHistogram);

    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM;
  }
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_SORTED_EVENTS)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataSortedEventsReadBuffer, 
      this->dataSortedEventsSizeBytes, this->dataSortedEvents);

    this->dataValid |= SYNAPTIC_EVENTS_VALID_SORTED_EVENTS;
  }
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataSortedEventsHistogramReadBuffer, 
      this->dataSortedEventsHistogramSizeBytes, this->dataSortedEventsHistogram);

    this->dataValid |= SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM;
  }
}
/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
