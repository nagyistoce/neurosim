
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
    swap2(dataUnsortedEventCounts, dataPastUnsortedEventCounts);
    swap2(dataUnsortedEventTargets, dataPastUnsortedEventTargets);
    swap2(dataUnsortedEventDelays, dataPastUnsortedEventDelays);
    swap2(dataUnsortedEventWeights, dataPastUnsortedEventWeights);
    swap2(dataHistogram, dataPastHistogram);
#endif

    memset(this->dataUnsortedEventCounts, 0, this->dataUnsortedEventCountsSizeBytes);
    memset(this->dataUnsortedEventTargets, 0, this->dataUnsortedEventTargetsSizeBytes);
    memset(this->dataUnsortedEventDelays, 0, this->dataUnsortedEventDelaysSizeBytes);
    memset(this->dataUnsortedEventWeights, 0, this->dataUnsortedEventWeightsSizeBytes);
    memset(this->dataHistogram, 0, this->dataHistogramSizeBytes);
    
    this->storeData(queue, block, SYNAPTIC_EVENTS_VALID_ALL);
  }
  
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  if(clearBitMask & 0x2)
  {
    memset(this->dataPastUnsortedEventCounts, 0, this->dataUnsortedEventCountsSizeBytes);
    memset(this->dataPastUnsortedEventTargets, 0, this->dataUnsortedEventTargetsSizeBytes);
    memset(this->dataPastUnsortedEventDelays, 0, this->dataUnsortedEventDelaysSizeBytes);
    memset(this->dataPastUnsortedEventWeights, 0, this->dataUnsortedEventWeightsSizeBytes);
    memset(this->dataPastHistogram, 0, this->dataHistogramSizeBytes);
  }
#endif
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
  swap2(dataUnsortedEventCounts, dataPastUnsortedEventCounts);
  swap2(dataUnsortedEventTargets, dataPastUnsortedEventTargets);
  swap2(dataUnsortedEventDelays, dataPastUnsortedEventDelays);
  swap2(dataUnsortedEventWeights, dataPastUnsortedEventWeights);
  swap2(dataHistogram, dataPastHistogram);
#endif
  
  memset(dataUnsortedEventCounts, 0, dataUnsortedEventCountsSizeBytes);
  memset(dataUnsortedEventTargets, 0, dataUnsortedEventTargetsSizeBytes);
  memset(dataUnsortedEventDelays, 0, dataUnsortedEventDelaysSizeBytes);
  memset(dataUnsortedEventWeights, 0, dataUnsortedEventWeightsSizeBytes);
  memset(dataHistogram, 0, dataHistogramSizeBytes);

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

  this->storeData(queue, block, SYNAPTIC_EVENTS_VALID_ALL);
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
#else
  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS);
#endif

  cl_uint ptr = 
    /*Event data buffers*/
    bufferID * (this->timeSlots) * (this->eventBufferSize) +
    /*Current event data buffer*/
    timeSlot * (this->eventBufferSize) +
    /*Current event*/
    eventID;
  
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



cl_uint
Data_SynapticEvents::getHistogramItem
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
    throw SimException("Data_SynapticEvents::getHistogramItem: time slot exceeds total time slots");
  }
  
  if(binID >= this->histogramBinCount)
  {
    throw SimException("Data_SynapticEvents::getHistogramItem: buffer ID exceeds buffer count");
  }
  
  if(bufferID >= this->histogramBinSize)
  {
    throw SimException("Data_SynapticEvents::getHistogramItem: buffer ID exceeds buffer count");
  }
#endif

  this->getData(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM);

  /*Pointer is based on time slot, bin, WG*/
  cl_uint ptr = 
  /*Time slot space = (time slot #) x (bins per WG) x WGs*/
  timeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
  /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
  binID*(this->histogramBinSize) + bufferID;
  
  if(type == Data_SynapticEvents::RECENT)
  {
    return this->dataHistogram[ptr];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == Data_SynapticEvents::PREVIOUS)
  {
    return this->dataPastHistogram[ptr];
  }
#endif
  else
  {
    throw SimException("Data_SynapticEvents::getHistogramItem: incorrect item type");
  }
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getPreviousHistogramItem
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
  return mySelf->getHistogramItem(queue, timeSlot, binID, bufferID, Data_SynapticEvents::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Data_SynapticEvents::getCurrentHistogramItem
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
  return mySelf->getHistogramItem(queue, timeSlot, binID, bufferID, Data_SynapticEvents::RECENT);
}
/**************************************************************************************************/



void 
Data_SynapticEvents::invalidateEvents
()
/**************************************************************************************************/
{
  this->dataValid &= 
    (
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS  ^ 0xFFFFFFFF) & 
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS        ^ 0xFFFFFFFF) & 
      (SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM         ^ 0xFFFFFFFF)
    );
}
/**************************************************************************************************/



void 
Data_SynapticEvents::invalidateHistogram
()
/**************************************************************************************************/
{
  this->dataValid &= (SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
Data_SynapticEvents::refresh
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
    SYNAPTIC_EVENTS_VALID_ALL
  );
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventBufferCount
()
/**************************************************************************************************/
{
  return eventBufferCount;
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventTimeSlots
()
/**************************************************************************************************/
{
  return timeSlots;
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventBufferSize
()
/**************************************************************************************************/
{
  return eventBufferSize;
}
/**************************************************************************************************/



cl_uint 
Data_SynapticEvents::getEventHistogramBinCount
()
/**************************************************************************************************/
{
  return histogramBinCount;
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
  
  CALLOC(dataHistogram, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataHistogram, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastHistogram = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventCountsBuffer, 
    this->dataUnsortedEventCountsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventTargetsBuffer, 
    this->dataUnsortedEventTargetsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventDelaysBuffer, 
    this->dataUnsortedEventDelaysSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataUnsortedEventWeightsBuffer, 
    this->dataUnsortedEventWeightsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataHistogramBuffer, 
    this->dataHistogramSizeBytes);
  
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
      cl_uint event_total = this->eventBufferSize;
      
      if(unsortedEventTimeSlotDelta >= 0)
      {
        event_total = cl_uint(unsortedEventTimeSlotDelta*
          (this->timeSlots - s - 1)*(1 - percentTimeSlotDeltaDeviation/100.0) +
          abs(unsortedEventTimeSlotDelta*(this->timeSlots - s - 1)*
          (percentTimeSlotDeltaDeviation/100.0)*(double)rand()/((double)RAND_MAX)));
      }
      else
      {
        event_total = cl_uint(this->eventBufferSize*(1 - percentTimeSlotDeltaDeviation/100.0) +
          abs(this->eventBufferSize*(percentTimeSlotDeltaDeviation/100.0)*
          (double)rand()/((double)RAND_MAX)));
      }
      
      if(event_total >= this->eventBufferSize){event_total = this->eventBufferSize;}

      this->dataUnsortedEventCounts[this->timeSlots*b + s] = event_total;
      timeSlotTotal += event_total;
      
      cl_uint count_inhibitory = cl_uint(event_total*(percentInhibitory/100.0));
      
      for(cl_uint e = 0; e < event_total; e++)
      {
        /*Compute pointer to event data*/
        cl_uint ptr = 
          /*Event data buffers*/
          b * this->timeSlots * 
          (this->eventBufferSize) +
          /*Current event data buffer*/
          s * 
          (this->eventBufferSize) +
          /*Current event*/
          e;

        /*Compute event data*/
        /*target neuron*/
        cl_uint target_neuron = 
          cl_uint(abs(((neuronCount)-1)*((double)rand()/((double)RAND_MAX))));
        /*weight*/
        cl_float weight = 6.0f/1.4f;
        if (count_inhibitory != 0)
        {
          count_inhibitory--;
          weight = -67.0f/1.4f;
        }
        cl_float ev_weight = cl_float(weight*((double)rand()/((double)RAND_MAX)));
        /*delay: avoid close proximity to the boarders because of decrement on the host*/
        cl_float ev_delay = 
          cl_float(minDelay*(percentDelayBoarder/100.0) + 
            abs(minDelay*(1.0-2*percentDelayBoarder/100.0)*((double)rand()/((double)RAND_MAX))));

        /*Store event*/
        this->dataUnsortedEventTargets[ptr] = target_neuron;
        *((cl_float *)(&(this->dataUnsortedEventWeights[ptr]))) = ev_weight;
        /*This reduction in accuracy is needed to match decrement operation on the host*/
        ev_delay += (s+1);
        *((cl_float *)(&(this->dataUnsortedEventDelays[ptr]))) = ev_delay;
        *((cl_float *)(&(this->dataUnsortedEventDelays[ptr]))) -= (s+1);
        
        /*Compute histogram key for target neuron based on MSBs*/
        cl_uint bin = ((this->dataUnsortedEventDelays[ptr]) >> (histogramBitShift)) & 
          histogramBinMask;
        /*Offset is based on time slot, bin, WG*/
        cl_uint offset = 
        /*Time slot*/
        s*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
        /*WG offset*/
        b +
        /*time slot + bin with histogramBinSize as a pitch*/
        bin*(this->histogramBinSize);

        /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
        this->dataHistogram[offset]++;
      }
    }
    
    ss << s << "(" << timeSlotTotal << "), ";
  }
  
  LOG_SIM("Data_SynapticEvents::initializeEventBuffers: allocation of events to time slots " 
    << ss.str());
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
      cl_uint temp = this->dataHistogram[offset + j];
      this->dataHistogram[offset + j] = runningSum;
      runningSum += temp;

      if(j%(this->histogramBinSize) == 0 && j != 0)
      {
        print_bins ? std::cout << runningSum-runningSize << ", ", true : false;
        runningSize = runningSum;
      }
    }
    
    this->dataHistogram[offset + (this->histogramBinCount)*(this->histogramBinSize)] = 
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
    swap2(dataUnsortedEventCounts, dataPastUnsortedEventCounts);
#endif
    
    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventCountsBuffer, 
      this->dataUnsortedEventCountsSizeBytes, this->dataUnsortedEventCounts);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS;
  }
  
  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(dataUnsortedEventTargets, dataPastUnsortedEventTargets);
    swap2(dataUnsortedEventDelays, dataPastUnsortedEventDelays);
    swap2(dataUnsortedEventWeights, dataPastUnsortedEventWeights);
#endif

    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventTargetsBuffer, 
      this->dataUnsortedEventTargetsSizeBytes, this->dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventDelaysBuffer, 
      this->dataUnsortedEventDelaysSizeBytes, this->dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(block, queue, this->dataUnsortedEventWeightsBuffer, 
      this->dataUnsortedEventWeightsSizeBytes, this->dataUnsortedEventWeights);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS;
  }
  

  IF_HIT_READ(selectBitMask, SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM, this->dataValid)
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(dataHistogram, dataPastHistogram);
#endif

    ENQUEUE_READ_BUFFER(block, queue, this->dataHistogramBuffer, 
      this->dataHistogramSizeBytes, this->dataHistogram);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM;
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
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataHistogramBuffer, 
      this->dataHistogramSizeBytes, this->dataHistogram);

    this->dataValid |= SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM;
  }
}
/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
