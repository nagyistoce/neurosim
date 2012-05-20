
/* ===============================================================================================



  =============================================================================================== */


  
/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "SynapticEvents.hpp"

/**************************************************************************************************/



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



SynapticEvents::SynapticEvents()
/**************************************************************************************************/
{
  this->reset(false);
}
/**************************************************************************************************/



SynapticEvents::~SynapticEvents()
/**************************************************************************************************/
{
  this->reset(true);
}
/**************************************************************************************************/



void
SynapticEvents::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  cl_uint                             timeSlots,
  cl_uint                             eventBufferCount,
  cl_uint                             eventBufferSize,
  cl_uint                             histogramBinCount,
  cl_uint                             histogramBinSize,
  struct kernelStatistics             *kernelStats,
  std::stringstream                   *dataToSimulationLogFile,
  std::stringstream                   *dataToReportLogFile
)
/**************************************************************************************************/
{
  if(!this->resetObject)
  {
    throw SimException("SynapticEvents::initialize: attemp to initialize without reset");
  }
  
  this->resetObject = false;
  this->dataValid = 0;
  this->timeSlots = timeSlots;
  this->eventBufferCount = eventBufferCount;
  this->eventBufferSize = eventBufferSize;
  this->histogramBinCount = histogramBinCount;
  this->histogramBinSize = histogramBinSize;
  this->dataToSimulationLogFile = dataToSimulationLogFile;
  this->dataToReportLogFile = dataToReportLogFile;
  
  size_t size = (this->eventBufferCount)*(this->timeSlots);
  
  CALLOC_O(dataUnsortedEventCounts, cl_uint, (this->eventBufferCount)*(this->timeSlots));
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventCounts, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventCounts = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  size = (this->eventBufferCount)*(this->timeSlots)*(this->eventBufferSize);
  
  CALLOC_O(dataUnsortedEventTargets, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventTargets, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventTargets = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif

  CALLOC_O(dataUnsortedEventDelays, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventDelays, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventDelays = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  CALLOC_O(dataUnsortedEventWeights, cl_uint, (this->eventBufferCount)*(this->timeSlots)*
    (this->eventBufferSize));
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventWeights, kernelStats);
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  this->dataPastUnsortedEventWeights = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  size = (this->timeSlots)*((this->histogramBinCount)*(this->histogramBinSize) + 1);
  
  CALLOC_O(dataHistogram, cl_uint, size);
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
  
  this->storeBuffers(queue, block, SYNAPTIC_EVENTS_VALID_ALL);
}
/**************************************************************************************************/



void 
SynapticEvents::clearUnsortedEvents
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           clearBitMask
)
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif
  
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
    
    this->storeBuffers(queue, block, SYNAPTIC_EVENTS_VALID_ALL);
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
SynapticEvents::setUnsortedEvents
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_bool             initHistogram,
  cl_uint             unsortedEventTimeSlotDelta,
  cl_uint             eventDestinationBufferSize,
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
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif
  
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
      eventDestinationBufferSize
    );
  }

  this->storeBuffers(queue, block, SYNAPTIC_EVENTS_VALID_ALL);
}
/**************************************************************************************************/



cl_uint
SynapticEvents::getEventCount
(
  cl::CommandQueue  &queue,
  cl_uint           bufferID,
  cl_uint           timeSlot,
  const int         type
)
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
  
  if(bufferID >= this->eventBufferCount)
  {
    throw SimException("SynapticEvents::getEventCount: buffer ID exceeds buffer count");
  }
  
  if(timeSlot >= this->timeSlots)
  {
    throw SimException("SynapticEvents::getEventCount: time slot exceeds total time slots");
  }
#endif

  this->getEventBuffers(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS);
  
  if(type == SynapticEvents::RECENT)
  {
    return this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == SynapticEvents::PREVIOUS)
  {
    return this->dataPastUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot];
  }
#endif
  else
  {
    throw SimException("SynapticEvents::getEventCount: incorrect event count type");
  }
}
/**************************************************************************************************/



void
SynapticEvents::getEvent
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
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
  
  if(bufferID >= this->eventBufferCount)
  {
    throw SimException("SynapticEvents::getEvent: buffer ID exceeds buffer count");
  }
  
  if(timeSlot >= this->timeSlots)
  {
    throw SimException("SynapticEvents::getEvent: time slot exceeds total time slots");
  }
#endif

#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->getEventBuffers(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS | 
    SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS);
  
  if
  (
    (eventID >= this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot]) &&
     (eventID != 0) &&
     (this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot] != 0)
  )
  {
    std::stringstream ss;
    ss << "SynapticEvents::getEvent: event ID " << eventID << " exceeds event count " 
      << this->dataUnsortedEventCounts[(this->timeSlots)*bufferID + timeSlot] << std::endl;
    throw SimException(ss.str());
  }
#else
  this->getEventBuffers(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS);
#endif

  cl_uint ptr = 
    /*Event data buffers*/
    bufferID * (this->timeSlots) * (this->eventBufferSize) +
    /*Current event data buffer*/
    timeSlot * (this->eventBufferSize) +
    /*Current event*/
    eventID;
  
  if(type == SynapticEvents::RECENT)
  {
    targetNeuron = this->dataUnsortedEventTargets[ptr];
    eventTime = this->dataUnsortedEventDelays[ptr];
    weight = this->dataUnsortedEventWeights[ptr];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == SynapticEvents::PREVIOUS)
  {
    targetNeuron = this->dataPastUnsortedEventTargets[ptr];
    eventTime = this->dataPastUnsortedEventDelays[ptr];
    weight = this->dataPastUnsortedEventWeights[ptr];
  }
#endif
  else
  {
    throw SimException("SynapticEvents::getEvent: incorrect event type");
  }
}
/**************************************************************************************************/



cl_uint
SynapticEvents::getHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID,
  const int         type
)
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();

  if(timeSlot >= this->timeSlots)
  {
    throw SimException("SynapticEvents::getHistogramItem: time slot exceeds total time slots");
  }
  
  if(binID >= this->histogramBinCount)
  {
    throw SimException("SynapticEvents::getHistogramItem: buffer ID exceeds buffer count");
  }
  
  if(bufferID >= this->histogramBinSize)
  {
    throw SimException("SynapticEvents::getHistogramItem: buffer ID exceeds buffer count");
  }
#endif

  this->getEventBuffers(queue, CL_TRUE, SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM);

  /*Pointer is based on time slot, bin, WG*/
  cl_uint ptr = 
  /*Time slot space = (time slot #) x (bins per WG) x WGs*/
  timeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
  /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
  binID*(this->histogramBinSize) + bufferID;
  
  if(type == SynapticEvents::RECENT)
  {
    return this->dataHistogram[ptr];
  }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  else if(type == SynapticEvents::PREVIOUS)
  {
    return this->dataPastHistogram[ptr];
  }
#endif
  else
  {
    throw SimException("SynapticEvents::getHistogramItem: incorrect item type");
  }
}
/**************************************************************************************************/



void 
SynapticEvents::invalidateEvents
()
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= 
    (
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS  ^ 0xFFFFFFFF) & 
      (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS        ^ 0xFFFFFFFF) & 
      (SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM         ^ 0xFFFFFFFF)
    );
}
/**************************************************************************************************/



void 
SynapticEvents::invalidateHistogram
()
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= (SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
SynapticEvents::refresh
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->getEventBuffers
  (
    queue, 
    block, 
    SYNAPTIC_EVENTS_VALID_ALL
  );
}
/**************************************************************************************************/



void
SynapticEvents::randomizeHistogram
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             timeSlot
)
/**************************************************************************************************/
{
  cl_uint maxElementSize = (0xFFFFFFFF/((this->histogramBinCount)*(this->histogramBinSize)));
  
  for(cl_uint j = 0; j < (this->histogramBinCount); j++)
  {
    cl_uint offset = timeSlot*((this->histogramBinCount)*(this->histogramBinSize) + 1) + 
      j*(this->histogramBinSize);
    
    for(cl_uint k = 0; k < (this->histogramBinSize); k++)
    {
      this->dataHistogram[offset + k] = 
        cl_uint(abs(maxElementSize*((double)rand()/((double)RAND_MAX))));
    }
  }
  
  this->storeBuffers(queue, block, SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM);
}
/**************************************************************************************************/



cl_uint 
SynapticEvents::getEventBufferCount
()
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  return eventBufferCount;
}
/**************************************************************************************************/



cl_uint 
SynapticEvents::getEventTimeSlots
()
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  return timeSlots;
}
/**************************************************************************************************/



cl_uint 
SynapticEvents::getEventBufferSize
()
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  return eventBufferSize;
}
/**************************************************************************************************/



cl_uint 
SynapticEvents::getEventHistogramBinCount
()
/**************************************************************************************************/
{
#if SYNAPTIC_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  return histogramBinCount;
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void 
SynapticEvents::initializeEventBuffers
(
  cl_uint   unsortedEventTimeSlotDelta,
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
  SET_RANDOM_SEED(srandSeed, srandCounter);
  LOG("SynapticEvents::initializeEventBuffers: set srand seed to " << srandSeed, 0);
  
  std::stringstream ss;

  /*Initialize syn events and neuron counters*/
  for(cl_uint s = 0; s < this->timeSlots; s++)
  {
    cl_uint timeSlotTotal = 0;
    
    for(cl_uint b = 0; b < this->eventBufferCount; b++)
    {
      cl_uint event_total = this->eventBufferSize;
      
      if(unsortedEventTimeSlotDelta != NULL)
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
  
  LOG("SynapticEvents::initializeEventBuffers: allocation of events to time slots " << ss.str(), 0);
}
/**************************************************************************************************/



void 
SynapticEvents::initializeHistogram
(
  cl_uint   eventDestinationBufferSize
)
/**************************************************************************************************/
{
  cl_uint print_bins = 0;
  
  /*Compute offsets based on counters.*/
  cl_uint max_offset = 0;
  
  print_bins ? std::cout << "SynapticEvents::initializeHistogram: " << 
    "Number of synaptic events in bins: " << std::endl, true : false;
  
  for(cl_uint i = 0; i < this->timeSlots; i++)
  {
    cl_uint offset = i*((this->histogramBinCount)*(this->histogramBinSize) + 1);
    cl_uint runningSum = 0;
    cl_uint runningSize = 0;
    
    print_bins ? std::cout << "SynapticEvents::initializeHistogram: Time slot " << i 
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

  if(max_offset > eventDestinationBufferSize)
  {
    std::stringstream ss;
    ss << "SynapticEvents::initializeHistogram: Destination event buffer overflow. " << 
      "Need to increase eventDestinationBufferSize, which is currently " << 
      eventDestinationBufferSize << " above " << max_offset << std::endl;
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/



void
SynapticEvents::getEventBuffers
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  if
  (
    (selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS) && 
    !((this->dataValid) & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS)
  )
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(dataUnsortedEventCounts, dataPastUnsortedEventCounts);
#endif
    
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataUnsortedEventCountsBuffer, 
      this->dataUnsortedEventCountsSizeBytes, this->dataUnsortedEventCounts);
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS;
  }
  
  if
  (
    (selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS) && 
    !((this->dataValid) & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS)
  )
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(dataUnsortedEventTargets, dataPastUnsortedEventTargets);
    swap2(dataUnsortedEventDelays, dataPastUnsortedEventDelays);
    swap2(dataUnsortedEventWeights, dataPastUnsortedEventWeights);
#endif

    ENQUEUE_READ_BUFFER_O(block, queue, this->dataUnsortedEventTargetsBuffer, 
      this->dataUnsortedEventTargetsSizeBytes, this->dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataUnsortedEventDelaysBuffer, 
      this->dataUnsortedEventDelaysSizeBytes, this->dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataUnsortedEventWeightsBuffer, 
      this->dataUnsortedEventWeightsSizeBytes, this->dataUnsortedEventWeights);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS;
  }
  
  if
  (
    (selectBitMask & SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM) && 
    !((this->dataValid) & SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM)
  )
  {
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    swap2(dataHistogram, dataPastHistogram);
#endif

    ENQUEUE_READ_BUFFER_O(block, queue, this->dataHistogramBuffer, 
      this->dataHistogramSizeBytes, this->dataHistogram);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM;
  }
}
/**************************************************************************************************/



void
SynapticEvents::reset
(
  bool checkForNull
)
/**************************************************************************************************/
{
  try
  {
  this->resetObject = true;
  this->dataValid = 0;
  this->srandSeed = 0;
  this->srandCounter = 0;
  this->dataToSimulationLogFile = NULL;
  this->dataToReportLogFile = NULL;
  this->eventBufferCount = 0;
  this->eventBufferSize = 0;
  this->timeSlots = 0;
  this->histogramBinCount = 0;
  this->histogramBinSize = 0;
  
  this->dataUnsortedEventTargetsSize;
  this->dataUnsortedEventTargetsSizeBytes;
  
  this->dataUnsortedEventDelaysSize;
  this->dataUnsortedEventDelaysSizeBytes;
  
  this->dataUnsortedEventWeightsSize;
  this->dataUnsortedEventWeightsSizeBytes;
  
  this->dataUnsortedEventCountsSize;
  this->dataUnsortedEventCountsSizeBytes;
  
  this->dataHistogramSize;
  this->dataHistogramSizeBytes;
  
  if(checkForNull)
  {
    if(this->dataUnsortedEventCounts)
    {
      free(this->dataUnsortedEventCounts);
      this->dataUnsortedEventCounts = NULL;
    }
    if(this->dataUnsortedEventWeights)
    {
      free(this->dataUnsortedEventWeights);
      this->dataUnsortedEventWeights = NULL;
    }
    if(this->dataUnsortedEventDelays)
    {
      free(this->dataUnsortedEventDelays);
      this->dataUnsortedEventDelays = NULL;
    }
    if(this->dataUnsortedEventTargets)
    {
      free(this->dataUnsortedEventTargets);
      this->dataUnsortedEventTargets = NULL;
    }
    if(this->dataHistogram)
    {
      free(this->dataHistogram);
      this->dataHistogram = NULL;
    }
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    if(this->dataPastUnsortedEventCounts)
    {
      free(this->dataPastUnsortedEventCounts);
      this->dataPastUnsortedEventCounts = NULL;
    }
    if(this->dataPastUnsortedEventWeights)
    {
      free(this->dataPastUnsortedEventWeights);
      this->dataPastUnsortedEventWeights = NULL;
    }
    if(this->dataPastUnsortedEventDelays)
    {
      free(this->dataPastUnsortedEventDelays);
      this->dataPastUnsortedEventDelays = NULL;
    }
    if(this->dataPastUnsortedEventTargets)
    {
      free(this->dataPastUnsortedEventTargets);
      this->dataPastUnsortedEventTargets = NULL;
    }
    if(this->dataPastHistogram)
    {
      free(this->dataPastHistogram);
      this->dataPastHistogram = NULL;
    }
#endif
  }
  else
  {
    this->dataUnsortedEventCounts = NULL;
    this->dataUnsortedEventWeights = NULL;
    this->dataUnsortedEventDelays = NULL;
    this->dataUnsortedEventTargets = NULL;
    this->dataHistogram = NULL;
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    this->dataPastUnsortedEventCounts = NULL;
    this->dataPastUnsortedEventWeights = NULL;
    this->dataPastUnsortedEventDelays = NULL;
    this->dataPastUnsortedEventTargets = NULL;
    this->dataPastHistogram = NULL;
#endif
  }
  }
  CATCH(std::cerr, SynapticEvents::reset, throw SimException("SynapticEvents::reset: failed.");)
}
/**************************************************************************************************/



void
SynapticEvents::storeBuffers
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataUnsortedEventCountsBuffer, 
      this->dataUnsortedEventCountsSizeBytes, this->dataUnsortedEventCounts);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS;
  }
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataUnsortedEventTargetsBuffer, 
      this->dataUnsortedEventTargetsSizeBytes, this->dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataUnsortedEventDelaysBuffer, 
      this->dataUnsortedEventDelaysSizeBytes, this->dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataUnsortedEventWeightsBuffer, 
      this->dataUnsortedEventWeightsSizeBytes, this->dataUnsortedEventWeights);
      
    this->dataValid |= SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS;
  }
  
  if(selectBitMask & SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataHistogramBuffer, 
      this->dataHistogramSizeBytes, this->dataHistogram);

    this->dataValid |= SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM;
  }
}
/**************************************************************************************************/



void
SynapticEvents::isInitialized
()
/**************************************************************************************************/
{
  if(this->resetObject)
  {
    throw SimException("SynapticEvents::isInitialized: the object was not initialized");
  }
}
/**************************************************************************************************/
