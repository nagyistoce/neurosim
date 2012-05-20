
/* ===============================================================================================



  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "SpikeEvents.hpp"

/**************************************************************************************************/



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



SpikeEvents::SpikeEvents()
/**************************************************************************************************/
{
  this->reset(false);
}
/**************************************************************************************************/



SpikeEvents::~SpikeEvents()
/**************************************************************************************************/
{
  this->reset(true);
}
/**************************************************************************************************/



void
SpikeEvents::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  cl_uint                             simulationStepSize,
  cl_uint                             spikeDatumSize,
  cl_uint                             spikePacketSizeWords,
  cl_uint                             spikePacketSize,
  cl_uint                             spikePackets,
  cl_uint                             neuronCount,
  struct kernelStatistics             *kernelStats,
  std::stringstream                   *dataToSimulationLogFile,
  std::stringstream                   *dataToReportLogFile
)
/**************************************************************************************************/
{
  if(!this->resetObject)
  {
    throw SimException("SpikeEvents::initialize: attemp to initialize without reset");
  }
  
  this->resetObject = false;
  this->dataValid = 0;
  this->neuronCount = neuronCount;
  this->simulationStepSize = simulationStepSize;
  this->spikeDatumSize = spikeDatumSize;
  this->spikePackets = spikePackets;
  this->spikePacketSize = spikePacketSize;
  this->spikePacketSizeWords = spikePacketSizeWords;
  this->dataToSimulationLogFile = dataToSimulationLogFile;
  this->dataToReportLogFile = dataToReportLogFile;
  
  this->dataPastSpikePacketCounts = (cl_uint *)calloc(this->spikePackets, sizeof(cl_uint));
  this->dataPastSpikePackets = (cl_uint *)calloc((this->spikePackets)*(this->spikePacketSizeWords), 
    sizeof(cl_uint));

  CALLOC_O(dataSpikePackets, cl_uint, (this->spikePackets)*(this->spikePacketSizeWords));
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSpikePackets, kernelStats);
  CALLOC_O(dataSpikePacketCounts, cl_uint, this->spikePackets);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSpikePacketCounts, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSpikePacketsBuffer, 
    this->dataSpikePacketsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSpikePacketCountsBuffer, 
    this->dataSpikePacketCountsSizeBytes);
    
  this->storeBuffers(queue, block, SPIKE_EVENTS_VALID_SPIKES);
    
#if EXPAND_EVENTS_ENABLE
  this->setKernelArguments = true;
  
  /* register device local memory buffer for stats */
  cl_uint lmExpandEvents = sizeof(cl_uint)*EXPAND_EVENTS_CACHE_SIZE_WORDS; 
  lmExpandEventsSizeBytes = lmExpandEvents;
  REGISTER_MEMORY_O(device, EXPAND_EVENTS_KERNEL_NAME, MEM_LOCAL, lmExpandEvents, kernelStats);
  
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC_O(dataExpandEventsDebugHost, cl_uint, EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataExpandEventsDebugHost, kernelStats);
  /* allocate memory for debug device buffer */
  CALLOC_O(dataExpandEventsDebugDevice, cl_uint, EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataExpandEventsDebugDevice, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataExpandEventsDebugHostBuffer, 
    this->dataExpandEventsDebugHostSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataExpandEventsDebugDeviceBuffer, 
    this->dataExpandEventsDebugDeviceSizeBytes);
    
  this->storeBuffers(queue, block, SPIKE_EVENTS_VALID_EXPAND_DEBUG);
#endif

#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC_O(dataExpandEventsError, cl_uint, EXPAND_EVENTS_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataExpandEventsError, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataExpandEventsErrorBuffer, 
    this->dataExpandEventsErrorSizeBytes);
    
  this->storeBuffers(queue, block, SPIKE_EVENTS_VALID_EXPAND_ERROR);
#endif

  this->blockSizeX_KernelExpandEvents = EXPAND_EVENTS_WG_SIZE_WI;
  this->blockSizeY_KernelExpandEvents = 1;
  this->argNumExpandEvents = 0;
  
  this->globalThreadsExpandEvents = 
    new cl::NDRange(EXPAND_EVENTS_WG_SIZE_WI*EXPAND_EVENTS_GRID_SIZE_WG);
  this->localThreadsExpandEvents = 
    new cl::NDRange(this->blockSizeX_KernelExpandEvents, this->blockSizeY_KernelExpandEvents);
    
  createKernel
  (
    context,
    device,
    this->kernelExpandEvents,
    EXPAND_EVENTS_KERNEL_FILE_NAME,
    EXPAND_EVENTS_KERNEL_NAME,
    "",
    this->blockSizeX_KernelExpandEvents,
    this->blockSizeY_KernelExpandEvents
  );
#endif
}
/**************************************************************************************************/



void 
SpikeEvents::clearEvents
(
  cl::CommandQueue  &queue,
  cl_bool           block
)
/**************************************************************************************************/
{
  this->isInitialized();
  
  swap2(dataSpikePackets, dataPastSpikePackets);
  swap2(dataSpikePacketCounts, dataPastSpikePacketCounts);
  
  memset(this->dataSpikePackets, 0, this->dataSpikePacketsSizeBytes);
  memset(this->dataSpikePacketCounts, 0, this->dataSpikePacketCountsSizeBytes);
  
  this->storeBuffers(queue, block, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



void 
SpikeEvents::setEvents
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  double              spikeBufferMinPercentFill,
  double              spikeBufferMaxPercentFill,
  double              spikeNeuronsPermil
)
/**************************************************************************************************/
{
  this->isInitialized();
  
  swap2(dataSpikePackets, dataPastSpikePackets);
  swap2(dataSpikePacketCounts, dataPastSpikePacketCounts);
  
  memset(this->dataSpikePackets, 0, this->dataSpikePacketsSizeBytes);
  memset(this->dataSpikePacketCounts, 0, this->dataSpikePacketCountsSizeBytes);
  
  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG("SpikeEvents::setEvents: set srand seed to " << this->srandSeed, 0);

  cl_uint neuronsPerPacket = (this->neuronCount)/(this->spikePackets);
  cl_uint maxSpikesPerPacket = neuronsPerPacket;
  double targetPacketSpikes = 0;
  cl_uint totalNetworkSpikes = 0;

  if(spikeNeuronsPermil != NULL)
  {
    maxSpikesPerPacket = cl_uint(ceil(1.5*(double)maxSpikesPerPacket*(spikeNeuronsPermil/1000)));
    targetPacketSpikes = (double)neuronsPerPacket*(spikeNeuronsPermil/1000);
    if(maxSpikesPerPacket < 1){maxSpikesPerPacket = 1;}
    
    if(maxSpikesPerPacket > this->spikePacketSize)
    {
      std::stringstream ss;
      ss << "SpikeEvents::setEvents: can't fit spikes within " << this->spikePacketSize 
        << " boundary" << std::endl;
      throw SimException(ss.str());
    }
  }

  LOG("SpikeEvents::setEvents: max spikes per packet is " << maxSpikesPerPacket, 0);

  /* init spike data */
  for(cl_uint packet = 0; packet < this->spikePackets; packet++)
  {
    cl_uint packet_index = packet * (this->spikePacketSizeWords);
    
    int totalSpikes = this->spikePacketSize;
    
    if(spikeNeuronsPermil != NULL)
    {
      GET_RANDOM_INT(totalSpikes, maxSpikesPerPacket, spikeBufferMinPercentFill, 
        spikeBufferMaxPercentFill);
        
      /* correct spikes if not achieving the target */
      int spikeDifference =  int(targetPacketSpikes*packet - (double)totalNetworkSpikes);
      
      if(spikeDifference < 0)
      {
        totalSpikes = 0;
      }
      else if(spikeDifference > 0)
      {
        totalSpikes = ((cl_uint)spikeDifference < maxSpikesPerPacket)? 
          spikeDifference : maxSpikesPerPacket;
      }
    }
    else
    {
      GET_RANDOM_INT(totalSpikes, this->spikePacketSize, 
        spikeBufferMinPercentFill, spikeBufferMaxPercentFill);
    }
    
    totalNetworkSpikes += totalSpikes;
    
    if(totalSpikes == -1)
    {
      throw SimException("SpikeEvents::setEvents: not able to initialize spikes");
    }
    
    this->dataSpikePacketCounts[packet] = totalSpikes;
    if(totalSpikes <= 0){continue;}
    
    cl_uint packetNeuronsPerSpikeCount = neuronsPerPacket/totalSpikes;
    
    if(packetNeuronsPerSpikeCount > 2)
    {
      cl_uint currentNeuronIdStart = (neuronsPerPacket*packet);
      for(int i = 0; i < totalSpikes; i++)
      {
        cl_uint spiked_neuron = currentNeuronIdStart +
          cl_uint(abs((packetNeuronsPerSpikeCount-1)*((double)rand()/((double)RAND_MAX))));
          
        CL_DATA_TYPE spike_time = CL_DATA_TYPE(abs(CL_DATA_TYPE(this->simulationStepSize)*
          ((double)rand()/((double)RAND_MAX))));
          
        this->dataSpikePackets[packet_index + (this->spikeDatumSize) * i] = spiked_neuron;
          
        *((CL_DATA_TYPE *)(&(this->dataSpikePackets[packet_index + (this->spikeDatumSize)*i + 1]))) 
          = spike_time;
          
        currentNeuronIdStart += packetNeuronsPerSpikeCount;
      }
    }
    else
    {
      for(int i = 0; i < totalSpikes; i++)
      {
        cl_uint spiked_neuron = (neuronsPerPacket*packet) + i;
          
        CL_DATA_TYPE spike_time = CL_DATA_TYPE(abs(CL_DATA_TYPE(this->simulationStepSize)*
          ((double)rand()/((double)RAND_MAX))));

        this->dataSpikePackets[packet_index + (this->spikeDatumSize) * i] = spiked_neuron;
          
        *((CL_DATA_TYPE *)(&(this->dataSpikePackets[packet_index + (this->spikeDatumSize)*i + 1]))) 
          = spike_time;
      }
    }
  }
  
  LOG("SpikeEvents::setEvents: total spikes in the network " << totalNetworkSpikes, 0);

  this->storeBuffers(queue, block, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



cl_uint
SpikeEvents::getSpikeCount
(
  cl::CommandQueue  &queue,
  cl_uint           packet
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
  
  if(packet >= this->spikePackets)
  {
    throw SimException("SpikeEvents::getSpikeCount: packet ID exceeds packet count");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  return this->dataSpikePacketCounts[packet];
}
/**************************************************************************************************/



cl_uint
SpikeEvents::getPastSpikeCount
(
  cl::CommandQueue  &queue,
  cl_uint           packet
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
  
  if(packet >= this->spikePackets)
  {
    throw SimException("SpikeEvents::getPastSpikeCount: packet ID exceeds packet count");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  return this->dataPastSpikePacketCounts[packet];
}
/**************************************************************************************************/



void
SpikeEvents::getSpike
(
  cl::CommandQueue  &queue,
  cl_uint           packet,
  cl_uint           spike,
  cl_uint           &spikedNeuron,
  CL_DATA_TYPE      &spikeTime
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
  
  if(packet >= this->spikePackets)
  {
    throw SimException("SpikeEvents::getSpikeCount: packet ID exceeds packet count");
  }
  
  if(spike >= this->spikePacketSize)
  {
    throw SimException("SpikeEvents::getSpikeCount: spike ID exceeds spike packet boundary");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  cl_uint ptr = (this->spikePacketSizeWords)*packet + (this->spikeDatumSize)*spike;
  spikedNeuron = this->dataSpikePackets[ptr];
  spikeTime = *((CL_DATA_TYPE *)(&(this->dataSpikePackets[ptr + 1])));
}
/**************************************************************************************************/



void
SpikeEvents::getPastSpike
(
  cl::CommandQueue  &queue,
  cl_uint           packet,
  cl_uint           spike,
  cl_uint           &spikedNeuron,
  CL_DATA_TYPE      &spikeTime
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
  
  if(packet >= this->spikePackets)
  {
    throw SimException("SpikeEvents::getPastSpike: packet ID exceeds packet count");
  }
  
  if(spike >= this->spikePacketSize)
  {
    throw SimException("SpikeEvents::getPastSpike: spike ID exceeds spike packet boundary");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  cl_uint ptr = (this->spikePacketSizeWords)*packet + (this->spikeDatumSize)*spike;
  spikedNeuron = this->dataPastSpikePackets[ptr];
  spikeTime = *((CL_DATA_TYPE *)(&(this->dataPastSpikePackets[ptr + 1])));
}
/**************************************************************************************************/



void 
SpikeEvents::invalidateEvents
()
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= (SPIKE_EVENTS_VALID_SPIKES ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
void 
SpikeEvents::invalidateDebug
()
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= (SPIKE_EVENTS_VALID_EXPAND_DEBUG ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
void 
SpikeEvents::invalidateError
()
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= (SPIKE_EVENTS_VALID_EXPAND_ERROR ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



void 
SpikeEvents::refresh
(
  cl::CommandQueue    &queue
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



#if EXPAND_EVENTS_ENABLE
void 
SpikeEvents::expand
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  SynapticEvents                      &synapticEvents,
  Connectome                          &connectome,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  bool expandEventsReset = false;
  cl_uint expandEventsTimeStep = currentTimeStep % timeSlotCount;

#if (EXPAND_EVENTS_DEBUG_ENABLE)
  this->storeBuffers(queue, CL_TRUE, SPIKE_EVENTS_VALID_EXPAND_DEBUG);
#endif

  if(this->setKernelArguments)
  {
    this->setKernelArguments = false;
    
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    SET_KERNEL_ARG_O(this->kernelExpandEvents, this->dataExpandEventsDebugHostBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, this->dataExpandEventsDebugDeviceBuffer, 
      this->argNumExpandEvents++);
#endif

#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG_O(this->kernelExpandEvents, this->dataExpandEventsErrorBuffer, 
      this->argNumExpandEvents++);
#endif

#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
    SET_KERNEL_ARG_O(this->kernelExpandEvents, synapticEvents.dataHistogramBuffer, 
      this->argNumExpandEvents++);
#endif

    SET_KERNEL_ARG_O(this->kernelExpandEvents, this->dataSpikePacketCountsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, this->dataSpikePacketsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, synapticEvents.dataUnsortedEventCountsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, synapticEvents.dataUnsortedEventTargetsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, synapticEvents.dataUnsortedEventDelaysBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, synapticEvents.dataUnsortedEventWeightsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, connectome.dataSynapseTargetsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, connectome.dataSynapseDelaysBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, connectome.dataSynapseWeightsBuffer, 
      this->argNumExpandEvents++);
    SET_KERNEL_ARG_O(this->kernelExpandEvents, connectome.dataSynapsePointerBuffer, 
      this->argNumExpandEvents++);
  }
  
  SET_KERNEL_ARG_O(this->kernelExpandEvents, expandEventsTimeStep, this->argNumExpandEvents);
  
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  double startAppTime = 0, endAppTime = 0;
  if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif

  ENQUEUE_KERNEL(this->kernelExpandEvents, *globalThreadsExpandEvents, *localThreadsExpandEvents, 
    NULL, ndrEvt);
  
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  if(currentTimeStep >= START_PROFILING_AT_STEP)
  {
    cl_uint status = ndrEvt.wait();
    endAppTime = timeStampNs();
    
    if(status != CL_SUCCESS)
    {
      std::stringstream ss;
      ss << "SpikeEvents::expand: Failed cl:Event.wait() due to error code " << status << "\n";
      throw SimException(ss.str());
    }
    REGISTER_TIME(kernelExpandEvents, (endAppTime-startAppTime), 1.0)
  }
#endif

#if DEVICE_HOST_DATA_COHERENCE
  synapticEvents.invalidateEvents();
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  this->invalidateDebug();
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  this->invalidateError();
#endif
#endif

#if (EXPAND_EVENTS_DEBUG_ENABLE)
  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_EXPAND_DEBUG);
#endif

#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
  {
    this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_EXPAND_ERROR);
    
    if(this->dataExpandEventsError[0] != 0)
    {
      std::stringstream ss;
      ss << "SpikeEvents::expand: received error code from the device: " 
        << this->dataExpandEventsError[0] << std::endl;
      throw SimException(ss.str());
    }
  }
#endif

#if EXPAND_EVENTS_VERIFY_ENABLE
  this->verifyExpand
  (
    expandEventsTimeStep, 
    expandEventsReset,
    synapticEvents,
    connectome,
    queue
  );
#endif
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ENABLE
void 
SpikeEvents::expand
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
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
  SynapticEvents                      &synapticEvents,
  Connectome                          &connectome,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  /*Reset the data with new values every resetAtSlot steps for unit test for better represenation*/
  cl_uint expandEventsTimeStep = currentTimeStep % timeSlotCount;
  expandEventsTimeStep = expandEventsTimeStep % resetTimeSlot;
  bool expandEventsReset = (!expandEventsTimeStep);
  
  if(expandEventsReset)
  {
    connectome.setConnections
    (
      queue,
      CL_TRUE,
      gabaPercent,
      minDelay,
      maxDelay
    );
  
    this->setEvents
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
      if(expandEventsTimeStep == 1)
      {
        this->setEvents
        (
          queue,
          CL_TRUE,
          0.0, 0.0, NULL
        );
      }
      else if(expandEventsTimeStep == resetTimeSlot - 1)
      {
        this->setEvents
        (
          queue,
          CL_TRUE,
          100.0, 100.0, NULL
        );
      }
      else
      {
        this->setEvents
        (
          queue,
          CL_TRUE,
          0.0, 
          (double)(100*expandEventsTimeStep/resetTimeSlot), 
          NULL
        );
      }
    }
    else if((testMode >= 0) && (testMode <= 1000))
    {
      this->setEvents
      (
        queue,
        CL_TRUE,
        0.0, 100.0, (double)testMode
      );
    }
    else
    {
      std::stringstream ss;
      ss << "SpikeEvents::expand: invalid test mode: " << testMode << std::endl;
      throw SimException(ss.str());
    }
  }

  this->expand
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    kernelStats,
#endif
    currentTimeStep,
    timeSlotCount,
    synapticEvents,
    connectome,
    queue,
    ndrEvt
  );
}
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ENABLE
void 
SpikeEvents::expand
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  cl_uint                             overWriteSpikesUntilStep,
  double                              spikeBufferMinPercentFill,
  double                              spikeBufferMaxPercentFill,
  double                              spikeNeuronsPermil,
  SynapticEvents                      &synapticEvents,
  Connectome                          &connectome,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

  bool expandEventsReset = false;
  cl_uint expandEventsTimeStep = currentTimeStep % timeSlotCount;

  if(currentTimeStep < overWriteSpikesUntilStep)
  {
    this->setEvents
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
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    kernelStats,
#endif
    currentTimeStep,
    timeSlotCount,
    synapticEvents,
    connectome,
    queue,
    ndrEvt
  );
}
/**************************************************************************************************/
#endif



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_VERIFY_ENABLE)
void 
SpikeEvents::verifyExpand
(
  cl_uint           timeStep,
  bool              reset,
  SynapticEvents    &synapticEvents,
  Connectome        &connectome,
  cl::CommandQueue  &queue
)
/**************************************************************************************************/
{
#if SPIKE_EVENTS_VALIDATION_ENABLE
  this->isInitialized();
#endif

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
        eventCount = synapticEvents.getEventCount(queue, b, s, SynapticEvents::PREVIOUS);
        dataUnsortedEventCountsVerify[eventTimeSlots*b + s] = eventCount;

        for(cl_uint j = 0; j < eventCount; j++)
        {
          synapticEvents.getEvent(queue, b, j, s, SynapticEvents::PREVIOUS, 
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
  for(cl_uint packet = 0; packet < this->spikePackets; packet++)
  {
    cl_uint total_spikes = this->getSpikeCount(queue, packet);

    /*Iterate through spikes in a current packet*/
    for(cl_uint i = 0; i < total_spikes; i++)
    {
      cl_uint spiked_neuron = 0;
      cl_float spike_time = 0;
      this->getSpike(queue, packet, i, spiked_neuron, spike_time);

      cl_uint synapse_count = connectome.getSynapseCount(queue, spiked_neuron);
      
      /*Iterate through synapses of spiked neuron*/
      for(cl_uint j = 0; j < synapse_count; j++)
      {
#if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
        cl_uint buffer = (packet/(EXPAND_EVENTS_SPIKE_PACKETS_PER_WF));
#else
        cl_uint buffer = (packet/(EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF));
#endif
        cl_uint target_neuron = 0; cl_float weight = 0; cl_float delay = 0;
        
        connectome.getSynapse(queue, spiked_neuron, j, target_neuron, delay, weight);
        
        /*Add delay to spike time, decrement*/
        cl_float event_time = spike_time + delay; event_time = event_time - 1.0f;
        /*Make it relative to its time slot*/
        cl_float event_time_binned = event_time - (int)event_time;
        /*Events with 0.0 time bounce back to the previous time slot*/
        if(event_time_binned == 0.0f){event_time_binned = 1.0f;}
        cl_uint bin_correction = (int)event_time_binned;
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
        cl_uint bin = (event_time_uint>>EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT) & 
          EXPAND_EVENTS_HISTOGRAM_BIN_MASK;
          
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
    for(cl_uint packet = 0; packet < this->spikePackets; packet++)
    {
      for(cl_uint time_slot = 0; time_slot < eventTimeSlots; time_slot++)
      {
#if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
        cl_uint buffer = (packet/(EXPAND_EVENTS_SPIKE_PACKETS_PER_WF));
#else
        cl_uint buffer = (packet/(EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF));
#endif
        if(dataUnsortedEventCountsVerify[eventTimeSlots * buffer + time_slot] >= 
          eventBufferSize)
        {
          dataUnsortedEventCountsVerify[eventTimeSlots * buffer + time_slot] = 
            eventBufferSize-1;
        }
      }
    }
    
    ss << "SpikeEvents::verifyExpand: Event buffer overflow. " << 
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
        cl_uint eventCount = synapticEvents.getEventCount(queue, i, s, SynapticEvents::RECENT);
        
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
        ss << "\nSpikeEvents::verifyExpand: Failed to match synaptic event time slot counters " 
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
          synapticEvents.getEvent(queue, p, j, s, SynapticEvents::RECENT, 
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
          if(timeA < 0.0f || timeA > (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY))
          {
            ss << "SpikeEvents::verifyExpand: Found an event in actual data with time outside "
              << "of valid range: value " << timeA << ", range " << 0.0f
              << " - " << (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY) << std::endl;
            result = false; break;
          }
          if(timeV< 0.0f || timeV > (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY))
          {
            ss << "SpikeEvents::verifyExpand: Found an event in verification data with delay "
              << "outside of valid range: value " << timeV << ", range " << 0.0f
              << " - " << (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY) << std::endl;
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
          ss << "SpikeEvents::verifyExpand: Failed to verify sunaptic data checksum(s), time slot " 
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
            (synapticEvents.getHistogramItem(queue, s, j, k, SynapticEvents::RECENT) != 
            dataHistogramVerify[offset]);
        }
      }
      
      if(error_event_totals)
      {
        ss << "SpikeEvents::verifyExpand: Failed to match histogram " << error_event_totals 
          << " times." << std::endl;
        result = false; break;
      }
    }
  }
#endif
  }
  CATCH(ss, SpikeEvents::verifyExpand, result = false;)
  
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
#endif
/**************************************************************************************************/



void
SpikeEvents::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  if
  (
    (selectBitMask & SPIKE_EVENTS_VALID_SPIKES) && 
    !((this->dataValid) & SPIKE_EVENTS_VALID_SPIKES)
  )
  {
    swap2(dataSpikePackets, dataPastSpikePackets);
    swap2(dataSpikePacketCounts, dataPastSpikePacketCounts);
  
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSpikePacketCountsBuffer, 
      this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSpikePacketsBuffer, 
      this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
      
    this->dataValid |= SPIKE_EVENTS_VALID_SPIKES;
  }
  
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
  if
  (
    (selectBitMask & SPIKE_EVENTS_VALID_EXPAND_DEBUG) && 
    !((this->dataValid) & SPIKE_EVENTS_VALID_EXPAND_DEBUG)
  )
  {

    this->dataValid |= SPIKE_EVENTS_VALID_EXPAND_DEBUG;
  }
#endif

#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  if
  (
    (selectBitMask & SPIKE_EVENTS_VALID_EXPAND_ERROR) && 
    !((this->dataValid) & SPIKE_EVENTS_VALID_EXPAND_ERROR)
  )
  {

    this->dataValid |= SPIKE_EVENTS_VALID_EXPAND_ERROR;
  }
#endif
}
/**************************************************************************************************/



void
SpikeEvents::reset(bool checkForNull)
/**************************************************************************************************/
{
  try
  {
  this->resetObject = true;
  this->dataValid = 0;
  this->srandSeed = 0;
  this->srandCounter = 0;
  this->neuronCount = 0;
  this->spikePacketSize = 0;
  this->spikePacketSizeWords = 0;
  this->spikePackets = 0;
  this->simulationStepSize = 0;
  this->spikeDatumSize = 0;
  this->dataToSimulationLogFile = NULL;
  this->dataToReportLogFile = NULL;
    
  this->dataSpikePacketsSize = 0;
  this->dataSpikePacketsSizeBytes = 0;

  this->dataSpikePacketCountsSize = 0;
  this->dataSpikePacketCountsSizeBytes = 0;
  
#if EXPAND_EVENTS_ENABLE
  setKernelArguments = true;
#endif

  if(checkForNull)
  {
    if(this->dataSpikePacketCounts)
    {
      free(this->dataSpikePacketCounts);
      this->dataSpikePacketCounts = NULL;
    }
    if(this->dataSpikePackets)
    {
      free(this->dataSpikePackets);
      this->dataSpikePackets = NULL;
    }
    if(this->dataPastSpikePacketCounts)
    {
      free(this->dataPastSpikePacketCounts);
      this->dataPastSpikePacketCounts = NULL;
    }
    if(this->dataPastSpikePackets)
    {
      free(this->dataPastSpikePackets);
      this->dataPastSpikePackets = NULL;
    }
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
    if(this->dataExpandEventsDebugHost)
    {
      free(this->dataExpandEventsDebugHost);
      this->dataExpandEventsDebugHost = NULL;
    }
    if(this->dataExpandEventsDebugDevice)
    {
      free(this->dataExpandEventsDebugDevice);
      this->dataExpandEventsDebugDevice = NULL;
    }
#endif
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    if(this->dataExpandEventsError)
    {
      free(this->dataExpandEventsError);
      this->dataExpandEventsError = NULL;
    }
#endif 
  }
  else
  {
    this->dataSpikePacketCounts = NULL;
    this->dataSpikePackets = NULL;
    this->dataPastSpikePacketCounts = NULL;
    this->dataPastSpikePackets = NULL;
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
    this->dataExpandEventsDebugHost = NULL;
    this->dataExpandEventsDebugDevice = NULL;
#endif
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    this->dataExpandEventsError = NULL;
#endif
  }
  }
  CATCH(std::cerr, SpikeEvents::reset, throw SimException("SpikeEvents::reset: failed.");)
}
/**************************************************************************************************/



void
SpikeEvents::storeBuffers
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           selectBitMask
)
/**************************************************************************************************/
{
  if(selectBitMask & SPIKE_EVENTS_VALID_SPIKES)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSpikePacketsBuffer, 
      this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSpikePacketCountsBuffer, 
      this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
      
    this->dataValid |= SPIKE_EVENTS_VALID_SPIKES;
  }
  
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
  if(selectBitMask & SPIKE_EVENTS_VALID_EXPAND_DEBUG)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataExpandEventsDebugHostBuffer, 
      this->dataExpandEventsDebugHostSizeBytes, this->dataExpandEventsDebugHost);
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataExpandEventsDebugDeviceBuffer, 
      this->dataExpandEventsDebugDeviceSizeBytes, this->dataExpandEventsDebugDevice);
    
    this->dataValid |= SPIKE_EVENTS_VALID_EXPAND_DEBUG;
  }
#endif

#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  if(selectBitMask & SPIKE_EVENTS_VALID_EXPAND_ERROR)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataExpandEventsErrorBuffer, 
      this->dataExpandEventsErrorSizeBytes, this->dataExpandEventsError);
    
    this->dataValid |= SPIKE_EVENTS_VALID_EXPAND_ERROR;
  }
#endif
}
/**************************************************************************************************/



void
SpikeEvents::isInitialized
()
/**************************************************************************************************/
{
  if(this->resetObject)
  {
    throw SimException("SpikeEvents::isInitialized: the object was not initialized");
  }
}
/**************************************************************************************************/
