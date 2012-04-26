
/* ===============================================================================================



  =============================================================================================== */



#include "SpikeEvents.hpp"



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
  cl_uint                             simulationTimeSteps,
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
  this->dataValid = true;
  this->neuronCount = neuronCount;
  this->simulationTimeSteps = simulationTimeSteps;
  this->spikeDatumSize = spikeDatumSize;
  this->spikePackets = spikePackets;
  this->spikePacketSize = spikePacketSize;
  this->spikePacketSizeWords = spikePacketSizeWords;
  this->dataToSimulationLogFile = dataToSimulationLogFile;
  this->dataToReportLogFile = dataToReportLogFile;

  CALLOC_O(dataSpikePackets, cl_uint, (this->spikePackets)*(this->spikePacketSizeWords));
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSpikePackets, kernelStats);
  CALLOC_O(dataSpikePacketCounts, cl_uint, this->spikePackets);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSpikePacketCounts, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSpikePacketsBuffer, 
    this->dataSpikePacketsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSpikePacketCountsBuffer, 
    this->dataSpikePacketCountsSizeBytes);
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
  
  memset(this->dataSpikePackets, 0, this->dataSpikePacketsSizeBytes);
  memset(this->dataSpikePacketCounts, 0, this->dataSpikePacketCountsSizeBytes);
  
  this->storeBuffers(queue, block);
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
          
        CL_DATA_TYPE spike_time = CL_DATA_TYPE(abs(CL_DATA_TYPE(this->simulationTimeSteps)*
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
          
        CL_DATA_TYPE spike_time = CL_DATA_TYPE(abs(CL_DATA_TYPE(this->simulationTimeSteps)*
          ((double)rand()/((double)RAND_MAX))));

        this->dataSpikePackets[packet_index + (this->spikeDatumSize) * i] = spiked_neuron;
          
        *((CL_DATA_TYPE *)(&(this->dataSpikePackets[packet_index + (this->spikeDatumSize)*i + 1]))) 
          = spike_time;
      }
    }
  }
  
  LOG("SpikeEvents::setEvents: total spikes in the network " << totalNetworkSpikes, 0);

  this->storeBuffers(queue, block);
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

  getSpikeEvents(queue, CL_TRUE);
  
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

  this->getPastSpikeEvents(queue, CL_TRUE);
  
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

  getSpikeEvents(queue, CL_TRUE);
  
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

  this->getPastSpikeEvents(queue, CL_TRUE);
  
  cl_uint ptr = (this->spikePacketSizeWords)*packet + (this->spikeDatumSize)*spike;
  spikedNeuron = this->dataPastSpikePackets[ptr];
  spikeTime = *((CL_DATA_TYPE *)(&(this->dataPastSpikePackets[ptr + 1])));
}
/**************************************************************************************************/



void
SpikeEvents::deletePastEvents
()
/**************************************************************************************************/
{
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

  this->dataValid = false;
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
SpikeEvents::getSpikeEvents
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  if(!this->dataValid)
  {
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSpikePacketCountsBuffer, 
      this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSpikePacketsBuffer, 
      this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
    this->dataValid = true;
  }
}
/**************************************************************************************************/



void
SpikeEvents::getPastSpikeEvents
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  if(this->dataPastSpikePacketCounts == NULL)
  {
    this->dataPastSpikePacketCounts = (cl_uint *)calloc(this->spikePackets, sizeof(cl_uint));
  }

  memcpy(this->dataPastSpikePacketCounts, this->dataSpikePacketCounts, ((this->spikePackets)*
    sizeof(cl_uint)));
    
  if(this->dataPastSpikePackets == NULL)
  {
    this->dataPastSpikePackets = (cl_uint *)calloc((this->spikePackets)*
    (this->spikePacketSizeWords), sizeof(cl_uint));
  }

  memcpy(this->dataPastSpikePackets, this->dataSpikePackets, ((this->spikePackets)*
    (this->spikePacketSizeWords)*sizeof(cl_uint)));
    
  if(!this->dataValid)
  {
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSpikePacketCountsBuffer, 
      this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSpikePacketsBuffer, 
      this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
    this->dataValid = true;
  }
}
/**************************************************************************************************/



void
SpikeEvents::reset(bool checkForNull)
/**************************************************************************************************/
{
  this->resetObject = true;
  this->dataValid = false;
  this->srandSeed = 0;
  this->srandCounter = 0;
  this->neuronCount = 0;
  this->spikePacketSize = 0;
  this->spikePacketSizeWords = 0;
  this->spikePackets = 0;
  this->simulationTimeSteps = 0;
  this->spikeDatumSize = 0;
  this->dataToSimulationLogFile = NULL;
  this->dataToReportLogFile = NULL;
    
  this->dataSpikePacketsSize = 0;
  this->dataSpikePacketsSizeBytes = 0;

  this->dataSpikePacketCountsSize = 0;
  this->dataSpikePacketCountsSizeBytes = 0;
  
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
  }
  else
  {
    this->dataSpikePacketCounts = NULL;
    this->dataSpikePackets = NULL;
    this->dataPastSpikePacketCounts = NULL;
    this->dataPastSpikePackets = NULL;
  }
}
/**************************************************************************************************/



void
SpikeEvents::storeBuffers
(
  cl::CommandQueue  &queue,
  cl_bool           block
)
/**************************************************************************************************/
{
  ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSpikePacketsBuffer, 
    this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
  ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSpikePacketCountsBuffer, 
    this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
    
  this->dataValid = true;
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
