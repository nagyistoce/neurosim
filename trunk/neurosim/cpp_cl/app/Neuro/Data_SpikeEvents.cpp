
/* ===============================================================================================



  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Data_SpikeEvents.hpp"

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
Data_SpikeEvents::clearEvents
(
  cl::CommandQueue  &queue,
  cl_bool           block
)
/**************************************************************************************************/
{
  swap2(this->dataSpikePackets, this->dataPastSpikePackets);
  swap2(this->dataSpikePacketCounts, this->dataPastSpikePacketCounts);
  
  memset(this->dataSpikePackets, 0, this->dataSpikePacketsSizeBytes);
  memset(this->dataSpikePacketCounts, 0, this->dataSpikePacketCountsSizeBytes);
  
  this->storeData(queue, block, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



void 
Data_SpikeEvents::setEvents
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  double              spikeBufferMinPercentFill,
  double              spikeBufferMaxPercentFill,
  double              spikeNeuronsPermil
)
/**************************************************************************************************/
{
  swap2(this->dataSpikePackets, this->dataPastSpikePackets);
  swap2(this->dataSpikePacketCounts, this->dataPastSpikePacketCounts);
  
  memset(this->dataSpikePackets, 0, this->dataSpikePacketsSizeBytes);
  memset(this->dataSpikePacketCounts, 0, this->dataSpikePacketCountsSizeBytes);
  
  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG_SIM("Data_SpikeEvents::setEvents: set srand seed to " << this->srandSeed);

  cl_uint neuronsPerPacket = (this->neuronCount)/(this->spikePackets);
  cl_uint maxSpikesPerPacket = neuronsPerPacket;
  double targetPacketSpikes = 0;
  cl_uint totalNetworkSpikes = 0;

  if(spikeNeuronsPermil >= 0.0)
  {
    maxSpikesPerPacket = cl_uint(ceil(1.5*(double)maxSpikesPerPacket*(spikeNeuronsPermil/1000)));
    targetPacketSpikes = (double)neuronsPerPacket*(spikeNeuronsPermil/1000);
    if(maxSpikesPerPacket < 1){maxSpikesPerPacket = 1;}
    
#if CLASS_VALIDATION_ENABLE
    if(maxSpikesPerPacket > this->spikePacketSize)
    {
      THROW_SIMEX("Data_SpikeEvents::setEvents: can't fit spikes within " << this->spikePacketSize 
        << " boundary");
    }
#endif
  }

  LOG_SIM("Data_SpikeEvents::setEvents: max spikes per packet is " << maxSpikesPerPacket);

  /* init spike data */
  for(cl_uint packet = 0; packet < this->spikePackets; packet++)
  {
    cl_uint packet_index = packet * (this->spikePacketSizeWords);
    
    cl_uint totalSpikes = this->spikePacketSize;
    
    if(spikeNeuronsPermil >= 0.0)
    {
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
      GET_RANDOM_INT(totalSpikes, maxSpikesPerPacket, spikeBufferMinPercentFill, 
        spikeBufferMaxPercentFill);
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
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
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
      GET_RANDOM_INT(totalSpikes, this->spikePacketSize, 
        spikeBufferMinPercentFill, spikeBufferMaxPercentFill);
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
    }
    
    totalNetworkSpikes += totalSpikes;
    
    this->dataSpikePacketCounts[packet] = totalSpikes;
    if(totalSpikes <= 0){continue;}
    
    cl_uint packetNeuronsPerSpikeCount = neuronsPerPacket/totalSpikes;
    
    if(packetNeuronsPerSpikeCount > 2)
    {
      cl_uint currentNeuronIdStart = (neuronsPerPacket*packet);
      for(cl_uint i = 0; i < totalSpikes; i++)
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
      for(cl_uint i = 0; i < totalSpikes; i++)
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
  
  LOG_SIM("Data_SpikeEvents::setEvents: total spikes in the network " << totalNetworkSpikes);

  this->storeData(queue, block, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



cl_uint
Data_SpikeEvents::getSpikePacketCount
()
/**************************************************************************************************/
{
  return this->spikePackets;
}
/**************************************************************************************************/



cl_uint
Data_SpikeEvents::getSpikeCount
(
  cl::CommandQueue  &queue,
  cl_uint           packet
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(packet >= this->spikePackets)
  {
    throw SimException("Data_SpikeEvents::getSpikeCount: packet ID exceeds packet count");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  return this->dataSpikePacketCounts[packet];
}
/**************************************************************************************************/



cl_uint
Data_SpikeEvents::getPastSpikeCount
(
  cl::CommandQueue  &queue,
  cl_uint           packet
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(packet >= this->spikePackets)
  {
    throw SimException("Data_SpikeEvents::getPastSpikeCount: packet ID exceeds packet count");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  return this->dataPastSpikePacketCounts[packet];
}
/**************************************************************************************************/



void
Data_SpikeEvents::getSpike
(
  cl::CommandQueue  &queue,
  cl_uint           packet,
  cl_uint           spike,
  cl_uint           &spikedNeuron,
  CL_DATA_TYPE      &spikeTime
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(packet >= this->spikePackets)
  {
    throw SimException("Data_SpikeEvents::getSpikeCount: packet ID exceeds packet count");
  }
  
  if(spike >= this->spikePacketSize)
  {
    throw SimException("Data_SpikeEvents::getSpikeCount: spike ID exceeds spike packet boundary");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  cl_uint ptr = (this->spikePacketSizeWords)*packet + (this->spikeDatumSize)*spike;
  spikedNeuron = this->dataSpikePackets[ptr];
  spikeTime = *((CL_DATA_TYPE *)(&(this->dataSpikePackets[ptr + 1])));
}
/**************************************************************************************************/



void
Data_SpikeEvents::getPastSpike
(
  cl::CommandQueue  &queue,
  cl_uint           packet,
  cl_uint           spike,
  cl_uint           &spikedNeuron,
  CL_DATA_TYPE      &spikeTime
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(packet >= this->spikePackets)
  {
    throw SimException("Data_SpikeEvents::getPastSpike: packet ID exceeds packet count");
  }
  
  if(spike >= this->spikePacketSize)
  {
    throw SimException("Data_SpikeEvents::getPastSpike: spike ID exceeds spike packet boundary");
  }
#endif

  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
  
  cl_uint ptr = (this->spikePacketSizeWords)*packet + (this->spikeDatumSize)*spike;
  spikedNeuron = this->dataPastSpikePackets[ptr];
  spikeTime = *((CL_DATA_TYPE *)(&(this->dataPastSpikePackets[ptr + 1])));
}
/**************************************************************************************************/



void 
Data_SpikeEvents::invalidateEvents
()
/**************************************************************************************************/
{
  this->dataValid &= (SPIKE_EVENTS_VALID_SPIKES ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
Data_SpikeEvents::refresh
(
  cl::CommandQueue    &queue
)
/**************************************************************************************************/
{
  this->getData(queue, CL_TRUE, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Data_SpikeEvents::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
  this->dataPastSpikePacketCounts = (cl_uint *)calloc(this->spikePackets, sizeof(cl_uint));
  this->dataPastSpikePackets = (cl_uint *)calloc((this->spikePackets)*(this->spikePacketSizeWords), 
    sizeof(cl_uint));

  CALLOC(dataSpikePackets, cl_uint, (this->spikePackets)*(this->spikePacketSizeWords));
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSpikePackets, kernelStats);
  CALLOC(dataSpikePacketCounts, cl_uint, this->spikePackets);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSpikePacketCounts, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSpikePacketsBuffer, 
    this->dataSpikePacketsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataSpikePacketCountsBuffer, 
    this->dataSpikePacketCountsSizeBytes);
    
  this->storeData(queue, block, SPIKE_EVENTS_VALID_SPIKES);
}
/**************************************************************************************************/



void
Data_SpikeEvents::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  IF_HIT_READ(selectBitMask, SPIKE_EVENTS_VALID_SPIKES, this->dataValid)
  {
    swap2(this->dataSpikePackets, this->dataPastSpikePackets);
    swap2(this->dataSpikePacketCounts, this->dataPastSpikePacketCounts);
  
    ENQUEUE_READ_BUFFER(block, queue, this->dataSpikePacketCountsBuffer, 
      this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
    ENQUEUE_READ_BUFFER(block, queue, this->dataSpikePacketsBuffer, 
      this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
      
    this->dataValid |= SPIKE_EVENTS_VALID_SPIKES;
  }
}
/**************************************************************************************************/



void
Data_SpikeEvents::storeData
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           selectBitMask
)
/**************************************************************************************************/
{
  if(selectBitMask & SPIKE_EVENTS_VALID_SPIKES)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataSpikePacketsBuffer, 
      this->dataSpikePacketsSizeBytes, this->dataSpikePackets);
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataSpikePacketCountsBuffer, 
      this->dataSpikePacketCountsSizeBytes, this->dataSpikePacketCounts);
      
    this->dataValid |= SPIKE_EVENTS_VALID_SPIKES;
  }
}
/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
