
/* ===============================================================================================



  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Data_Connectome.hpp"

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
Data_Connectome::resetConnections
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  double              gabaPercent,
  double              minDelay,
  double              maxDelay
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(minDelay > maxDelay)
  {
    THROW_SIMEX("Data_Connectome::resetConnections: parameter minDelay " << minDelay << " > " 
      << maxDelay);
  }
#endif

  unsigned long long gabaSynapseCount = 0;
  unsigned long long ampaSynapseCount = 0;
  unsigned long long maxSynapsePerNeuron = 0;
  unsigned long long minSynapsePerNeuron = 0xFFFFFFFF;
  double averageSynapsePerNeuron = 0;
  
  memset(this->dataSynapseTargets, 0, this->dataSynapseTargetsSizeBytes);
  memset(this->dataSynapseDelays, 0, this->dataSynapseDelaysSizeBytes);
  memset(this->dataSynapseWeights, 0, this->dataSynapseWeightsSizeBytes);

  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG_SIM("Data_Connectome::resetConnections: set srand seed to " << this->srandSeed);
  
  /*init synaptic data*/
  for(cl_uint i = 0; i < this->neuronCount; i++)
  {
    cl_uint ptrStart = this->dataSynapsePointer[i];
    cl_uint ptrEnd = this->dataSynapsePointer[i+1];
    
#if CLASS_VALIDATION_ENABLE
    if(ptrStart > ptrEnd)
    {
      THROW_SIMEX("Data_Connectome::resetConnections: negative connection count (" 
        << ptrEnd - ptrStart << ") detected for neuron ID " << i);
    }
#endif

    cl_uint synapseCount = ptrEnd - ptrStart;
    
    maxSynapsePerNeuron = (maxSynapsePerNeuron > synapseCount) ? maxSynapsePerNeuron : synapseCount;
    minSynapsePerNeuron = (minSynapsePerNeuron < synapseCount) ? minSynapsePerNeuron : synapseCount;
    averageSynapsePerNeuron = averageSynapsePerNeuron + 
      (synapseCount - averageSynapsePerNeuron)/(i+1);

    /*Initialize synaptic structure*/
    for(cl_uint j = 0; j < synapseCount; j++)
    {
      cl_uint offset = (ptrStart + j);
      
#if CLASS_VALIDATION_ENABLE
      if(offset > this->dataSynapseTargetsSize)
      {
        THROW_SIMEX("Data_Connectome::resetConnections: synapse pointer is outside of "
          << "synapse data range for neuron ID " << i << ": " << offset << " > " 
          << this->dataSynapseTargetsSize);
      }
#endif
      
      cl_uint reinitCount = 100;
      this->dataSynapseTargets[offset] = i;
      
      /*target neuron (avoid direct feedback connections)*/
      while(this->dataSynapseTargets[offset] == i && reinitCount)
      {
        reinitCount--;
        this->dataSynapseTargets[offset] = 
          cl_uint(abs((this->neuronCount-1)*((double)rand()/((double)RAND_MAX))));
      }
      
#if CLASS_VALIDATION_ENABLE
      if(!reinitCount)
      {
        THROW_SIMEX("Data_Connectome::resetConnections: failed to generate a connection "
          << "without direct feedback for neuron ID " << i);
      }
#endif
      
      /*weight*/
      double weightType = abs(100.0*((double)rand()/((double)RAND_MAX)));
      cl_float weight = 6.0f/1.4f;
      if(weightType < gabaPercent)
      {
        weight = -67.0f/1.4f;
        gabaSynapseCount++;
      }
      else
      {
        ampaSynapseCount++;
      }
      this->dataSynapseWeights[offset] = cl_float(weight*((double)rand()/((double)RAND_MAX)));
      
      /*delay*/
      this->dataSynapseDelays[offset] = cl_float(minDelay + 
        abs((maxDelay-minDelay)*((double)rand()/((double)RAND_MAX))));
    }
  }

  LOG_REP("Gaba Synapse Count:" << gabaSynapseCount);
  LOG_REP("Ampa Synapse Count:" << ampaSynapseCount);
  LOG_REP("Total Synapse Count:" << (ampaSynapseCount+gabaSynapseCount));
  LOG_REP("Max Synapses/Neuron:" << maxSynapsePerNeuron);
  LOG_REP("Min Synapses/Neuron:" << minSynapsePerNeuron);
  LOG_REP("Average Synapses/Neuron:" << averageSynapsePerNeuron);

  this->storeData(queue, block);
}
/**************************************************************************************************/



cl_uint
Data_Connectome::getSynapseCount
(
  cl_uint           neuronID
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(neuronID >= this->neuronCount)
  {
    throw SimException("Data_Connectome::getSynapseCount: neuron ID exceeds neuron count");
  }
#endif
  /* Disable temporarely since device doesn't change connectome w/o STDP
  this->getData(queue, CL_TRUE);*/

  cl_uint ptrStart = this->dataSynapsePointer[neuronID];
  cl_uint ptrEnd = this->dataSynapsePointer[neuronID+1];
  
#if CLASS_VALIDATION_ENABLE
  if(ptrStart > ptrEnd)
  {
    THROW_SIMEX("Data_Connectome::getSynapseCount: negative connection count (" << ptrEnd - ptrStart 
      << ") detected for neuron ID " << neuronID);
  }
#endif

  return ptrEnd - ptrStart;
}
/**************************************************************************************************/



void
Data_Connectome::getSynapse
(
  cl_uint           neuronID,
  cl_uint           synapse,
  cl_uint           &targetNeuron,
  CL_DATA_TYPE      &delay,
  CL_DATA_TYPE      &weight
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(neuronID >= this->neuronCount)
  {
    throw SimException("Data_Connectome::getSynapse: neuron ID exceeds neuron count");
  }
#endif
  /* Disable temporarely since device doesn't change connectome w/o STDP
  this->getData(queue, CL_TRUE);*/
  
  cl_uint ptrStart = this->dataSynapsePointer[neuronID];

#if CLASS_VALIDATION_ENABLE
  cl_uint ptrEnd = this->dataSynapsePointer[neuronID+1];
  
  if(ptrStart > ptrEnd)
  {
    THROW_SIMEX("Data_Connectome::getSynapse: negative connection count (" << ptrEnd - ptrStart 
      << ") detected for neuron ID " << neuronID);
  }
  
  cl_uint synapseCount = ptrEnd - ptrStart;
  
  if(synapse >= synapseCount)
  {
    THROW_SIMEX("Data_Connectome::getSynapse: synapse ID " << synapse << " exceeds synapse count " 
      << synapseCount << "for neuron ID " << neuronID);
  }
#endif
  
  cl_uint ptr = ptrStart + synapse;
  targetNeuron = this->dataSynapseTargets[ptr];
  delay = this->dataSynapseDelays[ptr];
  weight = this->dataSynapseWeights[ptr];
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Data_Connectome::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             &kernelStats,
  cl_uint                             maxConnectionsPerNeuron,
  double                              connectionDeviationRatio
)
/**************************************************************************************************/
{
  /* allocate memory for synaptic pointer */
  CALLOC(dataSynapsePointer, cl_uint, (this->neuronCount + 1));
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSynapsePointer, kernelStats);
  
  /*init synaptic pointer*/
  this->dataSynapsePointer[0] = 0;
  for(cl_uint i = 1; i < (this->neuronCount + 1); i++)
  {
    this->dataSynapsePointer[i] = this->dataSynapsePointer[i-1] + 
      cl_uint(((double)maxConnectionsPerNeuron)*((1.0-connectionDeviationRatio)+
      abs((connectionDeviationRatio*((double)rand()/((double)RAND_MAX))))));
  }

  /* allocate memory for synaptic data */
  size_t size = this->dataSynapsePointer[this->neuronCount];
    
  CALLOC(dataSynapseTargets, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSynapseTargets, kernelStats);
  
  CALLOC(dataSynapseDelays, cl_float, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSynapseDelays, kernelStats);
  
  CALLOC(dataSynapseWeights, cl_float, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSynapseWeights, kernelStats);
  
  /*create device buffers*/
  CREATE_BUFFER_O(context, CL_MEM_READ_ONLY, this->dataSynapseTargetsBuffer, 
    this->dataSynapseTargetsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_ONLY, this->dataSynapseDelaysBuffer, 
    this->dataSynapseDelaysSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_ONLY, this->dataSynapseWeightsBuffer, 
    this->dataSynapseWeightsSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_ONLY, this->dataSynapsePointerBuffer, 
    this->dataSynapsePointerSizeBytes);
    
  this->storeData(queue, block);
}
/**************************************************************************************************/



void
Data_Connectome::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  if(!this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataSynapsePointerBuffer, 
      this->dataSynapsePointerSizeBytes, this->dataSynapsePointer);
    ENQUEUE_READ_BUFFER(block, queue, this->dataSynapseTargetsBuffer, 
      this->dataSynapseTargetsSizeBytes, this->dataSynapseTargets);
    ENQUEUE_READ_BUFFER(block, queue, this->dataSynapseDelaysBuffer, 
      this->dataSynapseDelaysSizeBytes, this->dataSynapseDelays);
    ENQUEUE_READ_BUFFER(block, queue, this->dataSynapseWeightsBuffer, 
      this->dataSynapseWeightsSizeBytes, this->dataSynapseWeights);

    this->dataValid = true;
  }
}
/**************************************************************************************************/



void
Data_Connectome::storeData
(
  cl::CommandQueue  &queue,
  cl_bool           block
)
/**************************************************************************************************/
{
  ENQUEUE_WRITE_BUFFER(block, queue, this->dataSynapsePointerBuffer, 
    this->dataSynapsePointerSizeBytes, this->dataSynapsePointer);
  ENQUEUE_WRITE_BUFFER(block, queue, this->dataSynapseTargetsBuffer, 
    this->dataSynapseTargetsSizeBytes, this->dataSynapseTargets);
  ENQUEUE_WRITE_BUFFER(block, queue, this->dataSynapseDelaysBuffer, 
    this->dataSynapseDelaysSizeBytes, this->dataSynapseDelays);
  ENQUEUE_WRITE_BUFFER(block, queue, this->dataSynapseWeightsBuffer, 
    this->dataSynapseWeightsSizeBytes, this->dataSynapseWeights);
    
  this->dataValid = true;
}
/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
