
/* ===============================================================================================



  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Connectome.hpp"

/**************************************************************************************************/



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



Connectome::Connectome()
/**************************************************************************************************/
{
  this->reset(false);
}
/**************************************************************************************************/



Connectome::~Connectome()
/**************************************************************************************************/
{
  this->reset(true);
}
/**************************************************************************************************/



void
Connectome::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  cl_uint                             neuronCount,
  cl_uint                             maxConnectionsPerNeuron,
  double                              connectionDeviationRatio,
  struct kernelStatistics             *kernelStats,
  std::stringstream                   *dataToSimulationLogFile,
  std::stringstream                   *dataToReportLogFile
)
/**************************************************************************************************/
{
  if(!this->resetObject)
  {
    throw SimException("Connectome::initialize: attemp to initialize without reset");
  }
  
  this->resetObject = false;
  this->dataValid = false;
  this->neuronCount = neuronCount;
  this->dataToSimulationLogFile = dataToSimulationLogFile;
  this->dataToReportLogFile = dataToReportLogFile;

  /* allocate memory for synaptic pointer */
  CALLOC_O(dataSynapsePointer, cl_uint, (this->neuronCount + 1));
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
    
  CALLOC_O(dataSynapseTargets, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSynapseTargets, kernelStats);
  
  CALLOC_O(dataSynapseDelays, cl_float, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataSynapseDelays, kernelStats);
  
  CALLOC_O(dataSynapseWeights, cl_float, size);
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
    
  this->storeBuffers(queue, block);
}
/**************************************************************************************************/



void 
Connectome::setConnections
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  double              gabaPercent,
  double              minDelay,
  double              maxDelay
)
/**************************************************************************************************/
{
  this->isInitialized();
  
  if(minDelay > maxDelay)
  {
    std::stringstream ss;
    ss << "Connectome::setEvents: parameter minDelay " << minDelay << " > " << maxDelay 
      << std::endl;
    throw SimException(ss.str());
  }

  unsigned long long gabaSynapseCount = 0;
  unsigned long long ampaSynapseCount = 0;
  unsigned long long maxSynapsePerNeuron = 0;
  unsigned long long minSynapsePerNeuron = 0xFFFFFFFF;
  double averageSynapsePerNeuron = 0;
  
  memset(this->dataSynapseTargets, 0, this->dataSynapseTargetsSizeBytes);
  memset(this->dataSynapseDelays, 0, this->dataSynapseDelaysSizeBytes);
  memset(this->dataSynapseWeights, 0, this->dataSynapseWeightsSizeBytes);

  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG("Connectome::setEvents: set srand seed to " << this->srandSeed, 0);
  
  /*init synaptic data*/
  for(cl_uint i = 0; i < this->neuronCount; i++)
  {
    cl_uint ptrStart = this->dataSynapsePointer[i];
    cl_uint ptrEnd = this->dataSynapsePointer[i+1];
    
    if(ptrStart > ptrEnd)
    {
      std::stringstream ss;
      ss << "Connectome::setEvents: negative connection count (" << ptrEnd - ptrStart 
        << ") detected for neuron ID " << i << std::endl;
      throw SimException(ss.str());
    }
    
    cl_uint synapseCount = ptrEnd - ptrStart;
    
    maxSynapsePerNeuron = max(maxSynapsePerNeuron,synapseCount);
    minSynapsePerNeuron = min(minSynapsePerNeuron,synapseCount);
    averageSynapsePerNeuron = averageSynapsePerNeuron + 
      (synapseCount - averageSynapsePerNeuron)/(i+1);

    /*Initialize synaptic structure*/
    for(cl_uint j = 0; j < synapseCount; j++)
    {
      cl_uint offset = (ptrStart + j);
      
      if(offset > this->dataSynapseTargetsSize)
      {
        std::stringstream ss;
        ss << "Connectome::setEvents: synapse pointer is outside of "
          << "synapse data range for neuron ID " << i << ": " << offset << " > " 
          << this->dataSynapseTargetsSize << std::endl;
        throw SimException(ss.str());
      }
      
      cl_uint reinitCount = 100;
      this->dataSynapseTargets[offset] = i;
      
      /*target neuron (avoid direct feedback connections)*/
      while(this->dataSynapseTargets[offset] == i && reinitCount)
      {
        reinitCount--;
        this->dataSynapseTargets[offset] = 
          cl_uint(abs((this->neuronCount-1)*((double)rand()/((double)RAND_MAX))));
      }
      
      if(!reinitCount)
      {
        std::stringstream ss;
        ss << "Connectome::setEvents: failed to generate a connection "
          << "without direct feedback for neuron ID " << i << std::endl;
        throw SimException(ss.str());
      }
      
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

  LOG("Gaba Synapse Count:" << gabaSynapseCount, 1);
  LOG("Ampa Synapse Count:" << ampaSynapseCount, 1);
  LOG("Total Synapse Count:" << (ampaSynapseCount+gabaSynapseCount), 1);
  LOG("Max Synapses/Neuron:" << maxSynapsePerNeuron, 1);
  LOG("Min Synapses/Neuron:" << minSynapsePerNeuron, 1);
  LOG("Average Synapses/Neuron:" << averageSynapsePerNeuron, 1);

  this->storeBuffers(queue, block);
}
/**************************************************************************************************/



cl_uint
Connectome::getSynapseCount
(
  cl::CommandQueue  &queue,
  cl_uint           neuronID
)
/**************************************************************************************************/
{
#if CONNECTOME_VALIDATION_ENABLE
  this->isInitialized();
  
  if(neuronID >= this->neuronCount)
  {
    throw SimException("Connectome::getSynapseCount: neuron ID exceeds neuron count");
  }
#endif
  // Disable temporarely since device doesn't change connectome w/o STDP
  //this->getConnections(queue, CL_TRUE);

  cl_uint ptrStart = this->dataSynapsePointer[neuronID];
  cl_uint ptrEnd = this->dataSynapsePointer[neuronID+1];
  
#if CONNECTOME_VALIDATION_ENABLE
  if(ptrStart > ptrEnd)
  {
    std::stringstream ss;
    ss << "Connectome::getSynapseCount: negative connection count (" << ptrEnd - ptrStart 
      << ") detected for neuron ID " << neuronID << std::endl;
    throw SimException(ss.str());
  }
#endif

  return ptrEnd - ptrStart;
}
/**************************************************************************************************/



void
Connectome::getSynapse
(
  cl::CommandQueue  &queue,
  cl_uint           neuronID,
  cl_uint           synapse,
  cl_uint           &targetNeuron,
  CL_DATA_TYPE      &delay,
  CL_DATA_TYPE      &weight
)
/**************************************************************************************************/
{
#if CONNECTOME_VALIDATION_ENABLE
  this->isInitialized();
  
  if(neuronID >= this->neuronCount)
  {
    throw SimException("Connectome::getSynapse: neuron ID exceeds neuron count");
  }
#endif
  // Disable temporarely since device doesn't change connectome w/o STDP
  //this->getConnections(queue, CL_TRUE);
  
  cl_uint ptrStart = this->dataSynapsePointer[neuronID];
  cl_uint ptrEnd = this->dataSynapsePointer[neuronID+1];
  
#if CONNECTOME_VALIDATION_ENABLE
  if(ptrStart > ptrEnd)
  {
    std::stringstream ss;
    ss << "Connectome::getSynapse: negative connection count (" << ptrEnd - ptrStart 
      << ") detected for neuron ID " << neuronID << std::endl;
    throw SimException(ss.str());
  }
  
  cl_uint synapseCount = ptrEnd - ptrStart;
  
  if(synapse >= synapseCount)
  {
    std::stringstream ss;
    ss << "Connectome::getSynapse: synapse ID " << synapse << " exceeds synapse count " 
      << synapseCount << "for neuron ID " << neuronID << std::endl;
    throw SimException(ss.str());
  }
#endif
  
  cl_uint ptr = ptrStart + synapse;
  targetNeuron = this->dataSynapseTargets[ptr];
  delay = this->dataSynapseDelays[ptr];
  weight = this->dataSynapseWeights[ptr];
}
/**************************************************************************************************/



void 
Connectome::refresh
(
  cl::CommandQueue    &queue
)
/**************************************************************************************************/
{
#if CONNECTOME_VALIDATION_ENABLE
  this->isInitialized();
#endif

  // Disable temporarely since device doesn't change connectome w/o STDP
  //this->getConnections(queue, CL_TRUE);
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Connectome::getConnections
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  if(!this->dataValid)
  {
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSynapsePointerBuffer, 
      this->dataSynapsePointerSizeBytes, this->dataSynapsePointer);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSynapseTargetsBuffer, 
      this->dataSynapseTargetsSizeBytes, this->dataSynapseTargets);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSynapseDelaysBuffer, 
      this->dataSynapseDelaysSizeBytes, this->dataSynapseDelays);
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataSynapseWeightsBuffer, 
      this->dataSynapseWeightsSizeBytes, this->dataSynapseWeights);

    this->dataValid = true;
  }
}
/**************************************************************************************************/



void
Connectome::reset(bool checkForNull)
/**************************************************************************************************/
{
  try
  {
  this->resetObject = true;
  this->dataValid = false;
  this->srandSeed = 0;
  this->srandCounter = 0;
  this->neuronCount = 0;
  this->dataToSimulationLogFile = NULL;
  this->dataToReportLogFile = NULL;
    
  this->dataSynapsePointerSize = 0;
  this->dataSynapsePointerSizeBytes = 0;

  this->dataSynapseTargetsSize = 0;
  this->dataSynapseTargetsSizeBytes = 0;
  
  this->dataSynapseDelaysSize = 0;
  this->dataSynapseDelaysSizeBytes = 0;
  
  this->dataSynapseWeightsSize = 0;
  this->dataSynapseWeightsSizeBytes = 0;
  
  if(checkForNull)
  {
    if(this->dataSynapsePointer)
    {
      free(this->dataSynapsePointer);
      this->dataSynapsePointer = NULL;
    }
    if(this->dataSynapseTargets)
    {
      free(this->dataSynapseTargets);
      this->dataSynapseTargets = NULL;
    }
    if(this->dataSynapseDelays)
    {
      free(this->dataSynapseDelays);
      this->dataSynapseDelays = NULL;
    }
    if(this->dataSynapseWeights)
    {
      free(this->dataSynapseWeights);
      this->dataSynapseWeights = NULL;
    }
  }
  else
  {
    this->dataSynapsePointer = NULL;
    this->dataSynapseTargets = NULL;
    this->dataSynapseDelays = NULL;
    this->dataSynapseWeights = NULL;
  }
  }
  CATCH(std::cerr, Connectome::reset, throw SimException("Connectome::reset: failed.");)
}
/**************************************************************************************************/



void
Connectome::storeBuffers
(
  cl::CommandQueue  &queue,
  cl_bool           block
)
/**************************************************************************************************/
{
  ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSynapsePointerBuffer, 
    this->dataSynapsePointerSizeBytes, this->dataSynapsePointer);
  ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSynapseTargetsBuffer, 
    this->dataSynapseTargetsSizeBytes, this->dataSynapseTargets);
  ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSynapseDelaysBuffer, 
    this->dataSynapseDelaysSizeBytes, this->dataSynapseDelays);
  ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataSynapseWeightsBuffer, 
    this->dataSynapseWeightsSizeBytes, this->dataSynapseWeights);
    
  this->dataValid = true;
}
/**************************************************************************************************/



void
Connectome::isInitialized
()
/**************************************************************************************************/
{
  if(this->resetObject)
  {
    throw SimException("Connectome::isInitialized: the object was not initialized");
  }
}
/**************************************************************************************************/
