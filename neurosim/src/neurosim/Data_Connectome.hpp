
#ifndef CONNECTOME_H_
#define CONNECTOME_H_



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Common.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



/***************************************************************************************************
  Class Preprocessor Definitions
***************************************************************************************************/



/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/



/**************************************************************************************************/



/**
  @class Data_Connectome

  Models connectome (synaptic connections).

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/04/23
 */
class Data_Connectome 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/
  
  
  
  cl::Buffer dataSynapsePointerBuffer;
  cl::Buffer dataSynapseTargetsBuffer;
  cl::Buffer dataSynapseDelaysBuffer;
  cl::Buffer dataSynapseWeightsBuffer;
  
  
  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



  bool dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;
  
  cl_uint neuronCount;
  
  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;
  
  size_t dataSynapsePointerSize;
  size_t dataSynapsePointerSizeBytes;
  cl_uint* dataSynapsePointer;
  
  size_t dataSynapseTargetsSize;
  size_t dataSynapseTargetsSizeBytes;
  cl_uint* dataSynapseTargets;
  
  size_t dataSynapseDelaysSize;
  size_t dataSynapseDelaysSizeBytes;
  cl_float* dataSynapseDelays;
  
  size_t dataSynapseWeightsSize;
  size_t dataSynapseWeightsSizeBytes;
  cl_float* dataSynapseWeights;
  
  
  
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  Data_Connectome
  (
    cl::Context                         &context,
    cl::Device                          &device,
    cl::CommandQueue                    &queue,
    cl_bool                             block,
    struct kernelStatistics             &kernelStats,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile,
    cl_uint                             neuronCount,
    cl_uint                             maxConnectionsPerNeuron,
    double                              connectionDeviationRatio,
    double                              gabaPercent,
    double                              minDelay,
    double                              maxDelay
  ) : 
    dataSynapsePointerBuffer(),
    dataSynapseTargetsBuffer(),
    dataSynapseDelaysBuffer(),
    dataSynapseWeightsBuffer(),
    dataValid(0),
    srandSeed(1),
    srandCounter(0),
    neuronCount(neuronCount),
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile),
    dataSynapsePointerSize(0),
    dataSynapsePointerSizeBytes(0),
    dataSynapsePointer(NULL),
    dataSynapseTargetsSize(0),
    dataSynapseTargetsSizeBytes(0),
    dataSynapseTargets(NULL),
    dataSynapseDelaysSize(0),
    dataSynapseDelaysSizeBytes(0),
    dataSynapseDelays(NULL),
    dataSynapseWeightsSize(0),
    dataSynapseWeightsSizeBytes(0),
    dataSynapseWeights(NULL)
/**************************************************************************************************/
  {
    /* Must be called only from the constructor and only once*/
    this->initialize
    (
      context,
      device,
      queue,
      block,
      kernelStats,
      maxConnectionsPerNeuron,
      connectionDeviationRatio
    );
    
    this->resetConnections
    (
      queue,
      block,
      gabaPercent,
      minDelay,
      maxDelay
    );
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Data_Connectome()
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
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Sets all connections according to parameters
  */
  void 
  resetConnections
  (
    cl::CommandQueue&,
    cl_bool,
    double,
    double,
    double
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns synapse count.
  */
  cl_uint
  getSynapseCount
  (
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns synapse datum (target neuron ID, weight, delay) by reference for a given source neurons
    ID and synapse ID.
  */
  void
  getSynapse
  (
    cl_uint,
    cl_uint,
    cl_uint&,
    CL_DATA_TYPE&,
    CL_DATA_TYPE&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  private:  /*private methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes this object. Must be called only from the constructor and only once.
  */
  void
  initialize
  (
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    struct kernelStatistics&,
    cl_uint,
    double
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads data from the device if it is invalid on the host.
  */
  void
  getData
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Stores data on device.
  */
  void
  storeData
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/
};



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif /*CONNECTOME_H_*/
