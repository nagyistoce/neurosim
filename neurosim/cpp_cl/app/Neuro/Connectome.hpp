
#ifndef CONNECTOME_H_
#define CONNECTOME_H_



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Common.hpp"

/**************************************************************************************************/



/**
  @class Connectome

  Models connectome (synaptic connections).

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/04/23
 */
class Connectome 
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



  bool resetObject;
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
  Connectome
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Connectome
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes this object.
  */
  void
  initialize
  (
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    cl_uint, 
    cl_uint,
    double, 
    struct kernelStatistics*,
    std::stringstream*,
    std::stringstream*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Sets all connections according to parameters
  */
  void 
  setConnections
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
    cl::CommandQueue&,
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
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint&,
    CL_DATA_TYPE&,
    CL_DATA_TYPE&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current spike events (they become past spikes).
  */
  void 
  refresh
  (
    cl::CommandQueue&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  private:  /*private methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads connectome from the device if it is invalid on the host.
  */
  void
  getConnections
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Resets all variables to the default state.
  */
  void
  reset
  (
    bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Stores connectome to the device.
  */
  void
  storeBuffers
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies if this object was initialized.
  */
  void
  isInitialized
  ();
/**************************************************************************************************/
};

#endif
