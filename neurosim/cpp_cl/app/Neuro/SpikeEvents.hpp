
#ifndef SPIKE_EVENTS_H_
#define SPIKE_EVENTS_H_



#include "Definitions.hpp"



/**
  @class SpikeEvents

  Models spike events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/04/17
 */
class SpikeEvents 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/
  
  
  
  cl::Buffer dataSpikePacketsBuffer;
  cl::Buffer dataSpikePacketCountsBuffer;
  
  
  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



  bool resetObject;
  bool dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;
  
  cl_uint neuronCount;
  cl_uint spikePacketSize;
  cl_uint spikePacketSizeWords;
  cl_uint spikePackets;
  cl_uint simulationTimeSteps;
  cl_uint spikeDatumSize;
  
  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;
  
  cl_uint dataSpikePacketsSize;
  cl_uint dataSpikePacketsSizeBytes;
  cl_uint *dataSpikePackets;
  
  cl_uint dataSpikePacketCountsSize;
  cl_uint dataSpikePacketCountsSizeBytes;
  cl_uint *dataSpikePacketCounts;

  cl_uint *dataPastSpikePackets;
  cl_uint *dataPastSpikePacketCounts;
  
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  SpikeEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~SpikeEvents
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
    cl_uint, 
    cl_uint,
    cl_uint, 
    cl_uint,
    cl_uint, 
    cl_uint,
    struct kernelStatistics*,
    std::stringstream*,
    std::stringstream*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Clears all spike events.
  */
  void
  clearEvents
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Sets all spike events according to parameters
  */
  void 
  setEvents
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
    Returns current spike count.
  */
  cl_uint
  getSpikeCount
  (
    cl::CommandQueue&,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns the past spike count.
  */
  cl_uint
  getPastSpikeCount
  (
    cl::CommandQueue&,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns spike datum (spiked neuron and time) by reference for given spike and packet number.
  */
  void
  getSpike
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint&,
    CL_DATA_TYPE&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns by reference the past spike datum for given spike and packet number.
  */
  void
  getPastSpike
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint&,
    CL_DATA_TYPE&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Clears all past spike events.
  */
  void 
  deletePastEvents
  ();
/**************************************************************************************************/


/**************************************************************************************************/
  /**
    Invalidates current spike events (they become past spikes).
  */
  void 
  invalidateEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  private:  /*private methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads current spike events from the device if they are invalid on the host.
  */
  void
  getSpikeEvents
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Saves past spikes and loads current spikes from the device if they are invalid on the host.
  */
  void
  getPastSpikeEvents
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
    Stores spike data on device.
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
