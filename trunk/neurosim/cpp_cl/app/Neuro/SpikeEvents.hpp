
#ifndef SPIKE_EVENTS_H_
#define SPIKE_EVENTS_H_



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Common.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Class Preprocessor Definitions
***************************************************************************************************/

#define SPIKE_EVENTS_VALID_SPIKES         0x1
#define SPIKE_EVENTS_VALID_EXPAND_DEBUG   0x2
#define SPIKE_EVENTS_VALID_EXPAND_ERROR   0x4

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/

class Connectome;
class SynapticEvents;

/**************************************************************************************************/



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
  cl_uint dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;
  
  cl_uint neuronCount;
  cl_uint spikePacketSize;
  cl_uint spikePacketSizeWords;
  cl_uint spikePackets;
  cl_uint simulationStepSize;
  cl_uint spikeDatumSize;
  
  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;
  
  size_t dataSpikePacketsSize;
  size_t dataSpikePacketsSizeBytes;
  cl_uint *dataSpikePackets;
  
  size_t dataSpikePacketCountsSize;
  size_t dataSpikePacketCountsSizeBytes;
  cl_uint *dataSpikePacketCounts;

  cl_uint *dataPastSpikePackets;
  cl_uint *dataPastSpikePacketCounts;
  
#if EXPAND_EVENTS_ENABLE
  bool setKernelArguments;
  cl_uint lmExpandEventsSizeBytes;
  size_t blockSizeX_KernelExpandEvents;
  size_t blockSizeY_KernelExpandEvents;
  cl_uint argNumExpandEvents;
  cl::Kernel kernelExpandEvents;
  cl::NDRange *globalThreadsExpandEvents;
  cl::NDRange *localThreadsExpandEvents;
  
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
  size_t dataExpandEventsDebugHostSize;
  size_t dataExpandEventsDebugHostSizeBytes;
  cl_uint* dataExpandEventsDebugHost;
  cl::Buffer dataExpandEventsDebugHostBuffer;
  
  size_t dataExpandEventsDebugDeviceSize;
  size_t dataExpandEventsDebugDeviceSizeBytes;
  cl_uint* dataExpandEventsDebugDevice;
  cl::Buffer dataExpandEventsDebugDeviceBuffer;
#endif

#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  size_t dataExpandEventsErrorSize;
  size_t dataExpandEventsErrorSizeBytes;
  cl_uint* dataExpandEventsError;
  cl::Buffer dataExpandEventsErrorBuffer;
#endif
#endif



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
    cl::CommandQueue&,
    cl_bool,
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
    Invalidates current spike events (they become past spikes).
  */
  void 
  invalidateEvents
  ();
/**************************************************************************************************/



#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_DEBUG_ENABLE)
/**************************************************************************************************/
  /**
    Invalidates current spike events (they become past spikes).
  */
  void 
  invalidateDebug
  ();
/**************************************************************************************************/
#endif



#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ERROR_TRACK_ENABLE)
/**************************************************************************************************/
  /**
    Invalidates current spike events (they become past spikes).
  */
  void 
  invalidateError
  ();
/**************************************************************************************************/
#endif



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



#if EXPAND_EVENTS_ENABLE
/**************************************************************************************************/
  /**
    Loads current spike events from the device if they are invalid on the host.
  */
  void 
  expand
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl_uint,
    SynapticEvents&,
    Connectome&,
    cl::CommandQueue&,
    cl::Event&
  );
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ENABLE
/**************************************************************************************************/
  /**
    Loads current spike events from the device if they are invalid on the host.
  */
  void 
  expand
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl_uint,
    cl_uint,
    cl_int,
    double,
    double,
    double,
    double,
    double,
    double,
    SynapticEvents&,
    Connectome&,
    cl::CommandQueue&,
    cl::Event&
  );
/**************************************************************************************************/
#endif



#if EXPAND_EVENTS_ENABLE
/**************************************************************************************************/
  /**
    Loads current spike events from the device if they are invalid on the host.
  */
  void 
  expand
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    double,
    SynapticEvents&,
    Connectome&,
    cl::CommandQueue&,
    cl::Event&
  );
/**************************************************************************************************/
#endif



/**************************************************************************************************/
  private:  /*private methods*/
/**************************************************************************************************/



#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_VERIFY_ENABLE)
/**************************************************************************************************/
  /**
    Loads current spike events from the device if they are invalid on the host.
  */
  void 
  verifyExpand
  (
    cl_uint,
    bool,
    SynapticEvents&,
    Connectome&,
    cl::CommandQueue&
  );
/**************************************************************************************************/
#endif



/**************************************************************************************************/
  /**
    Loads current spike events from the device if they are invalid on the host.
  */
  void
  getData
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
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
    cl_bool,
    cl_uint
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
