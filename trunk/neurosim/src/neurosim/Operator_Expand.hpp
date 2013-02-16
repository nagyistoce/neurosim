
#ifndef OPERATOR_EXPAND_H_
#define OPERATOR_EXPAND_H_



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



#if ENABLE_OPERATOR_EXPAND



/***************************************************************************************************
  Class Preprocessor Definitions
***************************************************************************************************/

#define OPERATOR_EXPAND_VALID_DEBUG         0x1
#define OPERATOR_EXPAND_VALID_ERROR         0x2
#define OPERATOR_EXPAND_VALID_ALL           (OPERATOR_EXPAND_VALID_DEBUG |\
                                            OPERATOR_EXPAND_VALID_ERROR)

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/

class Data_Connectome;
class Data_SynapticEvents;
class Data_SpikeEvents;

/**************************************************************************************************/



/**
  @class Operator_Expand

  Models spike events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/09/15
 */
class Operator_Expand 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
  __int64 performanceCounter, performanceFrequency;
#endif

  bool setKernelArguments;
  size_t blockSizeX_KernelExpandEvents;
  size_t blockSizeY_KernelExpandEvents;
  cl_uint argNumExpandEvents;
  cl::Kernel kernelExpandEvents;
  cl::NDRange *globalThreadsExpandEvents;
  cl::NDRange *localThreadsExpandEvents;
  
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  size_t dataExpandEventsDebugHostSize;
  size_t dataExpandEventsDebugHostSizeBytes;
  cl_uint* dataExpandEventsDebugHost;
  cl::Buffer dataExpandEventsDebugHostBuffer;
  
  size_t dataExpandEventsDebugDeviceSize;
  size_t dataExpandEventsDebugDeviceSizeBytes;
  cl_uint* dataExpandEventsDebugDevice;
  cl::Buffer dataExpandEventsDebugDeviceBuffer;
#endif

#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  size_t dataExpandEventsErrorSize;
  size_t dataExpandEventsErrorSizeBytes;
  cl_uint* dataExpandEventsError;
  cl::Buffer dataExpandEventsErrorBuffer;
#endif

  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;
  
  cl_uint dataValid;
  
  
  
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  Operator_Expand
  (
#if (EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    cl_bool                             block,
    cl::CommandQueue                    &queue,
#endif
    cl::Context                         &context,
    cl::Device                          &device,
    struct kernelStatistics             &kernelStats,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile
  ) :
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
    performanceCounter(0),
    performanceFrequency(0),
#endif
    /* *** */
    setKernelArguments(true),
    blockSizeX_KernelExpandEvents(EXPAND_EVENTS_WG_SIZE_WI),
    blockSizeY_KernelExpandEvents(1),
    argNumExpandEvents(0),
    kernelExpandEvents(),
    globalThreadsExpandEvents(new cl::NDRange(EXPAND_EVENTS_WG_SIZE_WI*EXPAND_EVENTS_GRID_SIZE_WG)),
    localThreadsExpandEvents(new cl::NDRange(EXPAND_EVENTS_WG_SIZE_WI, 1)),
    /* *** */
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    dataExpandEventsDebugHostSize(0),
    dataExpandEventsDebugHostSizeBytes(0),
    dataExpandEventsDebugHost(NULL),
    dataExpandEventsDebugHostBuffer(),
    dataExpandEventsDebugDeviceSize(0),
    dataExpandEventsDebugDeviceSizeBytes(0),
    dataExpandEventsDebugDevice(NULL),
    dataExpandEventsDebugDeviceBuffer(),
#endif
    /* *** */
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    dataExpandEventsErrorSize(0),
    dataExpandEventsErrorSizeBytes(0),
    dataExpandEventsError(NULL),
    dataExpandEventsErrorBuffer(),
#endif
    /* *** */
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile),
    dataValid(0)
/**************************************************************************************************/
  {
    /* Must be called only from the constructor and only once*/
    this->initialize
    (
#if (EXPAND_EVENTS_DEBUG_ENABLE)
      EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS,
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
      EXPAND_EVENTS_ERROR_BUFFER_SIZE_WORDS,
#endif
#if (EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE)
      queue,
      block,
#endif
      context,
      device,
      kernelStats,
      EXPAND_EVENTS_CACHE_SIZE_WORDS
    );
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Operator_Expand()
  {
    if(this->globalThreadsExpandEvents)
    {
      delete(this->globalThreadsExpandEvents);
    }
    if(this->localThreadsExpandEvents)
    {
      delete(this->localThreadsExpandEvents);
    }
    /* *** */
#if (EXPAND_EVENTS_DEBUG_ENABLE)
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
    /* *** */
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    if(this->dataExpandEventsError)
    {
      free(this->dataExpandEventsError);
      this->dataExpandEventsError = NULL;
    }
#endif
    /* *** */
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates debug data on the host.
  */
#if EXPAND_EVENTS_DEBUG_ENABLE
  void 
  invalidateDebug
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates error data on the host.
  */
#if EXPAND_EVENTS_ERROR_TRACK_ENABLE
  void 
  invalidateError
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs expand operation.
  */
  void 
  expand
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    cl_uint,
    Data_SpikeEvents&,
    Data_SynapticEvents&,
    Data_Connectome&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs expand operation.
  */
  void 
  expand
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    double,
    Data_SpikeEvents&,
    Data_SynapticEvents&,
    Data_Connectome&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs expand operation with a unit test.
  */
#if ENABLE_UNIT_TEST_EXPAND_EVENTS
  void 
  expandUnitTest
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
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
    Data_SpikeEvents&,
    Data_SynapticEvents&,
    Data_Connectome&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies expand operation.
  */
#if EXPAND_EVENTS_VERIFY_ENABLE
  void 
  verifyExpand
  (
    cl::CommandQueue&,
    bool,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    Data_SpikeEvents&,
    Data_SynapticEvents&,
    Data_Connectome&
  );
#endif
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
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    size_t,
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    size_t,
#endif
#if (EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    cl::CommandQueue&,
    cl_bool,
#endif
    cl::Context&,
    cl::Device&,
    struct kernelStatistics&,
    size_t
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads data from the device if it is invalid on the host.
  */
#if EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE
  void
  getData
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Stores data on device.
  */
#if EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE
  void
  storeData
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
#endif
/**************************************************************************************************/
};



#endif /*ENABLE_OPERATOR_EXPAND*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif /*OPERATOR_EXPAND_H_*/
