
#ifndef OPERATOR_GROUP_H_
#define OPERATOR_GROUP_H_



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



#if ENABLE_OPERATOR_GROUP



/***************************************************************************************************
  Class Preprocessor Definitions
***************************************************************************************************/

#define OPERATOR_GROUP_VALID_DEBUG         0x1
#define OPERATOR_GROUP_VALID_ERROR         0x2
#define OPERATOR_GROUP_VALID_ALL           (OPERATOR_GROUP_VALID_DEBUG |\
                                            OPERATOR_GROUP_VALID_ERROR)

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/

class Data_SynapticEvents;
class Operator_Sort;

/**************************************************************************************************/



/**
  @class Operator_Group

  Models spike events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/09/15
 */
class Operator_Group 
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

#if (GROUP_EVENTS_DEBUG_ENABLE)
  size_t dataGroupEventsDebugHostSize;
  size_t dataGroupEventsDebugHostSizeBytes;
  cl_uint* dataGroupEventsDebugHost;
  cl::Buffer dataGroupEventsDebugHostBuffer;
  
  size_t dataGroupEventsDebugDeviceSize;
  size_t dataGroupEventsDebugDeviceSizeBytes;
  cl_uint* dataGroupEventsDebugDevice;
  cl::Buffer dataGroupEventsDebugDeviceBuffer;
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  size_t dataGroupEventsErrorSize;
  size_t dataGroupEventsErrorSizeBytes;
  cl_uint* dataGroupEventsError;
  cl::Buffer dataGroupEventsErrorBuffer;
#endif

#if GROUP_EVENTS_ENABLE_V00
  cl_uint histogramBinShift_v00;
  cl_uint histogramOutBinShift_v00;
  cl_uint valuesMode_v00;
  cl_uint sourceEventsDataStructureType_v00;
  bool relocateValues_v00;
  bool enableStepShift_v00;
  bool setKernelArguments_v00;
  size_t blockSizeX_kernelGroupEventsV00;
  size_t blockSizeY_kernelGroupEventsV00;
  cl_uint argNumGroupEventsV00;
  cl::Kernel kernelGroupEventsV00;
  cl::NDRange *globalThreadsGroupEventsV00;
  cl::NDRange *localThreadsGroupEventsV00;
#endif

#if GROUP_EVENTS_ENABLE_V01
  cl_uint histogramBinShift_v01;
  cl_uint histogramOutBinShift_v01;
  cl_uint valuesMode_v01;
  cl_uint sourceEventsDataStructureType_v01;
  bool relocateValues_v01;
  bool enableStepShift_v01;
  bool setKernelArguments_v01;
  size_t blockSizeX_kernelGroupEventsV01;
  size_t blockSizeY_kernelGroupEventsV01;
  cl_uint argNumGroupEventsV01;
  cl::Kernel kernelGroupEventsV01;
  cl::NDRange *globalThreadsGroupEventsV01;
  cl::NDRange *localThreadsGroupEventsV01;
#endif

#if GROUP_EVENTS_ENABLE_V02
  cl_uint histogramBinShift_v02;
  cl_uint histogramOutBinShift_v02;
  cl_uint valuesMode_v02;
  cl_uint sourceEventsDataStructureType_v02;
  bool relocateValues_v02;
  bool enableStepShift_v02;
  bool setKernelArguments_v02;
  size_t blockSizeX_kernelGroupEventsV02;
  size_t blockSizeY_kernelGroupEventsV02;
  cl_uint argNumGroupEventsV02;
  cl::Kernel kernelGroupEventsV02;
  cl::NDRange *globalThreadsGroupEventsV02;
  cl::NDRange *localThreadsGroupEventsV02;
#endif

#if GROUP_EVENTS_ENABLE_V03
  cl_uint histogramBinShift_v03;
  cl_uint histogramOutBinShift_v03;
  cl_uint valuesMode_v03;
  cl_uint sourceEventsDataStructureType_v03;
  bool relocateValues_v03;
  bool enableStepShift_v03;
  bool setKernelArguments_v03;
  size_t blockSizeX_kernelGroupEventsV03;
  size_t blockSizeY_kernelGroupEventsV03;
  cl_uint argNumGroupEventsV03;
  cl::Kernel kernelGroupEventsV03;
  cl::NDRange *globalThreadsGroupEventsV03;
  cl::NDRange *localThreadsGroupEventsV03;
#endif

  cl_uint timeSlots;
  cl_uint neuronCount;
  cl_uint histogramBinCount;
  cl_uint histogramBinSize;
  cl_uint histogramBinMask;
  cl_uint histogramOutBinCount;
  cl_uint histogramOutBinMask;
  cl_uint histogramOutGridSize;
  cl_uint eventDataSrcBufCount;
  cl_uint eventDataSrcBufMaxSize;
  cl_uint eventDataDstBufMaxSize;
  cl_uint eventDataUnitSize;
  cl_uint gridSizeWg;
  cl_uint wgSize;
  cl_uint wiDataCount;
  cl_uint wgEventBufCount;
  double  minimumDelay;
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
  Operator_Group
  (
#if (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
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
#if (GROUP_EVENTS_DEBUG_ENABLE)
    dataGroupEventsDebugHostSize(0),
    dataGroupEventsDebugHostSizeBytes(0),
    dataGroupEventsDebugHost(NULL),
    dataGroupEventsDebugHostBuffer(),
    dataGroupEventsDebugDeviceSize(0),
    dataGroupEventsDebugDeviceSizeBytes(0),
    dataGroupEventsDebugDevice(NULL),
    dataGroupEventsDebugDeviceBuffer(),
#endif
    /* *** */
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    dataGroupEventsErrorSize(0),
    dataGroupEventsErrorSizeBytes(0),
    dataGroupEventsError(NULL),
    dataGroupEventsErrorBuffer(),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V00
    histogramBinShift_v00(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00),
    histogramOutBinShift_v00(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V00),
    valuesMode_v00(GROUP_EVENTS_VALUES_MODE_V00),
    sourceEventsDataStructureType_v00(GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V00),
    relocateValues_v00(GROUP_EVENTS_RELOCATE_VALUES_V00),
    enableStepShift_v00(GROUP_EVENTS_ENABLE_STEP_SHIFT_V00),
    setKernelArguments_v00(true),
    blockSizeX_kernelGroupEventsV00(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV00(1),
    argNumGroupEventsV00(0),
    kernelGroupEventsV00(),
    globalThreadsGroupEventsV00(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG)),
    localThreadsGroupEventsV00(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI, 1)),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V01
    histogramBinShift_v01(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01),
    histogramOutBinShift_v01(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01),
    valuesMode_v01(GROUP_EVENTS_VALUES_MODE_V01),
    sourceEventsDataStructureType_v01(GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V01),
    relocateValues_v01(GROUP_EVENTS_RELOCATE_VALUES_V01),
    enableStepShift_v01(GROUP_EVENTS_ENABLE_STEP_SHIFT_V01),
    setKernelArguments_v01(true),
    blockSizeX_kernelGroupEventsV01(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV01(1),
    argNumGroupEventsV01(0),
    kernelGroupEventsV01(),
    globalThreadsGroupEventsV01(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG)),
    localThreadsGroupEventsV01(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI, 1)),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V02
    histogramBinShift_v02(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02),
    histogramOutBinShift_v02(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V02),
    valuesMode_v02(GROUP_EVENTS_VALUES_MODE_V02),
    sourceEventsDataStructureType_v02(GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V02),
    relocateValues_v02(GROUP_EVENTS_RELOCATE_VALUES_V02),
    enableStepShift_v02(GROUP_EVENTS_ENABLE_STEP_SHIFT_V02),
    setKernelArguments_v02(true),
    blockSizeX_kernelGroupEventsV02(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV02(1),
    argNumGroupEventsV02(0),
    kernelGroupEventsV02(),
    globalThreadsGroupEventsV02(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG)),
    localThreadsGroupEventsV02(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI, 1)),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V03
    histogramBinShift_v03(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03),
    histogramOutBinShift_v03(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V03),
    valuesMode_v03(GROUP_EVENTS_VALUES_MODE_V03),
    sourceEventsDataStructureType_v03(GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V03),
    relocateValues_v03(GROUP_EVENTS_RELOCATE_VALUES_V03),
    enableStepShift_v03(GROUP_EVENTS_ENABLE_STEP_SHIFT_V03),
    setKernelArguments_v03(true),
    blockSizeX_kernelGroupEventsV03(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV03(1),
    argNumGroupEventsV03(0),
    kernelGroupEventsV03(),
    globalThreadsGroupEventsV03(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG)),
    localThreadsGroupEventsV03(new cl::NDRange(GROUP_EVENTS_WG_SIZE_WI, 1)),
#endif
    /* *** */
    timeSlots(GROUP_EVENTS_TIME_SLOTS),
    neuronCount(GROUP_EVENTS_TOTAL_NEURONS),
    histogramBinCount(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS),
    histogramBinSize(GROUP_EVENTS_HISTOGRAM_BIN_SIZE),
    histogramBinMask(GROUP_EVENTS_HISTOGRAM_BIN_MASK),
    histogramOutBinCount(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT),
    histogramOutBinMask(GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT),
    histogramOutGridSize(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE),
    eventDataSrcBufCount(GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS),
    eventDataSrcBufMaxSize(GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE),
    eventDataDstBufMaxSize(GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE),
    eventDataUnitSize(GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS),
    gridSizeWg(GROUP_EVENTS_GRID_SIZE_WG),
    wgSize(GROUP_EVENTS_WG_SIZE_WI),
    wiDataCount(GROUP_EVENTS_ELEMENTS_PER_WI),
    wgEventBufCount(GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG),
    minimumDelay(GROUP_EVENTS_MIN_DELAY),
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile),
    dataValid(0)
/**************************************************************************************************/
  {
    /* Must be called only from the constructor and only once*/
    this->initialize
    (
#if (GROUP_EVENTS_DEBUG_ENABLE)
      GROUP_EVENTS_DEBUG_BUFFER_SIZE_WORDS,
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
      GROUP_EVENTS_ERROR_BUFFER_SIZE_WORDS,
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
      queue,
      block,
#endif
#if GROUP_EVENTS_ENABLE_V00
      GROUP_EVENTS_CACHE_SIZE_WORDS_V00,
#endif
#if GROUP_EVENTS_ENABLE_V01
      GROUP_EVENTS_CACHE_SIZE_WORDS_V01,
#endif
#if GROUP_EVENTS_ENABLE_V02
      GROUP_EVENTS_CACHE_SIZE_WORDS_V02,
#endif
#if GROUP_EVENTS_ENABLE_V03
      GROUP_EVENTS_CACHE_SIZE_WORDS_V03,
#endif
      context,
      device,
      kernelStats
    );
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Operator_Group()
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
    if(this->dataGroupEventsDebugHost)
    {
      free(this->dataGroupEventsDebugHost);
      this->dataGroupEventsDebugHost = NULL;
    }
    if(this->dataGroupEventsDebugDevice)
    {
      free(this->dataGroupEventsDebugDevice);
      this->dataGroupEventsDebugDevice = NULL;
    }
#endif
    /* *** */
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(this->dataGroupEventsError)
    {
      free(this->dataGroupEventsError);
      this->dataGroupEventsError = NULL;
    }
#endif
    /* *** */
#if (GROUP_EVENTS_ENABLE_V00)
    if(this->globalThreadsGroupEventsV00)
    {
      delete(this->globalThreadsGroupEventsV00);
    }
    if(this->localThreadsGroupEventsV00)
    {
      delete(this->localThreadsGroupEventsV00);
    }
#endif 
    /* *** */
#if (GROUP_EVENTS_ENABLE_V01)
    if(this->globalThreadsGroupEventsV01)
    {
      delete(this->globalThreadsGroupEventsV01);
    }
    if(this->localThreadsGroupEventsV01)
    {
      delete(this->localThreadsGroupEventsV01);
    }
#endif 
    /* *** */
#if (GROUP_EVENTS_ENABLE_V02)
    if(this->globalThreadsGroupEventsV02)
    {
      delete(this->globalThreadsGroupEventsV02);
    }
    if(this->localThreadsGroupEventsV02)
    {
      delete(this->localThreadsGroupEventsV02);
    }
#endif 
    /* *** */
#if (GROUP_EVENTS_ENABLE_V03)
    if(this->globalThreadsGroupEventsV03)
    {
      delete(this->globalThreadsGroupEventsV03);
    }
    if(this->localThreadsGroupEventsV03)
    {
      delete(this->localThreadsGroupEventsV03);
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
#if GROUP_EVENTS_DEBUG_ENABLE
  void 
  invalidateDebug
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates error data on the host.
  */
#if GROUP_EVENTS_ERROR_TRACK_ENABLE
  void 
  invalidateError
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs group operation.
  */
#if GROUP_EVENTS_ENABLE_V00
  void 
  group_v00
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs unit test of group operation.
  */
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V00
  void 
  groupUnitTest_v00
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    int,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies group operation.
  */
#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V00
  void 
  verifyGroup_v00
  (
    cl::CommandQueue&,
    bool,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs group operation.
  */
#if GROUP_EVENTS_ENABLE_V01
  void 
  group_v01
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs unit test of group operation.
  */
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V01
  void 
  groupUnitTest_v01
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    double,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies group operation.
  */
#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V01
  void 
  verifyGroup_v01
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs group operation.
  */
#if GROUP_EVENTS_ENABLE_V02
  void 
  group_v02
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs unit test of group operation.
  */
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V02
  void 
  groupUnitTest_v02
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    double,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies group operation.
  */
#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V02
  void 
  verifyGroup_v02
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs group operation.
  */
#if GROUP_EVENTS_ENABLE_V03
  void 
  group_v03
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs unit test of group operation.
  */
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V03
  void 
  groupUnitTest_v03
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    double,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies group operation.
  */
#if GROUP_EVENTS_VERIFY_ENABLE && GROUP_EVENTS_ENABLE_V03
  void 
  verifyGroup_v03
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    Data_SynapticEvents&
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
#if (GROUP_EVENTS_DEBUG_ENABLE)
    size_t,
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    size_t,
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
    cl::CommandQueue&,
    cl_bool,
#endif
#if GROUP_EVENTS_ENABLE_V00
    size_t,
#endif
#if GROUP_EVENTS_ENABLE_V01
    size_t,
#endif
#if GROUP_EVENTS_ENABLE_V02
    size_t,
#endif
#if GROUP_EVENTS_ENABLE_V03
    size_t,
#endif
    cl::Context&,
    cl::Device&,
    struct kernelStatistics&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes grouped events for verification
  */
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V01
  bool 
  setGrouppedEvents
  (
    std::stringstream*,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    size_t,
    cl_uint,
    cl_uint*,
    cl_uint*,
    cl_uint*,
    cl_uint*,
    Data_SynapticEvents&,
    cl_uint*,
    cl_uint*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Sets reference data and verifies group operation.
  */
#if GROUP_EVENTS_VERIFY_ENABLE
  void 
  setReferenceAndVerifyGroup
  (
    bool,
    bool,
    bool,
    bool,
    bool,
    bool,
    bool,
    int,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    Data_SynapticEvents&,
    cl::CommandQueue&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies group operation generic to all versions
  */
#if GROUP_EVENTS_VERIFY_ENABLE
  bool 
  verifyGroup
  (
    bool,
    bool,
    bool,
    bool,
    bool,
    bool,
    bool,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint*,
    cl_uint*,
    std::stringstream*,
    std::stringstream*,
    cl::CommandQueue&,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads data from the device if it is invalid on the host.
  */
#if GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
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
#if GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
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



#endif /*ENABLE_OPERATOR_GROUP*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif /*OPERATOR_GROUP_H_*/
