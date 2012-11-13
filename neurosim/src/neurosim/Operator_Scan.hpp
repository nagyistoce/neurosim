
#ifndef OPERATOR_SCAN_H_
#define OPERATOR_SCAN_H_



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



#if ENABLE_OPERATOR_SCAN



/***************************************************************************************************
  Class Preprocessor Definitions
***************************************************************************************************/

#define OPERATOR_SCAN_VALID_DEBUG             0x1
#define OPERATOR_SCAN_VALID_ERROR             0x2
#define OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM  0x4
#define OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM  0x8
#define OPERATOR_SCAN_VALID_ALL               (OPERATOR_SCAN_VALID_DEBUG |\
                                              OPERATOR_SCAN_VALID_ERROR |\
                                              OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM |\
                                              OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM)

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/



/**************************************************************************************************/



/**
  @class Operator_Scan

  Prefix sum operation.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/09/01
 */
class Operator_Scan 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



#if (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
  __int64 performanceCounter, performanceFrequency;
#endif
  
#if (SCAN_DEBUG_ENABLE)
  size_t dataScanDebugHostSize;
  size_t dataScanDebugHostSizeBytes;
  cl_uint* dataScanDebugHost;
  cl::Buffer dataScanDebugHostBuffer;
  
  size_t dataScanDebugDeviceSize;
  size_t dataScanDebugDeviceSizeBytes;
  cl_uint* dataScanDebugDevice;
  cl::Buffer dataScanDebugDeviceBuffer;
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  size_t dataScanErrorSize;
  size_t dataScanErrorSizeBytes;
  cl_uint* dataScanError;
  cl::Buffer dataScanErrorBuffer;
#endif
  
#if SCAN_ENABLE_V00
  bool setKernelArguments_v00;
  size_t blockSizeX_kernelScanHistogramV00;
  size_t blockSizeY_kernelScanHistogramV00;
  cl_uint argNumScan00;
  cl::Kernel kernelScanHistogramV00;
  cl::NDRange *globalThreadsScanV00;
  cl::NDRange *localThreadsScanV00;
#endif

#if SCAN_ENABLE_V01
  bool setKernelArguments_v01;
  size_t blockSizeX_kernelScanHistogramV01;
  size_t blockSizeY_kernelScanHistogramV01;
  cl_uint argNumScan01;
  cl::Kernel kernelScanHistogramV01;
  cl::NDRange *globalThreadsScanV01;
  cl::NDRange *localThreadsScanV01;
#endif

#if (ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01)
  unsigned int srandSeed;
  unsigned int srandCounter;
#endif

#if (ENABLE_UNIT_TEST_SCAN_V00)
  cl_uint timeSlots;
  cl_uint histogramBinCount_v00;
  cl_uint histogramBinSize_v00;
  size_t dataHistogramV00Size;
  size_t dataHistogramV00SizeBytes;
  cl_uint* dataHistogramV00;
  cl_uint* dataPastHistogramV00;
  cl::Buffer dataHistogramV00Buffer;
#endif

#if (ENABLE_UNIT_TEST_SCAN_V01)
  cl_uint histogramBinBackets;
  cl_uint histogramBinCount_v01;
  cl_uint histogramBinSize_v01;
  size_t dataHistogramV01Size;
  size_t dataHistogramV01SizeBytes;
  cl_uint* dataHistogramV01;
  cl_uint* dataPastHistogramV01;
  cl::Buffer dataHistogramV01Buffer;
#endif
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
  Operator_Scan
  (
#if SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || \
    ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01
    cl_bool                             block,
    cl::CommandQueue                    &queue,
#endif
    cl::Context                         &context,
    cl::Device                          &device,
    struct kernelStatistics             &kernelStats,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile
  ) : 
#if (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)
    /* *** */
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
    performanceCounter(0),
    performanceFrequency(0),
#endif
    /* *** */
#if (SCAN_DEBUG_ENABLE)
    dataScanDebugHostSize(0),
    dataScanDebugHostSizeBytes(0),
    dataScanDebugHost(NULL),
    dataScanDebugHostBuffer(),
    dataScanDebugDeviceSize(0),
    dataScanDebugDeviceSizeBytes(0),
    dataScanDebugDevice(NULL),
    dataScanDebugDeviceBuffer(),
#endif
    /* *** */
#if (SCAN_ERROR_TRACK_ENABLE)
    dataScanErrorSize(0),
    dataScanErrorSizeBytes(0),
    dataScanError(NULL),
    dataScanErrorBuffer(),
#endif
    /* *** */
#if SCAN_ENABLE_V00
    setKernelArguments_v00(true),
    blockSizeX_kernelScanHistogramV00(SCAN_WG_SIZE_WI),
    blockSizeY_kernelScanHistogramV00(1),
    argNumScan00(0),
    kernelScanHistogramV00(),
    globalThreadsScanV00(new cl::NDRange(SCAN_WG_SIZE_WI*SCAN_GRID_SIZE_WG)),
    localThreadsScanV00(new cl::NDRange(SCAN_WG_SIZE_WI, 1)),
#endif
    /* *** */
#if SCAN_ENABLE_V01
    setKernelArguments_v01(true),
    blockSizeX_kernelScanHistogramV01(SCAN_WG_SIZE_WI),
    blockSizeY_kernelScanHistogramV01(1),
    argNumScan01(0),
    kernelScanHistogramV01(),
    globalThreadsScanV01(new cl::NDRange(SCAN_WG_SIZE_WI*SCAN_GRID_SIZE_WG)),
    localThreadsScanV01(new cl::NDRange(SCAN_WG_SIZE_WI, 1)),
#endif
    /* *** */
#if (ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01)
    srandSeed(1),
    srandCounter(0),
#endif
    /* *** */
#if (ENABLE_UNIT_TEST_SCAN_V00)
    timeSlots(SCAN_HISTOGRAM_TIME_SLOTS),
    histogramBinCount_v00(SCAN_HISTOGRAM_TOTAL_BINS_V00),
    histogramBinSize_v00(SCAN_HISTOGRAM_BIN_SIZE_V00),
    dataHistogramV00Size(0),
    dataHistogramV00SizeBytes(0),
    dataHistogramV00(NULL),
    dataPastHistogramV00(NULL),
    dataHistogramV00Buffer(),
#endif
    /* *** */
#if (ENABLE_UNIT_TEST_SCAN_V01)
    histogramBinBackets(SCAN_HISTOGRAM_BIN_BACKETS),
    histogramBinCount_v01(SCAN_HISTOGRAM_TOTAL_BINS_V01),
    histogramBinSize_v01(SCAN_HISTOGRAM_BIN_SIZE_V01),
    dataHistogramV01Size(0),
    dataHistogramV01SizeBytes(0),
    dataHistogramV01(NULL),
    dataPastHistogramV01(NULL),
    dataHistogramV01Buffer(),
#endif
    /* *** */
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
#if (SCAN_DEBUG_ENABLE)
      SCAN_DEBUG_BUFFER_SIZE_WORDS,
#endif
      /* *** */
#if (SCAN_ERROR_TRACK_ENABLE)
      SCAN_ERROR_BUFFER_SIZE_WORDS,
#endif
      /* *** */
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE)
      queue,
      block,
#endif
      context,
      device,
      kernelStats,
      SCAN_CACHE_SIZE_WORDS
    );
    
#if ENABLE_UNIT_TEST_SCAN_V00
    /* Must be called only from the constructor and only once*/
    this->setupUnitTest_v00
    (
      context,
      device,
      queue,
      block,
      kernelStats
    );
#endif

#if ENABLE_UNIT_TEST_SCAN_V01
    /* Must be called only from the constructor and only once*/
    this->setupUnitTest_v01
    (
      context,
      device,
      queue,
      block,
      kernelStats
    );
#endif
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Operator_Scan()
  {
#if (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)
    /* *** */
#if (SCAN_DEBUG_ENABLE)
    if(this->dataScanDebugHost)
    {
      free(this->dataScanDebugHost);
      this->dataScanDebugHost = NULL;
    }
    if(this->dataScanDebugDevice)
    {
      free(this->dataScanDebugDevice);
      this->dataScanDebugDevice = NULL;
    }
#endif
    /* *** */
#if (SCAN_ERROR_TRACK_ENABLE)
    if(this->dataScanError)
    {
      free(this->dataScanError);
      this->dataScanError = NULL;
    }
#endif 
    /* *** */
#if (SCAN_ENABLE_V00)
    if(this->globalThreadsScanV00)
    {
      delete(this->globalThreadsScanV00);
    }
    if(this->localThreadsScanV00)
    {
      delete(this->localThreadsScanV00);
    }
#endif 
    /* *** */
#if (SCAN_ENABLE_V01)
    if(this->globalThreadsScanV01)
    {
      delete(this->globalThreadsScanV01);
    }
    if(this->localThreadsScanV01)
    {
      delete(this->localThreadsScanV01);
    }
#endif 
    /* *** */
#if (ENABLE_UNIT_TEST_SCAN_V00)
    if(this->dataHistogramV00)
    {
      free(this->dataHistogramV00);
      this->dataHistogramV00 = NULL;
    }
    if(this->dataPastHistogramV00)
    {
      free(this->dataPastHistogramV00);
      this->dataPastHistogramV00 = NULL;
    }
#endif
    /* *** */
#if (ENABLE_UNIT_TEST_SCAN_V01)
    if(this->dataHistogramV01)
    {
      free(this->dataHistogramV01);
      this->dataHistogramV01 = NULL;
    }
    if(this->dataPastHistogramV01)
    {
      free(this->dataPastHistogramV01);
      this->dataPastHistogramV01 = NULL;
    }
#endif
    /* *** */
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
#if SCAN_DEBUG_ENABLE
  void 
  invalidateDebug
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates error data on the host.
  */
#if SCAN_ERROR_TRACK_ENABLE
  void 
  invalidateError
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates unit test data on the host.
  */
#if (ENABLE_UNIT_TEST_SCAN_V00)
  void 
  invalidateUnitTestData_v00
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates unit test data on the host.
  */
#if (ENABLE_UNIT_TEST_SCAN_V01)
  void 
  invalidateUnitTestData_v01
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation.
  */
#if SCAN_ENABLE_V00
  void 
  scan_v00
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint,
    cl_uint,
    cl::Buffer&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation.
  */
#if SCAN_ENABLE_V01
  void 
  scan_v01
  (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
    cl_uint,
#endif
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl::Buffer&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation as a unit test
  */
#if (ENABLE_UNIT_TEST_SCAN_V00)
  void 
  scanUnitTest_v00
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation as a unit test
  */
#if (ENABLE_UNIT_TEST_SCAN_V01)
  void 
  scanUnitTest_v01
  (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
    cl_uint,
#endif
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
    cl::CommandQueue&,
    cl::Event&,
    cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies scan operation.
  */
#if SCAN_VERIFY_ENABLE
  void 
  verifyScan_v00
  (
    cl::CommandQueue&,
    void*,
    cl_uint(*)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
    cl_uint(*)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies scan operation.
  */
#if SCAN_VERIFY_ENABLE
  void 
  verifyScan_v01
  (
    cl::CommandQueue&,
    void*,
    cl_uint(*)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
    cl_uint(*)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
    cl_uint,
    cl_uint,
    cl_uint
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
#if (SCAN_DEBUG_ENABLE)
    size_t,
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
    size_t,
#endif
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE)
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
    Does unit test setup. Must be called only from the constructor and only once.
  */
#if (ENABLE_UNIT_TEST_SCAN_V00)
  void
  setupUnitTest_v00
  (
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    struct kernelStatistics&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Does unit test setup. Must be called only from the constructor and only once.
  */
#if (ENABLE_UNIT_TEST_SCAN_V01)
  void
  setupUnitTest_v01
  (
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    struct kernelStatistics&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes unit test data.
  */
#if (ENABLE_UNIT_TEST_SCAN_V00)
  void
  setUnitTestData_v00
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes unit test data.
  */
#if (ENABLE_UNIT_TEST_SCAN_V01)
  void
  setUnitTestData_v01
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
#if (ENABLE_UNIT_TEST_SCAN_V00)
  static cl_uint
  getPreviousHistogramItem_v00
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
#if (ENABLE_UNIT_TEST_SCAN_V00)
  static cl_uint
  getCurrentHistogramItem_v00
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
#if (ENABLE_UNIT_TEST_SCAN_V01)
  static cl_uint
  getPreviousHistogramItem_v01
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
#if (ENABLE_UNIT_TEST_SCAN_V01)
  static cl_uint
  getCurrentHistogramItem_v01
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads data from the device if it is invalid on the host.
  */
#if SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || \
    ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01
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
#if SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || \
    ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01
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

#endif  /*ENABLE_OPERATOR_SCAN*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif  /*OPERATOR_SCAN_H_*/
