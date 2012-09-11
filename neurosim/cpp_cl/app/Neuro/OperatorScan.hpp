
#ifndef OPERATOR_SCAN_H_
#define OPERATOR_SCAN_H_



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Common.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Class-scope Preprocessor Definitions
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



/**
  @class OperatorScan

  Prefix sum operation.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/09/01
 */
class OperatorScan 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/

#if (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)
  size_t lmScanDataSizeBytes;
  
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
#endif
  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/

  bool resetObject;
  cl_uint dataValid;
  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;
  
#if (SCAN_ENABLE_UNIT_TEST_V00 || SCAN_ENABLE_UNIT_TEST_V01)
  unsigned int srandSeed;
  unsigned int srandCounter;
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

#if (SCAN_ENABLE_UNIT_TEST_V00)
  bool initDataUnitTest_v00;
  cl_uint timeSlots;
  cl_uint histogramBinCount_v00;
  cl_uint histogramBinSize_v00;
  size_t dataHistogramV00Size;
  size_t dataHistogramV00SizeBytes;
  cl_uint* dataHistogramV00;
  cl_uint* dataPastHistogramV00;
  cl::Buffer dataHistogramV00Buffer;
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

#if (SCAN_ENABLE_UNIT_TEST_V01)
  bool initDataUnitTest_v01;
  cl_uint histogramBinBackets;
  cl_uint histogramBinCount_v01;
  cl_uint histogramBinSize_v01;
  size_t dataHistogramV01Size;
  size_t dataHistogramV01SizeBytes;
  cl_uint* dataHistogramV01;
  cl_uint* dataPastHistogramV01;
  cl::Buffer dataHistogramV01Buffer;
#endif



/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  OperatorScan
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~OperatorScan
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
    struct kernelStatistics*,
    std::stringstream*,
    std::stringstream*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates debug data on the host.
  */
#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_DEBUG_ENABLE)
  void 
  invalidateDebug
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates error data on the host.
  */
#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_ERROR_TRACK_ENABLE)
  void 
  invalidateError
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates histogram data on the host.
  */
#if (SCAN_ENABLE_UNIT_TEST_V00)
  void 
  invalidateUnitTestData_v00
  ();
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates histogram data on the host.
  */
#if (SCAN_ENABLE_UNIT_TEST_V01)
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
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl_uint,
    cl::Buffer&,
    cl::CommandQueue&,
    cl::Event&
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
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl::Buffer&,
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation as a unit test
  */
#if (SCAN_ENABLE_UNIT_TEST_V00)
  void 
  scanUnitTest_v00
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation as a unit test
  */
#if (SCAN_ENABLE_UNIT_TEST_V01)
  void 
  scanUnitTest_v01
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    struct kernelStatistics*,
#endif
    cl_uint,
    cl_uint,
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation.
  */
#if (SCAN_ENABLE_UNIT_TEST_V00)
  void
  setupUnitTest_v00
  (
    cl_uint,
    cl_uint,
    cl_uint,
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    struct kernelStatistics*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs scan operation.
  */
#if (SCAN_ENABLE_UNIT_TEST_V01)
  void
  setupUnitTest_v01
  (
    cl_uint,
    cl_uint,
    cl_uint,
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    struct kernelStatistics*
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies scan operation.
  */
#if (SCAN_ENABLE_V00 && SCAN_VERIFY_ENABLE)
  void 
  verifyScan_v00
  (
#if !(SCAN_ENABLE_UNIT_TEST_V00)
  void*,
  cl_uint(*)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
  cl_uint(*)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
#endif
  cl_uint,
  cl_uint,
  cl_uint,
  cl::CommandQueue&,
  cl_uint
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies scan operation.
  */
#if (SCAN_ENABLE_V01 && SCAN_VERIFY_ENABLE)
  void 
  verifyScan_v01
  (
#if !(SCAN_ENABLE_UNIT_TEST_V01)
    cl_uint*,
    cl_uint*,
#endif
    cl_uint,
    cl_uint,
    cl_uint,
    cl::CommandQueue&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  private:  /*private methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes histogram with random values.
  */
#if (SCAN_ENABLE_UNIT_TEST_V00)
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
    Initializes histogram with random values.
  */
#if (SCAN_ENABLE_UNIT_TEST_V01)
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
    Stores event data on device.
  */
  void
  storeData
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Saves past events and loads current events from the device if they are invalid on the host.
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
    Verifies if this object was initialized.
  */
  void
  isInitialized
  ();
/**************************************************************************************************/
};

#endif
