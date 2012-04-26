/* ===============================================================================================



  =============================================================================================== */

#ifndef INTEGRATION_TEST_H_
#define INTEGRATION_TEST_H_



#include "Definitions.hpp"



/**
* IntegrationTest 
* Class implements OpenCL Gaussian Noise sample
* Derived from SDKSample base class
*/

class IntegrationTest : public SDKSample
{
    cl::Context context;                    /**< Context */
    vector<cl::Device> devices;        /**< vector of devices */
    vector<cl::Device> device;         /**< device to be used */
    vector<cl::Platform> platforms;    /**< vector of platforms */

    cl::CommandQueue commandQueue;          /**< command queue */
    
    kernelStatistics kernelStats;                 /**<Storage for memory stats*/
    cl_double setupTime;                    /**< time taken to setup OpenCL resources and building kernel */
    cl_double kernelTime;                   /**< time taken to run kernel and read result back */
    int iterations;                         /**< Number of iterations for kernel execution */
    int warmupcount;
    cl_uint currentTimeStep;
    cl_uint currentTimeSlot;
    
    cl_uint* dummyPointer;
    
    /*Model variables and parameters*/
    neuron_iz_ps     *nrn_ps;
    int              *ne;
    DATA_TYPE        *te_ps;
    DATA_TYPE        **co;
    DATA_TYPE        tau_ampa_ps;
    DATA_TYPE        tau_gaba_ps;
    DATA_TYPE        co_g_ampa_ps;
    DATA_TYPE        co_g_gaba_ps;
    DATA_TYPE        E_ps;
    DATA_TYPE        a_ps;
    DATA_TYPE        dt_ps;
    DATA_TYPE        tol_ps;
    int              steps_ps;
    
    std::string startTimeStamp;
    unsigned int srandSeed, srandCounter;
    
#ifdef WIN32
  __int64 current, performanceFrequency;
#endif
    
#if (LOG_MODEL_VARIABLES)
    std::ofstream     *traceFile;
    std::stringstream *dataToTraceFile;
#endif

    std::ofstream     *simulationLogFile;
    std::stringstream *dataToSimulationLogFile;

    std::ofstream     *reportLogFile;
    std::stringstream *dataToReportLogFile;

#if (SIMULATION_SNAPSHOT)
    std::ofstream     *snapshotLogFile;
    std::stringstream *dataToSnapshotLogFile;
#endif

#if STATISTICS_ENABLE
    double averageEventsInNetwork;
    double averageEventsInNetworkCounter;
    double averageSpikesInNetwork;
    double averageSpikesInNetworkCounter;
#endif

/**************************************************************************************************/
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
    cl_uint dataHistogramGroupEventsTikSize;
    cl_uint dataHistogramGroupEventsTikSizeBytes;
    cl_uint* dataHistogramGroupEventsTik;
    cl::Buffer dataHistogramGroupEventsTikBuffer;
    cl_uint dataHistogramGroupEventsVerifySize;
    cl_uint dataHistogramGroupEventsVerifySizeBytes;
    cl_uint* dataHistogramGroupEventsVerify;
    cl_uint lmGroupEventsHistogramOutSizeBytes;
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
    cl_uint dataUnsortedEventTargetsSize;
    cl_uint dataUnsortedEventTargetsSizeBytes;
    cl_uint* dataUnsortedEventTargets;
    cl::Buffer dataUnsortedEventTargetsBuffer;
    
    cl_uint dataUnsortedEventDelaysSize;
    cl_uint dataUnsortedEventDelaysSizeBytes;
    cl_uint* dataUnsortedEventDelays;
    cl::Buffer dataUnsortedEventDelaysBuffer;
    
    cl_uint dataUnsortedEventWeightsSize;
    cl_uint dataUnsortedEventWeightsSizeBytes;
    cl_uint* dataUnsortedEventWeights;
    cl::Buffer dataUnsortedEventWeightsBuffer;
    
    cl_uint dataUnsortedEventCountsSize;
    cl_uint dataUnsortedEventCountsSizeBytes;
    cl_uint* dataUnsortedEventCounts;
    cl::Buffer dataUnsortedEventCountsBuffer;
#endif
/**************************************************************************************************/
#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
    cl_uint dataHistogramSize;
    cl_uint dataHistogramSizeBytes;
    cl_uint* dataHistogram;
    cl::Buffer dataHistogramBuffer;
    
    cl_uint dataHistogramVerifySize;
    cl_uint dataHistogramVerifySizeBytes;
    cl_uint* dataHistogramVerify;
#endif
/**************************************************************************************************/
#if SORT_VERIFY_ENABLE
    cl_uint dataUnsortedEventsSnapShotSize;
    cl_uint dataUnsortedEventsSnapShotSizeBytes;
    cl_uint* dataUnsortedEventsSnapShot;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    cl_uint dataGroupEventsTikSize;
    cl_uint dataGroupEventsTikVerifySize;
    
    cl_uint dataGroupEventsTikSizeBytes;
    cl_uint dataGroupEventsTikVerifySizeBytes;

    cl_uint* dataGroupEventsTik;
    cl_uint* dataGroupEventsTikVerify;

    cl::Buffer dataGroupEventsTikBuffer;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
    cl_uint lmCacheGroupEventsSizeBytes;
    cl_uint lmlocalHistogramReferenceSizeBytes;
#if (GROUP_EVENTS_DEBUG_ENABLE)
    cl_uint dataDebugHostGroupEventsSize;
    cl_uint dataDebugHostGroupEventsSizeBytes;
    cl_uint* dataDebugHostGroupEvents;
    cl::Buffer dataDebugHostGroupEventsBuffer;
    cl_uint dataDebugDeviceGroupEventsSize;
    cl_uint dataDebugDeviceGroupEventsSizeBytes;
    cl_uint* dataDebugDeviceGroupEvents;
    cl::Buffer dataDebugDeviceGroupEventsBuffer;
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    cl_uint dataErrorGroupEventsSize;
    cl_uint dataErrorGroupEventsSizeBytes;
    cl_uint* dataErrorGroupEvents;
    cl::Buffer dataErrorGroupEventsBuffer;
#endif
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    SpikeEvents spikeEvents;
    Connectome  connectome;
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
    cl_uint lmExpandEventsSizeBytes;

    cl_uint dataUnsortedEventsTargetsVerifySize;
    cl_uint dataUnsortedEventsTargetsVerifySizeBytes;
    cl_uint* dataUnsortedEventsTargetsVerify;
    
    cl_uint dataUnsortedEventsDelaysVerifySize;
    cl_uint dataUnsortedEventsDelaysVerifySizeBytes;
    cl_uint* dataUnsortedEventsDelaysVerify;
    
    cl_uint dataUnsortedEventsWeightsVerifySize;
    cl_uint dataUnsortedEventsWeightsVerifySizeBytes;
    cl_uint* dataUnsortedEventsWeightsVerify;

    cl_uint dataUnsortedEventCountsVerifySize;
    cl_uint dataUnsortedEventCountsVerifySizeBytes;
    cl_uint* dataUnsortedEventCountsVerify;

    size_t blockSizeX_KernelExpandEvents;
    size_t blockSizeY_KernelExpandEvents;
    cl::Kernel kernelExpandEvents;

#if (EXPAND_EVENTS_DEBUG_ENABLE)
    cl_uint dataExpandEventsDebugHostSize;
    cl_uint dataExpandEventsDebugHostSizeBytes;
    cl_uint* dataExpandEventsDebugHost;
    cl::Buffer dataExpandEventsDebugHostBuffer;
    cl_uint dataExpandEventsDebugDeviceSize;
    cl_uint dataExpandEventsDebugDeviceSizeBytes;
    cl_uint* dataExpandEventsDebugDevice;
    cl::Buffer dataExpandEventsDebugDeviceBuffer;
#endif

#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    cl_uint dataExpandEventsErrorSize;
    cl_uint dataExpandEventsErrorSizeBytes;
    cl_uint* dataExpandEventsError;
    cl::Buffer dataExpandEventsErrorBuffer;
#endif
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
    cl_uint lmScanDataSizeBytes;
#if (SCAN_DEBUG_ENABLE)
    cl_uint dataScanDebugHostSize;
    cl_uint dataScanDebugHostSizeBytes;
    cl_uint* dataScanDebugHost;
    cl::Buffer dataScanDebugHostBuffer;
    cl_uint dataScanDebugDeviceSize;
    cl_uint dataScanDebugDeviceSizeBytes;
    cl_uint* dataScanDebugDevice;
    cl::Buffer dataScanDebugDeviceBuffer;
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
    cl_uint dataScanErrorSize;
    cl_uint dataScanErrorSizeBytes;
    cl_uint* dataScanError;
    cl::Buffer dataScanErrorBuffer;
#endif
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V00
    cl::Kernel kernelScanHistogramV00;
    size_t blockSizeX_kernelScanHistogramV00;
    size_t blockSizeY_kernelScanHistogramV00;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00
    cl::Kernel kernelGroupEventsV00;
    size_t blockSizeX_kernelGroupEventsV00;
    size_t blockSizeY_kernelGroupEventsV00;
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V01
    cl::Kernel kernelScanHistogramV01;
    size_t blockSizeX_kernelScanHistogramV01;
    size_t blockSizeY_kernelScanHistogramV01;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
    cl_uint dataHistogramGroupEventsTokSize;

    cl_uint dataHistogramGroupEventsTokSizeBytes;

    cl_uint* dataHistogramGroupEventsTok;

    cl::Buffer dataHistogramGroupEventsTokBuffer;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
    cl_uint dataGroupEventsTokSize;
    cl_uint dataGroupEventsTokVerifySize;

    cl_uint dataGroupEventsTokSizeBytes;
    cl_uint dataGroupEventsTokVerifySizeBytes;

    cl_uint* dataGroupEventsTok;
    cl_uint* dataGroupEventsTokVerify;
    
    cl::Buffer dataGroupEventsTokBuffer;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01
    cl::Kernel kernelGroupEventsV01;
    size_t blockSizeX_kernelGroupEventsV01;
    size_t blockSizeY_kernelGroupEventsV01;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V02
    cl::Kernel kernelGroupEventsV02;
    size_t blockSizeX_kernelGroupEventsV02;
    size_t blockSizeY_kernelGroupEventsV02;
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V03
    cl::Kernel kernelGroupEventsV03;
    size_t blockSizeX_kernelGroupEventsV03;
    size_t blockSizeY_kernelGroupEventsV03;
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    cl_uint dataMakeEventPtrsStructSize;
    cl_uint dataMakeEventPtrsStructSizeBytes;
    cl_uint* dataMakeEventPtrsStruct;
    cl::Buffer dataMakeEventPtrsStructBuffer;
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
    cl_uint lmLastElementMakeEventPtrsSizeBytes;
    cl_uint lmGenericMakeEventPtrsSizeBytes;

    cl::Kernel kernelMakeEventPtrs;
    size_t blockSizeX_kernelMakeEventPtrs;
    size_t blockSizeY_kernelMakeEventPtrs;
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    cl_uint lmGenericGlueEventPtrsSizeBytes;
    cl::Kernel kernelGlueEventPtrs;
    size_t blockSizeX_kernelGlueEventPtrs;
    size_t blockSizeY_kernelGlueEventPtrs;
#endif
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    cl_uint dataMakeEventPtrsDebugHostSize;
    cl_uint dataMakeEventPtrsDebugHostSizeBytes;
    cl_uint* dataMakeEventPtrsDebugHost;
    cl::Buffer dataMakeEventPtrsDebugHostBuffer;
    cl_uint dataMakeEventPtrsDebugDeviceSize;
    cl_uint dataMakeEventPtrsDebugDeviceSizeBytes;
    cl_uint* dataMakeEventPtrsDebugDevice;
    cl::Buffer dataMakeEventPtrsDebugDeviceBuffer;
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    cl_uint dataMakeEventPtrsErrorSize;
    cl_uint dataMakeEventPtrsErrorSizeBytes;
    cl_uint* dataMakeEventPtrsError;
    cl::Buffer dataMakeEventPtrsErrorBuffer;
#endif
#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
    cl_uint lmGeneralPurposeSizeBytes;
    cl_uint lmSpikePacketsSizeBytes;
    
    cl_uint modelParametersSize;
    cl_uint modelParametersSizeBytes;
    cl_float *modelParameters;
    cl::Buffer modelParametersBuffer;
    
    cl_uint modelVariablesSize;
    cl_uint modelVariablesSizeBytes;
    cl_float *modelVariables;
    cl::Buffer modelVariablesBuffer;
    
    cl_uint constantCoefficientsSize;
    cl_uint constantCoefficientsSizeBytes;
    cl_float *constantCoefficients;
    cl::Buffer constantCoefficientsBuffer;

    cl::Kernel kernelUpdateNeuronsV00;
    cl::Kernel kernelUpdateSpikedNeuronsV00;
    size_t blockSizeX_kernelUpdateNeuronsV00;
    size_t blockSizeY_kernelUpdateNeuronsV00;
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    cl_uint psToleranceSize;
    cl_uint psToleranceSizeBytes;
    CL_DATA_TYPE *psTolerance;
    cl::Buffer psToleranceBuffer;
#endif
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    cl_uint dataUpdateNeuronsDebugHostSize;
    cl_uint dataUpdateNeuronsDebugHostSizeBytes;
    cl_uint *dataUpdateNeuronsDebugHost;
    cl::Buffer dataUpdateNeuronsDebugHostBuffer;
    cl_uint dataUpdateNeuronsDebugDeviceSize;
    cl_uint dataUpdateNeuronsDebugDeviceSizeBytes;
    cl_uint *dataUpdateNeuronsDebugDevice;
    cl::Buffer dataUpdateNeuronsDebugDeviceBuffer;
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    cl_uint dataUpdateNeuronsErrorSize;
    cl_uint dataUpdateNeuronsErrorSizeBytes;
    cl_uint *dataUpdateNeuronsError;
    cl::Buffer dataUpdateNeuronsErrorBuffer;
#endif
#endif
/**************************************************************************************************/

public:

    /**
    * Read bitmap image and allocate host memory
    * @param inputImageName name of the input file
    * @return 1 on success and 0 on failure
    */
    int allocateHostData();
    int registerLocalMemory();

    /** 
    * Constructor 
    * Initialize member variables
    * @param name name of sample (string)
    */
    IntegrationTest(std::string name) : SDKSample(name),
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
        dataHistogramGroupEventsTik(NULL),
        dataHistogramGroupEventsVerify(NULL),
#endif
/* *** */
#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
        dataUnsortedEventTargets(NULL),
        dataUnsortedEventDelays(NULL),
        dataUnsortedEventWeights(NULL),
        dataUnsortedEventCounts(NULL),
#endif
/* *** */
#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
        dataHistogram(NULL),
        dataHistogramVerify(NULL),
#endif
/* *** */
#if EXPAND_EVENTS_ENABLE
        dataUnsortedEventsTargetsVerify(NULL),
        dataUnsortedEventsDelaysVerify(NULL),
        dataUnsortedEventsWeightsVerify(NULL),
        dataUnsortedEventCountsVerify(NULL),
#if (EXPAND_EVENTS_DEBUG_ENABLE)
        dataExpandEventsDebugHost(NULL),
        dataExpandEventsDebugDevice(NULL),
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
        dataExpandEventsError(NULL),
#endif
#endif
/* *** */
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
#if (SCAN_DEBUG_ENABLE)
        dataScanDebugHost(NULL),
        dataScanDebugDevice(NULL),
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
        dataScanError(NULL),
#endif
#endif
/* *** */
#if SORT_VERIFY_ENABLE
        dataUnsortedEventsSnapShot(NULL),
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
        dataGroupEventsTik(NULL),
        dataGroupEventsTikVerify(NULL),
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
#if (GROUP_EVENTS_DEBUG_ENABLE)
        dataDebugHostGroupEvents(NULL),
        dataDebugDeviceGroupEvents(NULL),
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
        dataErrorGroupEvents(NULL),
#endif
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
        dataHistogramGroupEventsTok(NULL),
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
        dataGroupEventsTok(NULL),
        dataGroupEventsTokVerify(NULL),
#endif
/* *** */
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    dataMakeEventPtrsStruct(NULL),
#endif
/* *** */
#if MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    dataMakeEventPtrsDebugHost(NULL),
    dataMakeEventPtrsDebugDevice(NULL),
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    dataMakeEventPtrsError(NULL),
#endif
#endif
/* *** */
#if UPDATE_NEURONS_ENABLE_V00
    modelParameters(NULL),
    modelVariables(NULL),
    constantCoefficients(NULL),
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    psTolerance(NULL),
#endif
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    dataUpdateNeuronsDebugHost(NULL),
    dataUpdateNeuronsDebugDevice(NULL),
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    dataUpdateNeuronsError(NULL),
#endif
#endif
/* *** */
        dummyPointer(NULL)
    {
#if EXPAND_EVENTS_ENABLE
      blockSizeX_KernelExpandEvents = EXPAND_EVENTS_WG_SIZE_WI;
      blockSizeY_KernelExpandEvents = 1;
#endif
#if SCAN_ENABLE_V00
      blockSizeX_kernelScanHistogramV00 = SCAN_WG_SIZE_WI;
      blockSizeY_kernelScanHistogramV00 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V00
      blockSizeX_kernelGroupEventsV00 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV00 = 1;
#endif
#if SCAN_ENABLE_V01
      blockSizeX_kernelScanHistogramV01 = SCAN_WG_SIZE_WI;
      blockSizeY_kernelScanHistogramV01 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V01
      blockSizeX_kernelGroupEventsV01 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV01 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V02
      blockSizeX_kernelGroupEventsV02 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV02 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V03
      blockSizeX_kernelGroupEventsV03 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV03 = 1;
#endif
#if MAKE_EVENT_PTRS_ENABLE
      blockSizeX_kernelMakeEventPtrs = MAKE_EVENT_PTRS_WG_SIZE_WI;
      blockSizeY_kernelMakeEventPtrs = 1;
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    blockSizeX_kernelGlueEventPtrs = GLUE_EVENT_PTRS_WG_SIZE_WI;
    blockSizeY_kernelGlueEventPtrs = 1;
#endif
#endif
#if UPDATE_NEURONS_ENABLE_V00
      blockSizeX_kernelUpdateNeuronsV00 = UPDATE_NEURONS_WG_SIZE_WI_V00;
      blockSizeY_kernelUpdateNeuronsV00 = 1;
#endif
      iterations = 1;
      warmupcount = 0;
      currentTimeStep = 0;
      currentTimeSlot = 0;
    }

    /** 
    * Constructor 
    * Initialize member variables
    * @param name name of sample (const char*)
    */
    IntegrationTest(const char* name) : SDKSample(name),
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
        dataHistogramGroupEventsTik(NULL),
        dataHistogramGroupEventsVerify(NULL),
#endif
/* *** */
#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
        dataUnsortedEventTargets(NULL),
        dataUnsortedEventDelays(NULL),
        dataUnsortedEventWeights(NULL),
        dataUnsortedEventCounts(NULL),
#endif
/* *** */
#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
        dataHistogram(NULL),
        dataHistogramVerify(NULL),
#endif
/* *** */
#if EXPAND_EVENTS_ENABLE
        dataUnsortedEventsTargetsVerify(NULL),
        dataUnsortedEventsDelaysVerify(NULL),
        dataUnsortedEventsWeightsVerify(NULL),
        dataUnsortedEventCountsVerify(NULL),
#if (EXPAND_EVENTS_DEBUG_ENABLE)
        dataExpandEventsDebugHost(NULL),
        dataExpandEventsDebugDevice(NULL),
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
        dataExpandEventsError(NULL),
#endif
#endif
/* *** */
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
#if (SCAN_DEBUG_ENABLE)
        dataScanDebugHost(NULL),
        dataScanDebugDevice(NULL),
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
        dataScanError(NULL),
#endif
#endif
/* *** */
#if SORT_VERIFY_ENABLE
        dataUnsortedEventsSnapShot(NULL),
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
        dataGroupEventsTik(NULL),
        dataGroupEventsTikVerify(NULL),
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
#if (GROUP_EVENTS_DEBUG_ENABLE)
        dataDebugHostGroupEvents(NULL),
        dataDebugDeviceGroupEvents(NULL),
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
        dataErrorGroupEvents(NULL),
#endif
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
        dataHistogramGroupEventsTok(NULL),
#endif
/* *** */
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
        dataGroupEventsTok(NULL),
        dataGroupEventsTokVerify(NULL),
#endif
/* *** */
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    dataMakeEventPtrsStruct(NULL),
#endif
/* *** */
#if MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    dataMakeEventPtrsDebugHost(NULL),
    dataMakeEventPtrsDebugDevice(NULL),
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    dataMakeEventPtrsError(NULL),
#endif
#endif
/* *** */
#if UPDATE_NEURONS_ENABLE_V00
    modelParameters(NULL),
    modelVariables(NULL),
    constantCoefficients(NULL),
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    psTolerance(NULL),
#endif
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    dataUpdateNeuronsDebugHost(NULL),
    dataUpdateNeuronsDebugDevice(NULL),
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    dataUpdateNeuronsError(NULL),
#endif
#endif
/* *** */
        dummyPointer(NULL)
    {
#if EXPAND_EVENTS_ENABLE
      blockSizeX_KernelExpandEvents = EXPAND_EVENTS_WG_SIZE_WI;
      blockSizeY_KernelExpandEvents = 1;
#endif
#if SCAN_ENABLE_V00
      blockSizeX_kernelScanHistogramV00 = SCAN_WG_SIZE_WI;
      blockSizeY_kernelScanHistogramV00 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V00
      blockSizeX_kernelGroupEventsV00 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV00 = 1;
#endif
#if SCAN_ENABLE_V01
      blockSizeX_kernelScanHistogramV01 = SCAN_WG_SIZE_WI;
      blockSizeY_kernelScanHistogramV01 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V01
      blockSizeX_kernelGroupEventsV01 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV01 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V02
      blockSizeX_kernelGroupEventsV02 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV02 = 1;
#endif
#if GROUP_EVENTS_ENABLE_V03
      blockSizeX_kernelGroupEventsV03 = GROUP_EVENTS_WG_SIZE_WI;
      blockSizeY_kernelGroupEventsV03 = 1;
#endif
#if MAKE_EVENT_PTRS_ENABLE
      blockSizeX_kernelMakeEventPtrs = MAKE_EVENT_PTRS_WG_SIZE_WI;
      blockSizeY_kernelMakeEventPtrs = 1;
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    blockSizeX_kernelGlueEventPtrs = GLUE_EVENT_PTRS_WG_SIZE_WI;
    blockSizeY_kernelGlueEventPtrs = 1;
#endif
#endif
#if UPDATE_NEURONS_ENABLE_V00
      blockSizeX_kernelUpdateNeuronsV00 = UPDATE_NEURONS_WG_SIZE_WI_V00;
      blockSizeY_kernelUpdateNeuronsV00 = 1;
#endif
      iterations = 1;
      warmupcount = 0;
      currentTimeStep = 0;
      currentTimeSlot = 0;
    }

    ~IntegrationTest()
    {
    }

    /**
     * Override from SDKSample, Generate binary image of given kernel 
     * and exit application
     */
    int genBinaryImage();

    /**
    * OpenCL related initialisations. 
    * Set up Context, Device list, Command Queue, Memory buffers
    * Build CL kernel program executable
    * @return 1 on success and 0 on failure
    */
    int setupCL();
    
    int
    createKernel
    (
      cl::Kernel    &kernel,
      std::string   kernelFileName,
      const char    *kernelName,
      std::string   flagsStr,
      size_t        blockSizeX,
      size_t        blockSizeY
    );
    
    
    bool 
    findTargetDevice
    (
      vector<cl::Device>                  platformDevices,
      const char                          *targetDevices,
      std::vector<cl::Device>::iterator   *d
    );

    /**
    * Set values for kernels' arguments, enqueue calls to the kernels
    * on to the command queue, wait till end of kernel execution.
    * Get kernel start and end time if timing is enabled
    * @return 1 on success and 0 on failure
    */
    int runCLKernels();

    /**
    * Override from SDKSample. Print sample stats.
    */
    void printStats();

    /**
    * Override from SDKSample. Initialize 
    * command line parser, add custom options
    */
    int initialize();

    /**
    * Override from SDKSample, adjust width and height 
    * of execution domain, perform all sample setup
    */
    int setup();

    /**
    * Override from SDKSample
    * Run OpenCL Sobel Filter
    */
    int run();

    /**
    * Override from SDKSample
    * Cleanup memory allocations
    */
    int cleanup();

    /**
    * Override from SDKSample
    * Verify against reference implementation
    */
    int verifyResults();
    
    /**
    * TODO: Add description
    */
    int getPlatformStats();
    
    /**
    * Verification/initialization methods for some kernel tests
    */

    int
    initializeGrouppedEvents
    (
      cl_uint enableValues,
      cl_uint totalBuffers,
      cl_uint bufferSize,
      cl_uint destinationBufferSize,
      cl_uint elementSizeWords,
      cl_uint histogramBinSize,
      cl_uint histogramTotalBins,
      cl_uint histogramBitMask,
      cl_uint histogramBitShift,
      cl_uint histogramOutBitShift,
      cl_uint histogramOutTotalGroups,
      cl_uint histogramOutBitMask,
      cl_uint histogramOutSize,
      cl_uint keyOffset,
      cl_uint *dataUnsortedEventCounts,
      cl_uint *dataUnsortedEventTargets,
      cl_uint *dataUnsortedEventDelays,
      cl_uint *dataUnsortedEventWeights,
      cl_uint *dataGroupedEvents,
      cl_uint *dataHistogram,
      cl_uint *dataHistogramOut
    );

    int
    initializeUnsortedEvents
    (
      cl_uint totalBuffers,
      cl_uint bufferSize,
      cl_uint totalNeurons,
      cl_uint maxDelay,
      cl_float minDelay,
      cl_uint histogramBinSize,
      cl_uint histogramBitMask,
      cl_uint histogramBitShift,
      cl_uint keyOffset,
      cl_float percentInhibitory,
      double   percentBufferSizeDeviation,
      cl_uint *dataUnsortedEventCounts,
      cl_uint *dataUnsortedEventTargets,
      cl_uint *dataUnsortedEventDelays,
      cl_uint *dataUnsortedEventWeights,
      cl_uint *dataHistogram
    );
    
    int 
    verifyKernelGroupEvents
    (
      cl_uint verifyKeyBins,
      cl_uint timeSlot,
      cl_uint step,
      cl_uint value1CarryEnable,
      cl_uint value2CarryEnable,
      cl_uint stepShiftEnable,
      cl_uint elementSizeWords,
      cl_uint destinationBufferSize,
      cl_uint histogramBinSize,
      cl_uint histogramTotalBins,
      cl_uint histogramBitMask,
      cl_uint histogramBitShift,
      cl_uint histogramOutTotalBins,
      cl_uint histogramOutTotalGroups,
      cl_uint groupSize,
      cl_uint *dataHistogram,
      cl_uint *dataHistogramOut,
      cl_uint *dataHistogramOutVerify,
      cl_uint *dataGroupedEvents,
      cl_uint *dataGroupedEventsVerify
    );
    
    int 
    captureUnsortedEvents
    (
      cl_uint *unsortedEventCounts,
      cl_uint *unsortedEventTargets,
      cl_uint *unsortedEventDelays,
      cl_uint *unsortedEventWeights
    );
    
    int 
    verifySortedEvents
    (
      cl_uint *sortedEvents, 
      cl_uint *pointerStruct, 
      cl_uint level
    );
    
    int initializeExpandEventsData();
    
    int verifyKernelExpandEvents
    (
      cl_uint,
      bool,
      SpikeEvents&,
      Connectome&,
      cl::CommandQueue&
    );
    
    int initializeDataForKernelScanHistogramV00(cl_uint timeSlot, cl_uint mode);
    int verifyKernelScanHistogramV00(cl_uint timeSlot);
    
    int initializeDataForKernelScanHistogramV01();
    int verifyKernelScanHistogramV01();
    
    int 
    initializeDataForKernelGroupEventsV00
    (
      cl_uint   unsortedEventTimeSlotDelta,
      double  percentTimeSlotDeltaDeviation,
      double  percentInh
    );
    
    int 
    initializeEventBuffers
    (
      cl_uint   timeSlots,
      cl_uint   buffers,
      cl_uint   bufferSize,
      cl_uint   neurons,
      cl_uint   unsortedEventTimeSlotDelta,
      cl_uint   histogramBitShift,
      cl_uint   histogramBinMask,
      cl_uint   histogramTotalBins,
      cl_uint   histogramBinSize,
      double    minDelay,
      double    percentDelayBoarder,
      double    percentTimeSlotDeltaDeviation,
      double    percentInh,
      cl_uint   *dataUnsortedEventCounts,
      cl_uint   *dataUnsortedEventTargets,
      cl_uint   *dataUnsortedEventDelays,
      cl_uint   *dataUnsortedEventWeights,
      cl_uint   *dataHistogram
    );
    
    int 
    initializeHistogram
    (
      cl_uint   timeSlots,
      cl_uint   histogramTotalBins,
      cl_uint   histogramBinSize,
      cl_uint   eventDestinationBufferSize,
      cl_uint   *dataHistogram
    );
    
    int verifyKernelGroupEventsV00(cl_uint keyOffset);
    
    int 
    initializeDataForKernelGroupEventsV01
    (
      int step, 
      cl_uint keyOffset,
      double percentBufferSizeDeviation
    );

    int verifyKernelGroupEventsV01(cl_uint step);

    int
    initializeDataForKernelGroupEventsV02_V03
    (
      cl_uint step, 
      cl_uint keyOffset,
      double percentBufferSizeDeviation
    );
    int verifyKernelGroupEventsV02(cl_uint step, cl_uint keyOffset);
    int verifyKernelGroupEventsV03(cl_uint step);
    
    int 
    initializeSortedEvents
    (
      cl_uint mode,
      cl_uint maxEventId,
      double eventStdDev,
      double structBufferSizeMargin,
      double gabaRatio,
      cl_uint totalEvents,
      cl_uint wfWorkSize,
      cl_uint structSize,
      cl_uint structPitch,
      cl_uint eventBufferSize,
      cl_uint *sortedEvents
    );
    
    int 
    initializeEventPointers
    (
      bool    verify,
      cl_uint totalSortedEvents,
      cl_uint totalPointers,
      cl_uint pointersPitch,
      cl_uint *sortedEvents,
      cl_uint *pointersToEvents
    );
    
    int
    initializeDataForKernelMakeEventPtrs
    (
      cl_uint mode,
      cl_uint step
    );
    
    int verifyKernelMakeEventPtrs();
    
    int
    psInit
    (
      cl_uint     totalNeurons,
      int         injectCurrentUntilStep,
      const char  *neuronVariablesSampleFile
    );

    void psClean();
    
    int 
    initializeDataForKernelUpdateNeurons
    (
      bool          resetEvents,
      bool          resetParameters,
      bool          resetVariables,
      double        gabaRatio,
      const char    *neuronVariablesSampleFile
    );

    int 
    propagateSpikes
    (
      unsigned int,
      unsigned int,
      neuron_iz_ps*,
      int*,
      DATA_TYPE*,
      Connectome&,
      cl::CommandQueue&
    );
    
    int 
    injectSortedEvents
    (
      bool          verify,
      bool          resetEventBuffer,
      cl_uint       totalNeurons,
      cl_uint       eventQueueSize,
      cl_uint       pointersPitch,
      cl_uint       sortedEventsSize,
      cl_uint       *sortedEvents,
      cl_uint       *pointersToEvents,
      neuron_iz_ps  *nrn
    );
    
    int 
    injectUnsortedEvents
    (
      cl_uint       timeSlots,
      cl_uint       buffers,
      cl_uint       bufferSize,
      cl_uint       totalNeurons,
      cl_uint       eventQueueSize,
      cl_uint       *dataUnsortedEventCounts,
      cl_uint       *dataUnsortedEventTargets,
      cl_uint       *dataUnsortedEventWeights,
      cl_uint       *dataUnsortedEventDelays,
      neuron_iz_ps  *nrn
    );
    
    int 
    stepIzPs
    (
      DATA_TYPE **yp, 
      DATA_TYPE **co, 
      DATA_TYPE *yold, 
      DATA_TYPE *ynew, 
      neuron_iz_ps *nrnp, 
      int *ne, 
      DATA_TYPE *te, 
      int *ip, 
      DATA_TYPE *fp, 
#if STATISTICS_ENABLE
      int *fcount,
      unsigned long long int *icount, 
      DATA_TYPE *mu, 
      int *max_order, 
#endif
      int ps_order_limit,
      int nr_order_limit,
      DATA_TYPE nrTolerance,
      bool variableDelaysEnalbe
    );
    
    int 
    updateStep
    (
      bool    ignoreFailures,
      int     injectCurrentUntilStep,
      cl_uint currentTimeStep,
      cl_uint totalNeurons,
      cl_uint psOrderLimit,
      cl_uint nrOrderLimit,
      double  nrTolerance
    );
    
    int
    verifyKernelUpdateNeurons
    (
      bool,
      bool,
      bool,
      cl_uint,
      cl_uint,
      unsigned int*,
      unsigned int*,
      DATA_TYPE*,
      SpikeEvents&,
      Connectome&,
      cl::CommandQueue&
    );
    
    int 
    verifyEvents
    (
      bool          ignoreWarnings,
      bool          correctWeightPositionMismatch,
      unsigned int  totalNeurons,
      unsigned int  structElementSize,
      unsigned int  sortedEventsSize,
      unsigned int  *pointerStruct,
      unsigned int  *sortedEvents,
      neuron_iz_ps  *nrn
    );
    
    int 
    initializeSpikeData
    (
      double spikeBufferMinPercentFill,
      double spikeBufferMaxPercentFill,
      double spikeNeuronsPercent
    );
    
    double timeStampNs();

#if SIMULATION_SNAPSHOT
    int 
    takeSimulationSnapshot
    (
      cl_uint,
      cl_uint,
      cl_uint*,
      cl_uint*,
      cl_uint*,
      cl_uint*,
      cl_float*,
      SpikeEvents&,
      cl::CommandQueue&
    );
#endif
};

#endif // INTEGRATION_TEST_H_
