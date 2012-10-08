
#ifndef NEUROSIM_H_
#define NEUROSIM_H_



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
  Forward declarations for cyclic dependency
***************************************************************************************************/

class Data_Connectome;
class Data_SpikeEvents;
class Data_SynapticEvents;
class Operator_Expand;
class Operator_Scan;
class Operator_Sort;

/**************************************************************************************************/



/**
  @class Neurosim

  Implements SNN simulation.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2011/09/17
 */
class Neurosim
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/



/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/

  /*OpenCL definitions*/
  cl::Context context;
  vector<cl::Device> devices;
  vector<cl::Platform> platforms;
  cl::CommandQueue commandQueue;
  
  /*Simulation control*/
  cl_uint currentTimeStep;
  cl_uint currentTimeSlot;
  
  /*Control of random number generation*/
  unsigned int srandSeed, srandCounter;

  /*Statistics, timing, logging*/
  std::string startTimeStamp;
  kernelStatistics kernelStats;
  /* *** */
#ifdef WIN32
  __int64 performanceCounter, performanceFrequency;
#endif
  /* *** */
#if (LOG_MODEL_VARIABLES)
  std::ofstream     traceFile;
  std::stringstream dataToTraceFile;
#endif
  /* *** */
  std::ofstream     simulationLogFile;
  std::stringstream *dataToSimulationLogFile;
  /* *** */
  std::ofstream     reportLogFile;
  std::stringstream *dataToReportLogFile;
  /* *** */
#if (SIMULATION_SNAPSHOT)
  std::ofstream     snapshotLogFile;
  std::stringstream dataToSnapshotLogFile;
#endif
  /* *** */
#if STATISTICS_ENABLE
  double averageEventsInNetwork;
  double averageEventsInNetworkCounter;
  double averageSpikesInNetwork;
  double averageSpikesInNetworkCounter;
  double setupTime;
  double runTime;
#endif

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

  /*Simulation objects*/
#if ENABLE_OPERATOR_SCAN
  Operator_Scan            *operatorScan;
#endif
  /* *** */
#if ENABLE_OPERATOR_SORT
  Operator_Sort            *operatorSort;
#endif
  /* *** */
#if ENABLE_OPERATOR_EXPAND
  Operator_Expand          *operatorExpand;
#endif
  /* *** */
  Data_SpikeEvents         *spikeEvents;
  Data_Connectome          *connectome;
  Data_SynapticEvents      *synapticEvents;
  
/**************************************************************************************************/
#if ((GROUP_EVENTS_ENABLE_V00) ||\
   (GROUP_EVENTS_ENABLE_V01) ||\
   (GROUP_EVENTS_ENABLE_V02) ||\
   (GROUP_EVENTS_ENABLE_V03))
  cl_uint dataHistogramGroupEventsVerifySize;
  cl_uint dataHistogramGroupEventsVerifySizeBytes;
  cl_uint* dataHistogramGroupEventsVerify;
  cl_uint lmGroupEventsHistogramOutSizeBytes;
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
#if GROUP_EVENTS_ENABLE_V00
  cl::Kernel kernelGroupEventsV00;
  size_t blockSizeX_kernelGroupEventsV00;
  size_t blockSizeY_kernelGroupEventsV00;
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



/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Constructor.
  */
  Neurosim() :
#if ((GROUP_EVENTS_ENABLE_V00) ||\
   (GROUP_EVENTS_ENABLE_V01) ||\
   (GROUP_EVENTS_ENABLE_V02) ||\
   (GROUP_EVENTS_ENABLE_V03))
    dataHistogramGroupEventsVerify(NULL),
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
#if (LOG_MODEL_VARIABLES)
    traceFile(LOG_MODEL_VARIABLES_FILE_NAME),
    dataToTraceFile("", std::ios::out | std::ios::app),
#endif
    /* *** */
#if (SIMULATION_SNAPSHOT)
    snapshotLogFile(LOG_SNAPSHOT_FILE_NAME),
    dataToSnapshotLogFile("", std::ios::out | std::ios::app),
#endif
    /* *** */
#if (LOG_SIMULATION)
    simulationLogFile(LOG_SIMULATION_FILE_NAME, std::ofstream::out | std::ofstream::app),
    dataToSimulationLogFile(new std::stringstream("", std::ios::out | std::ios::app)),
#else
    dataToSimulationLogFile(NULL),
#endif
    /* *** */
#if (LOG_REPORT)
    reportLogFile(LOG_REPORT_FILE_NAME),
    dataToReportLogFile(new std::stringstream("", std::ios::out | std::ios::app)),
#else
    dataToReportLogFile(NULL),
#endif
    /* *** */
    nrn_ps(NULL),
    ne(NULL),
    te_ps(NULL),
    co(NULL),
    /* *** */
#if GROUP_EVENTS_ENABLE_V00
    blockSizeX_kernelGroupEventsV00(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV00(1),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V01
    blockSizeX_kernelGroupEventsV01(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV01(1),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V02
    blockSizeX_kernelGroupEventsV02(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV02(1),
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V03
    blockSizeX_kernelGroupEventsV03(GROUP_EVENTS_WG_SIZE_WI),
    blockSizeY_kernelGroupEventsV03(1),
#endif
    /* *** */
#if MAKE_EVENT_PTRS_ENABLE
    blockSizeX_kernelMakeEventPtrs(MAKE_EVENT_PTRS_WG_SIZE_WI),
    blockSizeY_kernelMakeEventPtrs(1),
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    blockSizeX_kernelGlueEventPtrs(GLUE_EVENT_PTRS_WG_SIZE_WI),
    blockSizeY_kernelGlueEventPtrs(1),
#endif
#endif
    /* *** */
#if UPDATE_NEURONS_ENABLE_V00
    blockSizeX_kernelUpdateNeuronsV00(UPDATE_NEURONS_WG_SIZE_WI_V00),
    blockSizeY_kernelUpdateNeuronsV00(1),
#endif
    /* *** */
#if STATISTICS_ENABLE
    averageEventsInNetwork(0),
    averageEventsInNetworkCounter(0),
    averageSpikesInNetwork(0),
    averageSpikesInNetworkCounter(0),
    setupTime(0),
    runTime(0),
#endif
    /* *** */
#ifdef WIN32
    performanceCounter(0),
    performanceFrequency(0),
#endif
    /* *** */
#if ENABLE_OPERATOR_SCAN
    operatorScan(NULL),
#endif
    /* *** */
#if ENABLE_OPERATOR_SORT
    operatorSort(NULL),
#endif
    /* *** */
#if ENABLE_OPERATOR_EXPAND
    operatorExpand(NULL),
#endif
    /* *** */
    currentTimeStep(0),
    currentTimeSlot(0),
    srandSeed(1),
    srandCounter(0),
    kernelStats(),
    spikeEvents(NULL),
    connectome(NULL),
    synapticEvents(NULL)
/**************************************************************************************************/
  {
    time_t rawtime; 
    time ( &rawtime ); 
    this->startTimeStamp = std::string(ctime (&rawtime));
    
    SET_RANDOM_SEED_DIRECT(this->srandSeed);
    
#ifdef WIN32
    QueryPerformanceFrequency((LARGE_INTEGER *)&(this->performanceFrequency));
#endif

#if (LOG_SIMULATION)
    if(!simulationLogFile.is_open())
    {
      throw SimException("Neurosim::Constructor: Unable to open simulation log file");
    }
    
    *dataToSimulationLogFile << "-------\n" << startTimeStamp;
#endif

#if (LOG_REPORT)
    if(!reportLogFile.is_open())
    {
      throw SimException("Neurosim::Constructor: Unable to open report log file");
    }
#endif

#if (SIMULATION_SNAPSHOT)
    if(!snapshotLogFile.is_open())
    {
      throw SimException("Neurosim::Constructor: Unable to open snapshot log file");
    }
#endif

#if (LOG_MODEL_VARIABLES)
    if(!traceFile.is_open())
    {
      throw SimException("Neurosim::Constructor: Unable to open trace log file");
    }
#endif
  };
/**************************************************************************************************/
  
  
  
/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Neurosim()
/**************************************************************************************************/
  {
#if ((GROUP_EVENTS_ENABLE_V00) ||\
     (GROUP_EVENTS_ENABLE_V01) ||\
     (GROUP_EVENTS_ENABLE_V02) ||\
     (GROUP_EVENTS_ENABLE_V03))
    if(dataHistogramGroupEventsVerify)
        free(dataHistogramGroupEventsVerify);
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    if(dataGroupEventsTik)
        free(dataGroupEventsTik);
    if(dataGroupEventsTikVerify)
        free(dataGroupEventsTikVerify);
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
#if (GROUP_EVENTS_DEBUG_ENABLE)
    if(dataDebugHostGroupEvents)
        free(dataDebugHostGroupEvents);
    if(dataDebugDeviceGroupEvents)
        free(dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(dataErrorGroupEvents)
        free(dataErrorGroupEvents);
#endif
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
    if(dataHistogramGroupEventsTok)
        free(dataHistogramGroupEventsTok);
#endif
    /* *** */
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
    if(dataGroupEventsTok)
        free(dataGroupEventsTok);
    if(dataGroupEventsTokVerify)
        free(dataGroupEventsTokVerify);
#endif
    /* *** */
#if SORT_VERIFY_ENABLE
    if(dataUnsortedEventsSnapShot)
        free(dataUnsortedEventsSnapShot);
#endif
    /* *** */
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    if(dataMakeEventPtrsStruct)
        free(dataMakeEventPtrsStruct);
#endif
    /* *** */
#if MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    if(dataMakeEventPtrsDebugHost)
        free(dataMakeEventPtrsDebugHost);
    if(dataMakeEventPtrsDebugDevice)
        free(dataMakeEventPtrsDebugDevice);
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    if(dataMakeEventPtrsError)
        free(dataMakeEventPtrsError);
#endif
#endif
    /* *** */
#if UPDATE_NEURONS_ENABLE_V00
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    if(dataUpdateNeuronsDebugHost)
        free(dataUpdateNeuronsDebugHost);
    if(dataUpdateNeuronsDebugDevice)
        free(dataUpdateNeuronsDebugDevice);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(dataUpdateNeuronsError)
        free(dataUpdateNeuronsError);
#endif
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    if(psTolerance)
        free(psTolerance);
#endif
    if(modelParameters)
        free(modelParameters);
    if(modelVariables)
        free(modelVariables);
    if(constantCoefficients)
        free(constantCoefficients);
    psClean();
#endif
    /* *** */
#if (LOG_MODEL_VARIABLES)
    traceFile << dataToTraceFile.str();
    traceFile.close();
#endif
    /* *** */
#if (LOG_SIMULATION)
    simulationLogFile << (*dataToSimulationLogFile).str();
    simulationLogFile.close();
    delete dataToSimulationLogFile;
#endif
    /* *** */
#if (LOG_REPORT)
    reportLogFile << (*dataToReportLogFile).str();
    reportLogFile.close();
    delete dataToReportLogFile;
#endif
    /* *** */
#if (SIMULATION_SNAPSHOT)
    snapshotLogFile.close();
#endif
    /* *** */
#if ENABLE_OPERATOR_SCAN
    if(this->operatorScan)
    {
      delete(this->operatorScan);
    }
#endif
    /* *** */
#if ENABLE_OPERATOR_SORT
    if(this->operatorSort)
    {
      delete(this->operatorSort);
    }
#endif
    /* *** */
#if ENABLE_OPERATOR_EXPAND
    if(this->operatorExpand)
    {
      delete(this->operatorExpand);
    }
#endif
    /* *** */
    if(this->connectome)
    {
      delete(this->connectome);
    }
    
    if(this->synapticEvents)
    {
      delete(this->synapticEvents);
    }
    
    if(this->spikeEvents)
    {
      delete(this->spikeEvents);
    }
    /* *** */
  };
/**************************************************************************************************/

  
  
/**************************************************************************************************/
  /**
    Main execution call.
  */
  int execute();
/**************************************************************************************************/
  
  
  
/**************************************************************************************************/
  private:  /*private methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Copy constructor: non-copiable, C++98 style.
  */
  Neurosim(const Neurosim&);
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Assignment operator: non-assignable, C++98 style.
  */
  Neurosim& operator = (const Neurosim&);
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Call for execution setup.
  */
  void 
  setup
  ();
/**************************************************************************************************/
  
  
  
/**************************************************************************************************/
  /**
    Call to run simulation.
  */
  void 
  run
  ();
/**************************************************************************************************/

  
  
/**************************************************************************************************/
  /**
    Log platform statistics.
  */
  void 
  getPlatformStats
  ();
/**************************************************************************************************/
  
  
  
/**************************************************************************************************/
  /**
    Allocates data structures.
  */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  void
  allocateHostData
  (
    cl::Device&
  );
#endif
/**************************************************************************************************/
  
  
  
/**************************************************************************************************/
  /**
    Logs usage of local memory.
  */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  void 
  registerLocalMemory
  (
    cl::Device&
  );
#endif
/**************************************************************************************************/
  
  

/**************************************************************************************************/
  /**
    Logs mem allocations
  */
  void 
  getMemoryUsageStats
  ();
/**************************************************************************************************/



#if (STATISTICS_ENABLE)
/**************************************************************************************************/
  /**
    Logs some stats.
  */
  void 
  printStats
  ();
/**************************************************************************************************/
#endif


  
  /**
  * Verification/initialization methods for some kernel tests
  */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  int
  initializeGrouppedEvents
  (
    cl_uint enableValues,
    cl_uint totalBuffers,
    cl_uint bufferSize,
    cl_uint destinationBufferSize,
    cl_uint histogramBinSize,
    cl_uint histogramTotalBins,
    cl_uint histogramBitMask,
    cl_uint histogramBitShift,
    cl_uint histogramOutBitShift,
    cl_uint histogramOutTotalGroups,
    cl_uint histogramOutBitMask,
    size_t  histogramOutSize,
    cl_uint keyOffset,
    cl_uint *dataUnsortedEventCounts,
    cl_uint *dataUnsortedEventTargets,
    cl_uint *dataUnsortedEventDelays,
    cl_uint *dataUnsortedEventWeights,
    cl_uint *dataGroupedEvents,
    cl_uint *dataHistogram,
    cl_uint *dataHistogramOut
  );
#endif

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
  
#if SORT_VERIFY_ENABLE
  int 
  captureUnsortedEvents
  (
    cl_uint *unsortedEventCounts,
    cl_uint *unsortedEventTargets,
    cl_uint *unsortedEventDelays,
    cl_uint *unsortedEventWeights
  );
#endif
  
#if SORT_VERIFY_ENABLE
  int 
  verifySortedEvents
  (
    cl_uint *sortedEvents, 
    cl_uint *pointerStruct, 
    cl_uint level
  );
#endif
  
  int 
  verifyKernelGroupEventsV00
  ();
  
#if GROUP_EVENTS_ENABLE_V01
  int 
  initializeDataForKernelGroupEventsV01
  (
    int step, 
    cl_uint keyOffset,
    double percentBufferSizeDeviation
  );
#endif

#if GROUP_EVENTS_ENABLE_V01
  int verifyKernelGroupEventsV01
  (
    cl_uint step
  );
#endif

#if GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  int
  initializeDataForKernelGroupEventsV02_V03
  (
    int,
    cl_uint,
    double
  );
#endif
  
#if GROUP_EVENTS_ENABLE_V02
  int verifyKernelGroupEventsV02
  (
    cl_uint step, 
    cl_uint keyOffset
  );
#endif
  
#if GROUP_EVENTS_ENABLE_V03
  int verifyKernelGroupEventsV03
  (
    cl_uint step
  );
#endif
  
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
  
#if MAKE_EVENT_PTRS_ENABLE
  int
  initializeDataForKernelMakeEventPtrs
  (
    int mode,
    cl_uint step
  );
#endif
  
#if MAKE_EVENT_PTRS_ENABLE
  int 
  verifyKernelMakeEventPtrs
  ();
#endif
  
#if UPDATE_NEURONS_ENABLE_V00
  int
  psInit
  (
    cl_uint     totalNeurons,
    int         injectCurrentUntilStep,
    const char  *neuronVariablesSampleFile
  );
#endif

#if UPDATE_NEURONS_ENABLE_V00
  void 
  psClean
  ();
#endif
  
#if UPDATE_NEURONS_ENABLE_V00
  int 
  initializeDataForKernelUpdateNeurons
  (
    bool          resetEvents,
    bool          resetParameters,
    bool          resetVariables,
    double        gabaRatio,
    const char    *neuronVariablesSampleFile
  );
#endif

  int 
  propagateSpikes
  (
    unsigned int,
    unsigned int,
    neuron_iz_ps*,
    int*,
    DATA_TYPE*,
    Data_Connectome&
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
  
#if UPDATE_NEURONS_ENABLE_V00
  int 
  updateStep
  (
    bool    ignoreFailures,
    int     injectCurrentUntilStep,
    cl_uint currentTimeStep,
    cl_uint totalNeurons,
    int     psOrderLimit,
    int     nrOrderLimit,
    double  nrTolerance
  );
#endif
  
#if UPDATE_NEURONS_ENABLE_V00
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
    Data_SpikeEvents&,
    Data_Connectome&,
    cl::CommandQueue&
  );
#endif
  
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

#if SIMULATION_SNAPSHOT
  int 
  takeSimulationSnapshot
  (
    cl_uint,
    cl_uint,
    cl_uint*,
    cl_uint*,
    cl_float*,
    Data_SpikeEvents&,
    Data_SynapticEvents&,
    cl::CommandQueue&
  );
#endif
};



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif // NEUROSIM_H_
