
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

#if ENABLE_OPERATOR_EXPAND
class Operator_Expand;
#endif
/* *** */
#if ENABLE_OPERATOR_SORT
class Operator_Sort;
#endif

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



/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  cl_uint dataMakeEventPtrsStructSize;
  cl_uint dataMakeEventPtrsStructSizeBytes;
  cl_uint* dataMakeEventPtrsStruct;
  cl::Buffer dataMakeEventPtrsStructBuffer;
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
  cl::Kernel kernelMakeEventPtrs;
  size_t blockSizeX_kernelMakeEventPtrs;
  size_t blockSizeY_kernelMakeEventPtrs;
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
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
#if ENABLE_OPERATOR_SORT
  Operator_Sort *operatorSort;
#endif
/**************************************************************************************************/
#if ENABLE_OPERATOR_EXPAND
  Operator_Expand *operatorExpand;
#endif
/**************************************************************************************************/
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
  __int64 performanceCounter, performanceFrequency;
#endif
/**************************************************************************************************/
#if STATISTICS_ENABLE
  double averageEventsInNetwork;
  double averageEventsInNetworkCounter;
  double averageSpikesInNetwork;
  double averageSpikesInNetworkCounter;
  double setupTime;
  double runTime;
#endif
/**************************************************************************************************/
#if (LOG_MODEL_VARIABLES)
  std::ofstream     traceFile;
  std::stringstream dataToTraceFile;
#endif
/**************************************************************************************************/
#if (SIMULATION_SNAPSHOT)
  std::ofstream     snapshotLogFile;
  std::stringstream dataToSnapshotLogFile;
#endif
/**************************************************************************************************/
#if (LOG_SIMULATION)
  std::ofstream     simulationLogFile;
  std::stringstream *dataToSimulationLogFile;
#else
  std::stringstream *dataToSimulationLogFile;
#endif
/**************************************************************************************************/
#if (LOG_REPORT)
  std::ofstream     reportLogFile;
  std::stringstream *dataToReportLogFile;
#else
  std::stringstream *dataToReportLogFile;
#endif
/**************************************************************************************************/
  Data_SpikeEvents         *spikeEvents;
  Data_Connectome          *connectome;
  Data_SynapticEvents      *synapticEvents;
/**************************************************************************************************/
  cl::Context context;
  cl::CommandQueue commandQueue;
  vector<cl::Device> devices;
  vector<cl::Platform> platforms;
/**************************************************************************************************/
  cl_uint currentTimeStep; //TODO: convert to largerst positive int
  cl_uint currentTimeSlot;
  unsigned int srandSeed, srandCounter;
  std::string startTimeStamp;
  kernelStatistics kernelStats;
/**************************************************************************************************/



/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Constructor.
  */
  Neurosim() :
    /* *** */
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    dataMakeEventPtrsStructSize(0),
    dataMakeEventPtrsStructSizeBytes(0),
    dataMakeEventPtrsStruct(NULL),
    dataMakeEventPtrsStructBuffer(),
#endif
    /* *** */
#if MAKE_EVENT_PTRS_ENABLE
    kernelMakeEventPtrs(),
    blockSizeX_kernelMakeEventPtrs(MAKE_EVENT_PTRS_WG_SIZE_WI),
    blockSizeY_kernelMakeEventPtrs(1),
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    kernelGlueEventPtrs(),
    blockSizeX_kernelGlueEventPtrs(GLUE_EVENT_PTRS_WG_SIZE_WI),
    blockSizeY_kernelGlueEventPtrs(1),
#endif
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    dataMakeEventPtrsDebugHostSize(0),
    dataMakeEventPtrsDebugHostSizeBytes(0),
    dataMakeEventPtrsDebugHost(NULL),
    dataMakeEventPtrsDebugHostBuffer(),
    dataMakeEventPtrsDebugDeviceSize(0),
    dataMakeEventPtrsDebugDeviceSizeBytes(0),
    dataMakeEventPtrsDebugDevice(NULL),
    dataMakeEventPtrsDebugDeviceBuffer(),
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    dataMakeEventPtrsErrorSize(0),
    dataMakeEventPtrsErrorSizeBytes(0),
    dataMakeEventPtrsError(NULL),
    dataMakeEventPtrsErrorBuffer(),
#endif
#endif
    /* *** */
#if UPDATE_NEURONS_ENABLE_V00
  /*Model variables and parameters*/
    nrn_ps(NULL),
    ne(NULL),
    te_ps(NULL),
    co(NULL),
    tau_ampa_ps(0),
    tau_gaba_ps(0),
    co_g_ampa_ps(0),
    co_g_gaba_ps(0),
    E_ps(0),
    a_ps(0),
    dt_ps(0),
    tol_ps(0),
    steps_ps(0),
    modelParametersSize(0),
    modelParametersSizeBytes(0),
    modelParameters(NULL),
    modelParametersBuffer(),
    modelVariablesSize(0),
    modelVariablesSizeBytes(0),
    modelVariables(NULL),
    modelVariablesBuffer(),
    constantCoefficientsSize(0),
    constantCoefficientsSizeBytes(0),
    constantCoefficients(NULL),
    constantCoefficientsBuffer(),
    kernelUpdateNeuronsV00(),
    kernelUpdateSpikedNeuronsV00(),
    blockSizeX_kernelUpdateNeuronsV00(UPDATE_NEURONS_WG_SIZE_WI_V00),
    blockSizeY_kernelUpdateNeuronsV00(1),
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    psToleranceSize(0),
    psToleranceSizeBytes(0),
    psTolerance(NULL),
    psToleranceBuffer(),
#endif
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    dataUpdateNeuronsDebugHostSize(0),
    dataUpdateNeuronsDebugHostSizeBytes(0),
    dataUpdateNeuronsDebugHost(NULL),
    dataUpdateNeuronsDebugHostBuffer(),
    dataUpdateNeuronsDebugDeviceSize(0),
    dataUpdateNeuronsDebugDeviceSizeBytes(0),
    dataUpdateNeuronsDebugDevice(NULL),
    dataUpdateNeuronsDebugDeviceBuffer(),
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    dataUpdateNeuronsErrorSize(0),
    dataUpdateNeuronsErrorSizeBytes(0),
    dataUpdateNeuronsError(NULL),
    dataUpdateNeuronsErrorBuffer(),
#endif
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
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
    performanceCounter(0),
    performanceFrequency(0),
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
    spikeEvents(NULL),
    connectome(NULL),
    synapticEvents(NULL),
    /* *** */
    context(),
    commandQueue(),
    devices(),
    platforms(),
    /* *** */
    currentTimeStep(0),
    currentTimeSlot(0),
    srandSeed(1),
    srandCounter(0),
    startTimeStamp(std::string("")),
    kernelStats()
/**************************************************************************************************/
  {
    time_t rawtime; 
    time ( &rawtime ); 
    this->startTimeStamp = std::string(ctime (&rawtime));
    
    SET_RANDOM_SEED_DIRECT(this->srandSeed);
    
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
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
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
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
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
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

#if UPDATE_NEURONS_ENABLE_V00
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
#endif

#if UPDATE_NEURONS_ENABLE_V00
  int 
  injectSortedEvents
  (
    bool          verify,
    bool          resetEventBuffer,
    cl_uint       totalNeurons,
    cl_uint       eventQueueSize,
    cl_uint       pointersPitch,
    size_t        sortedEventsSize,
    cl_uint       *sortedEvents,
    cl_uint       *pointersToEvents,
    neuron_iz_ps  *nrn
  );
#endif

#if UPDATE_NEURONS_ENABLE_V00
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
#endif

#if UPDATE_NEURONS_ENABLE_V00
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
#endif

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
    size_t,
    unsigned int*,
    unsigned int*,
    DATA_TYPE*,
    Data_SpikeEvents&,
    Data_Connectome&,
    cl::CommandQueue&
  );
#endif

#if UPDATE_NEURONS_ENABLE_V00
  int 
  verifyEvents
  (
    bool          ignoreWarnings,
    bool          correctWeightPositionMismatch,
    unsigned int  totalNeurons,
    unsigned int  structElementSize,
    size_t        sortedEventsSize,
    unsigned int  *pointerStruct,
    unsigned int  *sortedEvents,
    neuron_iz_ps  *nrn
  );
#endif

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
