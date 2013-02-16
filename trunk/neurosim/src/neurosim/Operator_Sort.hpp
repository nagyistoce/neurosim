
#ifndef OPERATOR_SORT_H_
#define OPERATOR_SORT_H_



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



#if ENABLE_OPERATOR_SORT



/***************************************************************************************************
  Class Preprocessor Definitions
***************************************************************************************************/



/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/

class Data_SynapticEvents;

#if ENABLE_OPERATOR_SCAN
class Operator_Scan;
#endif
/* *** */
#if ENABLE_OPERATOR_GROUP
class Operator_Group;
#endif

/**************************************************************************************************/



/**
  @class Operator_Sort

  Models synaptic events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/09/15
 */
class Operator_Sort 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/
  

  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/

  cl_uint sortByTimeIterationCount;
  cl_uint sortByNeuronIterationCount;
  
#if ENABLE_UNIT_TEST_GROUP_EVENTS
  double groupTestMode;
#endif

#if ENABLE_OPERATOR_SCAN
  Operator_Scan *operatorScan;
#endif

#if ENABLE_OPERATOR_GROUP
  Operator_Group *operatorGroup;
#endif

#if SORT_VERIFY_ENABLE
  cl_uint dataUnsortedEventsSnapShotSize;
  cl_uint dataUnsortedEventsSnapShotSizeBytes;
  cl_uint* dataUnsortedEventsSnapShot;
#endif

  Data_SynapticEvents *synapticEvents;
  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;


 
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  Operator_Sort
  (
#if ENABLE_UNIT_TEST_GROUP_EVENTS
    double                              groupTestMode,
#endif
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || ENABLE_UNIT_TEST_SCAN_V00 || \
  ENABLE_UNIT_TEST_SCAN_V01) || (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
    cl_bool                             block,
    cl::CommandQueue                    &queue,
#endif
    cl_uint                             sortByTimeIterationCount,
    cl_uint                             sortByNeuronIterationCount,
    cl::Context                         &context,
    cl::Device                          &device,
    struct kernelStatistics             &kernelStats,
    Data_SynapticEvents                 *synapticEvents,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile
  ) :
    sortByTimeIterationCount(sortByTimeIterationCount),
    sortByNeuronIterationCount(sortByNeuronIterationCount),
    /* *** */
#if ENABLE_UNIT_TEST_GROUP_EVENTS
    groupTestMode(groupTestMode),
#endif
    /* *** */
#if ENABLE_OPERATOR_SCAN
    operatorScan(NULL),
#endif
    /* *** */
#if ENABLE_OPERATOR_GROUP
    operatorGroup(NULL),
#endif
    /* *** */
#if SORT_VERIFY_ENABLE
    dataUnsortedEventsSnapShotSize(0),
    dataUnsortedEventsSnapShotSizeBytes(0),
    dataUnsortedEventsSnapShot(NULL),
#endif
    /* *** */
    synapticEvents(synapticEvents),
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile)
/**************************************************************************************************/
  {
    /* Must be called only from the constructor and only once*/
    this->initialize
    (
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || ENABLE_UNIT_TEST_SCAN_V00 || \
  ENABLE_UNIT_TEST_SCAN_V01) || (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
      block,
      queue,
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
  ~Operator_Sort();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs sort of event data.
  */
#if (ENABLE_OPERATOR_SCAN || ENABLE_OPERATOR_GROUP)
  void
  sortEvents
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    cl_uint,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE || \
(ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Takes a snapshot of unsorted events and stores it in memory.
  */
#if SORT_VERIFY_ENABLE
  void 
  captureUnsortedEvents
  (
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    cl::CommandQueue&,
    Data_SynapticEvents&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Verifies sort operation.
  */
#if SORT_VERIFY_ENABLE
  void 
  verifySortedEvents
  (
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    cl_uint*,
    cl_uint*
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
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || ENABLE_UNIT_TEST_SCAN_V00 || \
  ENABLE_UNIT_TEST_SCAN_V01) || (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
    cl_bool,
    cl::CommandQueue&,
#endif
    cl::Context&,
    cl::Device&,
    struct kernelStatistics&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs first pass of sort of time data.
  */
#if (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
  void
  sortByTimeFirstPass
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    cl_uint,
#endif
    cl_uint,
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs sort of data based on time key.
  */
#if ((ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00) || \
    (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V02))
  void
  sortByTime
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    cl_uint,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE || \
(ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Performs sort of data based on neuron ID key.
  */
#if ((ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V03))
  void
  sortByNeuron
  (
#if KERNEL_LEVEL_PROFILING
    struct kernelStatistics&,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    cl_uint,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl::Event&
  );
#endif
/**************************************************************************************************/
};

#endif  /*ENABLE_OPERATOR_SORT*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif  /*OPERATOR_SORT_H_*/
