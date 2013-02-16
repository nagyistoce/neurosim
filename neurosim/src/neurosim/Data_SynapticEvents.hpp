
#ifndef SYNAPTIC_EVENTS_H_
#define SYNAPTIC_EVENTS_H_



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
  Class Preprocessor Definitions
***************************************************************************************************/

#define SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS                  1
#define SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS         0x1
#define SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS               0x2
#define SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM      0x4
#define SYNAPTIC_EVENTS_VALID_SORTED_EVENTS                 0x8
#define SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM       0x10
#define SYNAPTIC_EVENTS_VALID_ALL   (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS |\
                                     SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS |\
                                     SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_HISTOGRAM |\
                                     SYNAPTIC_EVENTS_VALID_SORTED_EVENTS |\
                                     SYNAPTIC_EVENTS_VALID_SORTED_EVENTS_HISTOGRAM)

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/

class Operator_Sort;

/**************************************************************************************************/



/**
  @class Data_SynapticEvents

  Models synaptic events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/05/17
 */
class Data_SynapticEvents 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/
  
  
  /*TODO: eventually all of them have to move to private*/
  cl_uint timeSlots;
  cl_uint histogramBinCount;
  cl_uint histogramBinSize;
  
  cl_uint sortedEventsHistogramBacketCount;
  cl_uint sortedEventsHistogramBinSize;
  cl_uint sortedEventsHistogramBinCount;
  
  static const int RECENT = 0;
  static const int PREVIOUS = 1;
  
  cl::Buffer dataUnsortedEventTargetsBuffer;
  cl::Buffer dataUnsortedEventDelaysBuffer;
  cl::Buffer dataUnsortedEventWeightsBuffer;
  cl::Buffer dataUnsortedEventCountsBuffer;
  cl::Buffer dataUnsortedEventsHistogramBuffer;
  cl::Buffer dataSortedEventsReadBuffer;
  cl::Buffer dataSortedEventsWriteBuffer;
  cl::Buffer dataSortedEventsHistogramWriteBuffer;
  cl::Buffer dataSortedEventsHistogramReadBuffer;

  size_t dataUnsortedEventTargetsSize;
  size_t dataUnsortedEventTargetsSizeBytes;
  cl_uint* dataUnsortedEventTargets;
  
  size_t dataUnsortedEventDelaysSize;
  size_t dataUnsortedEventDelaysSizeBytes;
  cl_uint* dataUnsortedEventDelays;
  
  size_t dataUnsortedEventWeightsSize;
  size_t dataUnsortedEventWeightsSizeBytes;
  cl_uint* dataUnsortedEventWeights;
  
  size_t dataUnsortedEventCountsSize;
  size_t dataUnsortedEventCountsSizeBytes;
  cl_uint* dataUnsortedEventCounts;

  size_t dataUnsortedEventsHistogramSize;
  size_t dataUnsortedEventsHistogramSizeBytes;
  cl_uint* dataUnsortedEventsHistogram;
  
  size_t dataSortedEventsSize;
  size_t dataSortedEventsSizeBytes;
  cl_uint* dataSortedEvents;

  size_t dataSortedEventsHistogramSize;
  size_t dataSortedEventsHistogramSizeBytes;
  cl_uint* dataSortedEventsHistogram;  
  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



  cl_uint dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;
  
  cl_uint eventBufferCount;
  cl_uint eventBufferSize;
  cl_uint sortedEventBufferCount;
  cl_uint sortedEventBufferSize;
  
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  cl_uint* dataPastUnsortedEventCounts;
  cl_uint* dataPastUnsortedEventTargets;
  cl_uint* dataPastUnsortedEventDelays;
  cl_uint* dataPastUnsortedEventWeights;
  cl_uint* dataPastUnsortedEventsHistogram;
  cl_uint* dataPastSortedEvents;
  cl_uint* dataPastSortedEventsHistogram;
#endif

  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;


 
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  Data_SynapticEvents
  (
    cl::Context                         &context,
    cl::Device                          &device,
    cl::CommandQueue                    &queue,
    cl_bool                             block,
    struct kernelStatistics             &kernelStats,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile,
    cl_uint                             timeSlots,
    cl_uint                             eventBufferCount,
    cl_uint                             eventBufferSize,
    cl_uint                             histogramBinCount,
    cl_uint                             histogramBinSize,
    cl_uint                             sortedEventsHistogramBacketCount,
    cl_uint                             sortedEventsHistogramBinSize,
    cl_uint                             sortedEventsHistogramBinCount,
    cl_uint                             sortedEventBufferCount,
    cl_uint                             sortedEventBufferSize
  ) :
    timeSlots(timeSlots),
    histogramBinCount(histogramBinCount),
    histogramBinSize(histogramBinSize),
    sortedEventsHistogramBacketCount(sortedEventsHistogramBacketCount),
    sortedEventsHistogramBinSize(sortedEventsHistogramBinSize),
    sortedEventsHistogramBinCount(sortedEventsHistogramBinCount),
    dataUnsortedEventTargetsBuffer(),
    dataUnsortedEventDelaysBuffer(),
    dataUnsortedEventWeightsBuffer(),
    dataUnsortedEventCountsBuffer(),
    dataUnsortedEventsHistogramBuffer(),
    dataSortedEventsReadBuffer(),
    dataSortedEventsWriteBuffer(),
    dataSortedEventsHistogramWriteBuffer(),
    dataSortedEventsHistogramReadBuffer(),
    dataUnsortedEventTargetsSize(0),
    dataUnsortedEventTargetsSizeBytes(0),
    dataUnsortedEventTargets(NULL),
    dataUnsortedEventDelaysSize(0),
    dataUnsortedEventDelaysSizeBytes(0),
    dataUnsortedEventDelays(NULL),
    dataUnsortedEventWeightsSize(0),
    dataUnsortedEventWeightsSizeBytes(0),
    dataUnsortedEventWeights(NULL),
    dataUnsortedEventCountsSize(0),
    dataUnsortedEventCountsSizeBytes(0),
    dataUnsortedEventCounts(NULL),
    dataUnsortedEventsHistogramSize(0),
    dataUnsortedEventsHistogramSizeBytes(0),
    dataUnsortedEventsHistogram(NULL),
    dataSortedEventsSize(0),
    dataSortedEventsSizeBytes(0),
    dataSortedEvents(NULL),    
    dataSortedEventsHistogramSize(0),
    dataSortedEventsHistogramSizeBytes(0),
    dataSortedEventsHistogram(NULL),    
    dataValid(0),
    srandSeed(1),
    srandCounter(0),
    eventBufferCount(eventBufferCount),
    eventBufferSize(eventBufferSize),
    sortedEventBufferCount(sortedEventBufferCount),
    sortedEventBufferSize(sortedEventBufferSize),
    /* *** */
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    dataPastUnsortedEventCounts(NULL),
    dataPastUnsortedEventTargets(NULL),
    dataPastUnsortedEventDelays(NULL),
    dataPastUnsortedEventWeights(NULL),
    dataPastUnsortedEventsHistogram(NULL),
    dataPastSortedEvents(NULL),
    dataPastSortedEventsHistogram(NULL),
#endif
    /* *** */
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile)
/**************************************************************************************************/
  {
    /* Must be called only from the constructor and only once*/
    this->initialize
    (
      context,
      device,
      queue,
      block,
      kernelStats
    );
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~Data_SynapticEvents()
  {
    if(this->dataUnsortedEventCounts)
    {
      free(this->dataUnsortedEventCounts);
      this->dataUnsortedEventCounts = NULL;
    }
    if(this->dataUnsortedEventWeights)
    {
      free(this->dataUnsortedEventWeights);
      this->dataUnsortedEventWeights = NULL;
    }
    if(this->dataUnsortedEventDelays)
    {
      free(this->dataUnsortedEventDelays);
      this->dataUnsortedEventDelays = NULL;
    }
    if(this->dataUnsortedEventTargets)
    {
      free(this->dataUnsortedEventTargets);
      this->dataUnsortedEventTargets = NULL;
    }
    if(this->dataUnsortedEventsHistogram)
    {
      free(this->dataUnsortedEventsHistogram);
      this->dataUnsortedEventsHistogram = NULL;
    }
    if(this->dataSortedEvents)
    {
      free(this->dataSortedEvents);
      this->dataSortedEvents = NULL;
    }    
    if(this->dataSortedEventsHistogram)
    {
      free(this->dataSortedEventsHistogram);
      this->dataSortedEventsHistogram = NULL;
    }    
    /* *** */
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    if(this->dataPastUnsortedEventCounts)
    {
      free(this->dataPastUnsortedEventCounts);
      this->dataPastUnsortedEventCounts = NULL;
    }
    if(this->dataPastUnsortedEventWeights)
    {
      free(this->dataPastUnsortedEventWeights);
      this->dataPastUnsortedEventWeights = NULL;
    }
    if(this->dataPastUnsortedEventDelays)
    {
      free(this->dataPastUnsortedEventDelays);
      this->dataPastUnsortedEventDelays = NULL;
    }
    if(this->dataPastUnsortedEventTargets)
    {
      free(this->dataPastUnsortedEventTargets);
      this->dataPastUnsortedEventTargets = NULL;
    }
    if(this->dataPastUnsortedEventsHistogram)
    {
      free(this->dataPastUnsortedEventsHistogram);
      this->dataPastUnsortedEventsHistogram = NULL;
    }
    if(this->dataPastSortedEvents)
    {
      free(this->dataPastSortedEvents);
      this->dataPastSortedEvents = NULL;
    }    
    if(this->dataPastSortedEventsHistogram)
    {
      free(this->dataPastSortedEventsHistogram);
      this->dataPastSortedEventsHistogram = NULL;
    }    
#endif
    /* *** */
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Clears all sorted synaptic events.
  */
  void
  clearSortedEvents
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Clears all unsorted synaptic events.
  */
  void
  clearUnsortedEvents
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes values for unsorted synaptic events in external buffers
  */
  void 
  initializeGrouppedEvents
  (
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
    cl_float,
    double,
    double,
    cl::CommandQueue&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes values for all unsorted synaptic events according to passed parameters
  */
  void 
  setUnsortedEvents
  (
#if CLASS_VALIDATION_ENABLE
    cl_uint,
#endif
    cl::CommandQueue&,
    cl_bool,
    cl_bool,
    int,
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    double,
    double
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Sets sorted histogram to speficied value.
  */
  void 
  setSortedHistogram
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint,
    int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns event count at specified time slot
  */
  cl_uint
  getEventCount
  (
    cl::CommandQueue&,
    cl_uint,
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns unsorted event count at specified buffer ID and time slot
  */
  cl_uint
  getEventCount
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns event data by reference from sorted events
  */
  cl_uint
  getEvent
  (
    cl::CommandQueue&,
    cl_uint,
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns event data by reference from unsorted events
  */
  void
  getEvent
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    const int,
    cl_uint&,
    cl_uint&,
    cl_uint&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns event data by reference from unsorted events
  */
  void
  getEvent
  (
    cl::CommandQueue&,
    cl_uint,
    const int,
    cl_uint&,
    cl_uint&,
    cl_uint&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram by reference
  */
  void
  getHistogram
  (
    cl::CommandQueue&,
    cl_uint**,
#if CLASS_VALIDATION_ENABLE
    size_t,
#endif
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  static cl_uint
  getPreviousUnsortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  cl_uint
  getPreviousUnsortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  static cl_uint
  getCurrentUnsortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  cl_uint
  getCurrentUnsortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  static cl_uint
  getPreviousSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  cl_uint
  getPreviousSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  cl_uint
  getPreviousSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  static cl_uint
  getCurrentSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    void*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  cl_uint
  getCurrentSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item
  */
  cl_uint
  getCurrentSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Instantiates and returns a pointer to a copy of past sorted events histogram
  */
  void
  getPreviousSortedHistogram
  (
#if CLASS_VALIDATION_ENABLE
    size_t,
#endif
    cl::CommandQueue&,
    cl_uint**
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Instantiates and returns a pointer to a copy of current sorted events histogram
  */
  void
  getCurrentSortedHistogram
  (
#if CLASS_VALIDATION_ENABLE
    size_t,
#endif
    cl::CommandQueue&,
    cl_uint**
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Refreshes host data with data from device.
  */
  void 
  refreshSortedEventsHistogram
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current sorted synaptic events.
  */
  void 
  invalidateSortedEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current unsorted synaptic events.
  */
  void 
  invalidateUnsortedEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current histogram for sorted synaptic events.
  */
  void 
  invalidateSortedHistogram
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current histogram for unsorted synaptic events.
  */
  void 
  invalidateUnsortedHistogram
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Refreshes host data with data from device for unsorted synaptic events histogram.
  */
  void 
  refreshUnsortedHistogram
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Refreshes host data with data from device for sorted synaptic events histogram.
  */
  void 
  refreshSortedHistogram
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Refreshes host data with data from device for sorted synaptic events.
  */
  void 
  refreshSortedEvents
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns number of event buffers in this object
  */
  cl_uint 
  getEventBufferCount
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns number of event time slots in this object
  */
  cl_uint 
  getEventTimeSlots
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns number of event buffer size in this object
  */
  cl_uint 
  getEventBufferSize
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns number of histogram bins in this object
  */
  cl_uint 
  getEventHistogramBinCount
  ();
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
    cl::Context&,
    cl::Device&,
    cl::CommandQueue&,
    cl_bool,
    struct kernelStatistics&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes event buffers
  */
  void 
  initializeEventBuffers
  (
    int,
    cl_uint,
    cl_uint,
    cl_uint,
    double,
    double,
    double,
    double
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes event buffer
  */
  void 
  initializeEventBuffer
  (
    cl_uint,
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
    cl_uint*,
    cl_uint*,
    cl_uint*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Group sorted events
  */
  void 
  groupEvents
  (
    cl::CommandQueue&,
    cl_bool,
    bool,
    size_t,
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
    cl_uint*,
    cl_uint*,
    cl_uint*,
    cl_uint*,
    cl_uint*,
    cl_uint*
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Initializes histogram
  */
  void 
  initializeHistogram
  (
#if CLASS_VALIDATION_ENABLE
    cl_uint
#endif
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item from histogram of unsorted events
  */
  cl_uint
  getUnsortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item from histogram of sorted events
  */
  cl_uint
  getSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint,
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns histogram item from histogram of sorted events
  */
  cl_uint
  getSortedHistogramItem
  (
    cl::CommandQueue&,
    cl_uint,
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Instantiates and returns a pointer to a copy of sorted events histogram
  */
  void
  getSortedHistogram
  (
    cl::CommandQueue&,
    cl_uint**,
#if CLASS_VALIDATION_ENABLE
    size_t,
#endif
    const int
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads data from the device if it is invalid on the host.
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
    Stores data on device.
  */
  void
  storeData
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
  );
/**************************************************************************************************/
};



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif  /*SYNAPTIC_EVENTS_H_*/
