
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

#define SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS 0x1
#define SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS       0x2
#define SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM        0x4
#define SYNAPTIC_EVENTS_VALID_ALL                   (SYNAPTIC_EVENTS_VALID_UNSORTED_EVENT_COUNTS |\
                                                     SYNAPTIC_EVENTS_VALID_UNSORTED_EVENTS |\
                                                     SYNAPTIC_EVENTS_VALID_HISTOGRAM_ITEM)
#define SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS          1

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/



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
  
  static const int RECENT = 0;
  static const int PREVIOUS = 1;
  
  cl::Buffer dataUnsortedEventTargetsBuffer;
  cl::Buffer dataUnsortedEventDelaysBuffer;
  cl::Buffer dataUnsortedEventWeightsBuffer;
  cl::Buffer dataUnsortedEventCountsBuffer;
  cl::Buffer dataHistogramBuffer;

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
  
  size_t dataHistogramSize;
  size_t dataHistogramSizeBytes;
  cl_uint* dataHistogram;
  

  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



  cl_uint dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;
  
  cl_uint eventBufferCount;
  cl_uint eventBufferSize;
  
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  cl_uint* dataPastUnsortedEventCounts;
  cl_uint* dataPastUnsortedEventTargets;
  cl_uint* dataPastUnsortedEventDelays;
  cl_uint* dataPastUnsortedEventWeights;
  cl_uint* dataPastHistogram;
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
    cl_uint                             histogramBinSize
  ) :
    timeSlots(timeSlots),
    histogramBinCount(histogramBinCount),
    histogramBinSize(histogramBinSize),
    dataUnsortedEventTargetsBuffer(),
    dataUnsortedEventDelaysBuffer(),
    dataUnsortedEventWeightsBuffer(),
    dataUnsortedEventCountsBuffer(),
    dataHistogramBuffer(),
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
    dataHistogramSize(0),
    dataHistogramSizeBytes(0),
    dataHistogram(NULL),
    dataValid(0),
    srandSeed(1),
    srandCounter(0),
    eventBufferCount(eventBufferCount),
    eventBufferSize(eventBufferSize),
    /* *** */
#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
    dataPastUnsortedEventCounts(NULL),
    dataPastUnsortedEventTargets(NULL),
    dataPastUnsortedEventDelays(NULL),
    dataPastUnsortedEventWeights(NULL),
    dataPastHistogram(NULL),
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
    if(this->dataHistogram)
    {
      free(this->dataHistogram);
      this->dataHistogram = NULL;
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
    if(this->dataPastHistogram)
    {
      free(this->dataPastHistogram);
      this->dataPastHistogram = NULL;
    }
#endif
    /* *** */
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
  };
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
    Returns event data by reference
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
    Returns histogram item
  */
  cl_uint
  getHistogramItem
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
    Returns histogram item
  */
  static cl_uint
  getPreviousHistogramItem
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
  static cl_uint
  getCurrentHistogramItem
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
    Invalidates current synaptic events.
  */
  void 
  invalidateEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current histogram.
  */
  void 
  invalidateHistogram
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Refreshes host data with data from device
  */
  void 
  refresh
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
