
#ifndef SYNAPTIC_EVENTS_H_
#define SYNAPTIC_EVENTS_H_



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Common.hpp"

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
  @class SynapticEvents

  Models synaptic events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/05/17
 */
class SynapticEvents 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/
  
  cl_uint timeSlots;
  cl_uint histogramBinCount;
  cl_uint histogramBinSize;
  
/*TODO: eventually they have to move to private*/
  static const int RECENT   = 0;
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

  bool resetObject;
  cl_uint dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;
  
  cl_uint eventBufferCount;
  cl_uint eventBufferSize;
  
  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;

#if SYNAPTIC_EVENTS_ENABLE_PAST_EVENTS
  cl_uint* dataPastUnsortedEventCounts;
  cl_uint* dataPastUnsortedEventTargets;
  cl_uint* dataPastUnsortedEventDelays;
  cl_uint* dataPastUnsortedEventWeights;
  cl_uint* dataPastHistogram;
#endif


  
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  SynapticEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Destructor.
  */
  ~SynapticEvents
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
    cl_uint, 
    cl_uint,
    cl_uint, 
    cl_uint, 
    cl_uint,
    struct kernelStatistics*,
    std::stringstream*,
    std::stringstream*
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
    Initializes values for all unsorted synaptic events according to passed parameters
  */
  void 
  setUnsortedEvents
  (
    cl::CommandQueue&,
    cl_bool,
    cl_bool,
    cl_uint,
    cl_uint,
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
  getRecentHistogramItem
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
    Assures that the date is up to date.
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
    Initializes event buffers
  */
  void 
  initializeEventBuffers
  (
    cl_uint,
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
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Saves past events and loads current events from the device if they are invalid on the host.
  */
  void
  getEventBuffers
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
    Stores event data on device.
  */
  void
  storeBuffers
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
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
