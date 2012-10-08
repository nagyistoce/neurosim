
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

#define OPERATOR_SORT_VALID_HISTOGRAM_TIK         0x1
#define OPERATOR_SORT_VALID_ALL                   (OPERATOR_SORT_VALID_HISTOGRAM_TIK)
#define OPERATOR_SORT_ENABLE_PAST_EVENTS          1

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/



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
  
  
  
  /*TODO: eventually all of them have to move to private*/
  size_t dataHistogramGroupEventsTikSize;
  size_t dataHistogramGroupEventsTikSizeBytes;
  cl_uint* dataHistogramGroupEventsTik;
  cl::Buffer dataHistogramGroupEventsTikBuffer;
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
  cl_uint* dataPastHistogramGroupEventsTik;
#endif
  
  
  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



  static const int RECENT = 0;
  static const int PREVIOUS = 1;
  
  cl_uint dataValid;
  cl_uint histogramBacketCount;
  cl_uint histogramBinSize;
  cl_uint histogramBinCount;
  
  unsigned int srandSeed;
  unsigned int srandCounter;

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
    cl::Context                         &context,
    cl::Device                          &device,
    cl::CommandQueue                    &queue,
    cl_bool                             block,
    struct kernelStatistics             &kernelStats,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile,
    cl_uint                             histogramBacketCount,
    cl_uint                             histogramBinSize,
    cl_uint                             histogramBinCount
  ) :
    dataHistogramGroupEventsTikSize(0),
    dataHistogramGroupEventsTikSizeBytes(0),
    dataHistogramGroupEventsTik(NULL),
    dataHistogramGroupEventsTikBuffer(),
    /* *** */
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
    dataPastHistogramGroupEventsTik(NULL),
#endif
    /* *** */
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile),
    dataValid(0),
    srandSeed(1),
    srandCounter(0),
    histogramBacketCount(histogramBacketCount),
    histogramBinSize(histogramBinSize),
    histogramBinCount(histogramBinCount)
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
  ~Operator_Sort()
  {
    if(this->dataHistogramGroupEventsTik)
    {
      free(this->dataHistogramGroupEventsTik);
      this->dataHistogramGroupEventsTik = NULL;
    }
    /* *** */
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
    if(this->dataPastHistogramGroupEventsTik)
    {
      free(this->dataPastHistogramGroupEventsTik);
      this->dataPastHistogramGroupEventsTik = NULL;
    }
#endif
    /* *** */
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Clears data.
  */
  void
  clearData
  (
    cl::CommandQueue&,
    cl_bool,
    cl_uint
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
    Invalidates current histogram.
  */
  void 
  invalidateHistogram
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Refreshes host data with data from device.
  */
  void 
  refreshHistogram
  (
    cl::CommandQueue&,
    cl_bool
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

};

#endif  /*ENABLE_OPERATOR_SORT*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif  /*OPERATOR_SORT_H_*/
