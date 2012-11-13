
#ifndef SPIKE_EVENTS_H_
#define SPIKE_EVENTS_H_



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

#define SPIKE_EVENTS_VALID_SPIKES         0x1
#define SPIKE_EVENTS_VALID_ALL            (SPIKE_EVENTS_VALID_SPIKES)

/**************************************************************************************************/



/***************************************************************************************************
  Forward declarations for cyclic dependency
***************************************************************************************************/



/**************************************************************************************************/



/**
  @class Data_SpikeEvents

  Models spike events.

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/04/17
 */
class Data_SpikeEvents 
{
/**************************************************************************************************/
  public:  /*public variables*/
/**************************************************************************************************/



  cl::Buffer dataSpikePacketsBuffer;
  cl::Buffer dataSpikePacketCountsBuffer;
  

  
/**************************************************************************************************/
  private:  /*private variables*/
/**************************************************************************************************/



  cl_uint neuronCount;
  cl_uint spikePacketSize;
  cl_uint spikePacketSizeWords;
  cl_uint spikePackets;
  cl_uint simulationStepSize;
  cl_uint spikeDatumSize;

  size_t dataSpikePacketsSize;
  size_t dataSpikePacketsSizeBytes;
  cl_uint *dataSpikePackets;
  
  size_t dataSpikePacketCountsSize;
  size_t dataSpikePacketCountsSizeBytes;
  cl_uint *dataSpikePacketCounts;

  cl_uint *dataPastSpikePackets;
  cl_uint *dataPastSpikePacketCounts;

  std::stringstream *dataToSimulationLogFile;
  std::stringstream *dataToReportLogFile;
  
  cl_uint dataValid;
  
  unsigned int srandSeed;
  unsigned int srandCounter;

  
  
/**************************************************************************************************/
  public: /*public methods*/
/**************************************************************************************************/


  
/**************************************************************************************************/
  /**
    Constructor.
  */
  Data_SpikeEvents
  (
    cl::Context                         &context,
    cl::Device                          &device,
    cl::CommandQueue                    &queue,
    cl_bool                             block,
    struct kernelStatistics             &kernelStats,
    std::stringstream                   *dataToSimulationLogFile,
    std::stringstream                   *dataToReportLogFile,
    cl_uint                             simulationStepSize,
    cl_uint                             spikeDatumSize,
    cl_uint                             spikePacketSizeWords,
    cl_uint                             spikePacketSize,
    cl_uint                             spikePackets,
    cl_uint                             neuronCount
  ) :
    dataSpikePacketsBuffer(),
    dataSpikePacketCountsBuffer(),
    neuronCount(neuronCount),
    spikePacketSize(spikePacketSize),
    spikePacketSizeWords(spikePacketSizeWords),
    spikePackets(spikePackets),
    simulationStepSize(simulationStepSize),
    spikeDatumSize(spikeDatumSize),
    dataSpikePacketsSize(0),
    dataSpikePacketsSizeBytes(0),
    dataSpikePackets(NULL),
    dataSpikePacketCountsSize(0),
    dataSpikePacketCountsSizeBytes(0),
    dataSpikePacketCounts(NULL),
    dataPastSpikePackets(NULL),
    dataPastSpikePacketCounts(NULL),
    dataToSimulationLogFile(dataToSimulationLogFile),
    dataToReportLogFile(dataToReportLogFile),
    dataValid(0),
    srandSeed(1),
    srandCounter(0)
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
  ~Data_SpikeEvents()
  {
    if(this->dataSpikePacketCounts)
    {
      free(this->dataSpikePacketCounts);
      this->dataSpikePacketCounts = NULL;
    }
    if(this->dataSpikePackets)
    {
      free(this->dataSpikePackets);
      this->dataSpikePackets = NULL;
    }
    if(this->dataPastSpikePacketCounts)
    {
      free(this->dataPastSpikePacketCounts);
      this->dataPastSpikePacketCounts = NULL;
    }
    if(this->dataPastSpikePackets)
    {
      free(this->dataPastSpikePackets);
      this->dataPastSpikePackets = NULL;
    }
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
  };
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Clears all spike events.
  */
  void
  clearEvents
  (
    cl::CommandQueue&,
    cl_bool
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Sets all spike events according to parameters
  */
  void 
  setEvents
  (
    cl::CommandQueue&,
    cl_bool,
    double,
    double,
    double
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns current spike packet count.
  */
  cl_uint
  getSpikePacketCount
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns spike count in a spike packet
  */
  cl_uint
  getSpikeCount
  (
    cl::CommandQueue&,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns the past spike count.
  */
  cl_uint
  getPastSpikeCount
  (
    cl::CommandQueue&,
    cl_uint
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns spike datum (spiked neuron and time) by reference for given spike and packet number.
  */
  void
  getSpike
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint&,
    CL_DATA_TYPE&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns by reference the past spike datum for given spike and packet number.
  */
  void
  getPastSpike
  (
    cl::CommandQueue&,
    cl_uint,
    cl_uint,
    cl_uint&,
    CL_DATA_TYPE&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Invalidates current spike events (they become past spikes).
  */
  void 
  invalidateEvents
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Loads data from device to host (if it is not valid on the host).
  */
  void 
  refresh
  (
    cl::CommandQueue&
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



#endif  /*SPIKE_EVENTS_H_*/
