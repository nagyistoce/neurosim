
/* ===============================================================================================



  =============================================================================================== */


  
  
/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Operator_Sort.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



#if ENABLE_OPERATOR_SORT



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



Operator_Sort::~Operator_Sort
()
/**************************************************************************************************/
{
  /* *** */
#if ENABLE_OPERATOR_SCAN
  if(this->operatorScan)
    delete(this->operatorScan);
#endif
/* *** */
#if ENABLE_OPERATOR_GROUP
  if(this->operatorGroup)
    delete(this->operatorGroup);
#endif
  /* *** */
#if SORT_VERIFY_ENABLE
  if(this->dataUnsortedEventsSnapShot)
    free(this->dataUnsortedEventsSnapShot);
#endif
  /* *** */
  this->synapticEvents = NULL;
  this->dataToSimulationLogFile = NULL;
  this->dataToReportLogFile = NULL;
}
/**************************************************************************************************/



#if (ENABLE_OPERATOR_SCAN || ENABLE_OPERATOR_GROUP)
void
Operator_Sort::sortEvents
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
  cl_uint                   currentTimeSlot,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE || \
(ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt
)
/**************************************************************************************************/
{
#if ((ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00) || \
    (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V02))
  this->sortByTime
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    currentTimeSlot,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE || \
(ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
    currentTimeStep,
#endif
    queue,
    ndrEvt
  );
#endif

#if ((ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V03))
  this->sortByNeuron
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    currentTimeSlot,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
    currentTimeStep,
#endif
    queue,
    ndrEvt
  );
#endif
}
/**************************************************************************************************/
#endif



#if SORT_VERIFY_ENABLE
void 
Operator_Sort::captureUnsortedEvents
(
  cl_uint                 currentTimeSlot,
  cl_uint                 eventBufferCount,
  cl_uint                 eventDataSize,
  double                  maxDelay,
  double                  minDelay,
  cl::CommandQueue        &queue,
  Data_SynapticEvents     &synapticEvents
)
/**************************************************************************************************/
{
  cl_uint ptr_store = 0;
  cl_uint event_total = 0;
  bool result = 1;
  std::stringstream ss;
  
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < eventBufferCount; b++)
  {
    event_total += synapticEvents.getEventCount(queue, b, currentTimeSlot, 
      Data_SynapticEvents::RECENT);
  }
  
  if(this->dataUnsortedEventsSnapShot){free(this->dataUnsortedEventsSnapShot);}
  CALLOC(dataUnsortedEventsSnapShot, cl_uint, event_total*eventDataSize);
  
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < eventBufferCount; b++)
  {
    event_total = synapticEvents.getEventCount(queue, b, currentTimeSlot, 
      Data_SynapticEvents::RECENT);

    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Check if valid and copy event*/
      cl_uint neuron = 0, time = 0, weight = 0;
      synapticEvents.getEvent(queue, b, e, currentTimeSlot, Data_SynapticEvents::RECENT, 
        neuron, time, weight);
      
      if(*((cl_float *)(&time)) < 0.0f || 
         *((cl_float *)(&time)) > (float)(maxDelay - minDelay))
      {
        ss << "Operator_Sort::captureUnsortedEvents: found an event with delay outside of"
          << " defined range: value " << *((cl_float *)(&time)) << ", range " << 0.0f 
          << " - " << (float)(maxDelay - minDelay) << std::endl;
        result = 0;
        break;
      }
      
      this->dataUnsortedEventsSnapShot[ptr_store] = neuron;
      this->dataUnsortedEventsSnapShot[ptr_store+1] = time;
      this->dataUnsortedEventsSnapShot[ptr_store+2] = weight;
      ptr_store += eventDataSize;
    }
  }
  
  if(!result)
  {
    free(this->dataUnsortedEventsSnapShot);
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/
#endif



#if SORT_VERIFY_ENABLE
#define PRINT_VERIFY_SORTED_EVENTS  0
void 
Operator_Sort::verifySortedEvents
(
  cl_uint currentTimeStep,
  cl_uint level,
  cl_uint eventDataSize,
  cl_uint eventDataDstBufMaxSize,
  cl_uint ptrStructCount,
  cl_uint ptrStructSize,
  cl_uint ptrStructItemSize,
  cl_uint ptrMode,
  double  maxDelay,
  double  minDelay,
  cl_uint *sortedEvents, 
  cl_uint *pointerStruct
)
/**************************************************************************************************/
{
  bool result = 1;
  cl_uint keyOffset = 0;
  cl_uint val1Offset = (keyOffset+1)%eventDataSize;
  cl_uint val2Offset = (keyOffset+2)%eventDataSize;
  char *sortedEventsCheck = NULL;
  std::stringstream ss;
  
  if(level > 0)
  {
    sortedEventsCheck = (char *)calloc(dataUnsortedEventsSnapShotSize/
      eventDataSize, sizeof(char));
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
    std::cout << "Operator_Sort::verifySortedEvents " << currentTimeStep << ": started" 
      << std::endl;
    std::cout << "Operator_Sort::verifySortedEvents " << currentTimeStep << ": found " 
      << dataUnsortedEventsSnapShotSize/eventDataSize 
      << " key-value(s) elements" << std::endl;
#endif
  
  /*Check sort order*/
  for
  (
    cl_uint i = 1; 
    i < dataUnsortedEventsSnapShotSize/eventDataSize; 
    i++
  ){
    /*Verify sorted order for keys*/
    cl_uint v1 = sortedEvents[keyOffset*eventDataDstBufMaxSize + i];
    cl_uint v2 = sortedEvents[keyOffset*eventDataDstBufMaxSize + i-1];
    
    if(v1 < v2)
    {
      ss << "Operator_Sort::verifySortedEvents " << currentTimeStep 
        << ": failed to verify sort order for key element " << i << ": " << v1 << "<" << v2 
        << std::endl;
      result = 0;
      break;
    }
    
    /*Verify sorted order for values*/
    if(v1 == v2)
    {
      cl_float time1 = 
        *((cl_float *)(&sortedEvents[val1Offset*eventDataDstBufMaxSize + i]));
      cl_float time2 = 
        *((cl_float *)(&sortedEvents[val1Offset*eventDataDstBufMaxSize + i-1]));
      
      if(time1 < 0.0f || time1 > (float)(maxDelay-minDelay))
      {
        ss << "Operator_Sort::verifySortedEvents: found an event in actual data with delay "
          << "outside of defined range: value " << time1 << ", range " << 0.0f << " - " 
          << (float)(maxDelay-minDelay) << std::endl;
        result = 0;
        break;
      }
      
      if(time1 < time2)
      {
        ss << "Operator_Sort::verifySortedEvents " << currentTimeStep 
          << ": failed to verify sort order for value element " << i << ": " << time1 << "<" 
          << time2 << std::endl;
        result = 0;
        break;
      }
    }
    /*
    std::cout  << v1 << "->" 
      << *((cl_float *)(&sortedEvents[val1Offset*eventDataDstBufMaxSize + i])) 
      << ", ";
    */
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
  if(result)
  {
    std::cout << "Operator_Sort::verifySortedEvents " << currentTimeStep 
      << ": verified sort order" << std::endl;
  }
#endif

  if(result && level > 0)
  {
#if PRINT_VERIFY_SORTED_EVENTS
    if(ptrMode == 1)
    {
      /*Determine how many avtive structs*/
      cl_uint activeStructsCount = 0xFFFFFFFF;
      for(cl_uint i = 0; i < ptrStructCount; i++)
      {
        if(pointerStruct[ptrStructSize*i] == 0)
        {
          activeStructsCount = i;
          break;
        }
      }

      cl_uint printFreq = 100000;
      std::cout << "Operator_Sort::verifySortedEvents " << currentTimeStep 
        << ": currently verifying key-value(s)/" << printFreq << ": ";
    }
#endif
    /*Each key-value(s) elements of this->dataUnsortedEventsSnapShot must be found in sortedEvents*/
    for
    (
      cl_uint p = 0; 
      p < dataUnsortedEventsSnapShotSize/eventDataSize; 
      p++
    ){
      cl_uint testedElementKey = 
        this->dataUnsortedEventsSnapShot[keyOffset + p*eventDataSize];
      cl_uint entryCount = 0;
      cl_uint entryAddress = 0xFFFFFFFF;
      cl_uint nextEntryAddress = 0xFFFFFFFF;

      if(ptrMode == 0)
      {
        entryAddress = *(pointerStruct + testedElementKey*ptrStructItemSize);
        entryCount = *(pointerStruct + testedElementKey*ptrStructItemSize + 1);
        nextEntryAddress = entryAddress + entryCount;
      }
      else if(ptrMode == 1)
      {
        cl_uint structPtr = 0;
        
        /*Determine to which struct testedElementKey belongs to*/
        /*TODO: binary search on sorted array*/
        for(structPtr = 0; structPtr < ptrStructCount; structPtr++)
        {
          entryCount = pointerStruct[ptrStructSize*structPtr];
          if(entryCount > 0)
          {
            cl_uint ptr = 
              /*base*/
              ptrStructSize*structPtr + 
              /*place for count*/
              1 + 
              /*last element*/
              (entryCount-1)*ptrStructItemSize;
            if(pointerStruct[ptr] >= testedElementKey){break;}
          }
        }

        /*Search for entry for this element*/
        /*TODO: binary search on sorted array*/
        for(cl_uint s = 0; s < entryCount; s++)
        {
          cl_uint ptr = 
            /*base*/
            ptrStructSize*structPtr + 
            /*place for count*/
            1 + 
            /*key of the element*/
            s*ptrStructItemSize;
          if(pointerStruct[ptr] == testedElementKey)
          {
            /*address*/
            entryAddress = pointerStruct[ptr+1];
            /*limit address for this element*/
            if(s < entryCount-1)
            {
              nextEntryAddress = pointerStruct[ptr+ptrStructItemSize + 1];
            }
            else
            {
              nextEntryAddress = pointerStruct[ptrStructSize*(structPtr + 1) + 
                ptrStructItemSize]; 
            }
            break;
          }
        }

        /*Verify that the entry was found*/
        if(entryAddress == 0xFFFFFFFF || nextEntryAddress == 0xFFFFFFFF)
        {
          ss << "Operator_Sort::verifySortedEvents " << currentTimeStep << ": not able to find "
            << "a struct for element " << testedElementKey << " at pointer " << p << std::endl;
          result = 0;
          break;
        }
      }

      /*Find the entry*/
      result = 0;
      /*TODO: binary search on sorted array*/
      for
      (
        cl_uint a = entryAddress; 
        a < nextEntryAddress; 
        a++
      ){
        if(testedElementKey == 
           sortedEvents[keyOffset*eventDataDstBufMaxSize + a])
        {
          if(this->dataUnsortedEventsSnapShot[val1Offset + p*eventDataSize] 
             == sortedEvents[val1Offset*eventDataDstBufMaxSize + a])
          {
            if(this->dataUnsortedEventsSnapShot[val2Offset + p*eventDataSize] 
               == sortedEvents[val2Offset*eventDataDstBufMaxSize + a])
            {
              result = 1;
              sortedEventsCheck[a] = 1;
              break;
            }
          }
        }
      }

      if(!result)
      {
        ss << "Operator_Sort::verifySortedEvents " << currentTimeStep << ": not able to find in "
          << "sorted data a combination of " << "key-value(s) for key " << p << "(" 
          << this->dataUnsortedEventsSnapShot[keyOffset + p*eventDataSize] 
          << ") from unsorted data" << std::endl;
        break;
      }
#if PRINT_VERIFY_SORTED_EVENTS
      if((p)%printFreq == 0)
      {
        std::cout  << ((p)/printFreq) << ",";
      }
#endif
    }
#if PRINT_VERIFY_SORTED_EVENTS
    if(result)
    {
      std::cout 
        << "\nverifySortedEvents" << currentTimeStep << ": verified mapping of unsorted "
        << "data in sorted data" << std::endl;
    }
#endif
  }
  
  /*Check for unrecognized keys*/
  if(result && level > 1)
  {
    for(cl_uint i = 0; i < dataUnsortedEventsSnapShotSize/eventDataSize; 
      i++)
    {
      if(!sortedEventsCheck[i])
      {
        ss << "Operator_Sort::verifySortedEvents " << currentTimeStep 
        << ": found unrecognized element (" << i << "): (" 
        << sortedEvents[keyOffset*eventDataDstBufMaxSize + i] 
        << ")->(" << sortedEvents[val1Offset*eventDataDstBufMaxSize + i] 
        << "," << sortedEvents[val2Offset*eventDataDstBufMaxSize + i] << ")" 
        << std::endl;
        result = 0;
        break;
      }
    }
#if PRINT_VERIFY_SORTED_EVENTS
    if(result)
    {
      std::cout << "Operator_Sort::verifySortedEvents " << currentTimeStep 
        << ": checked for unrecognized keys" << std::endl;
    }
#endif
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
  std::cout << "Operator_Sort::verifySortedEvents " << currentTimeStep << ": finished" 
    << std::endl;
#endif

  if(sortedEventsCheck){free(sortedEventsCheck);}
  if(this->dataUnsortedEventsSnapShot){free(this->dataUnsortedEventsSnapShot);} 
  this->dataUnsortedEventsSnapShot = NULL;

  if(!result)
  {
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/
#endif



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Operator_Sort::initialize
(
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || ENABLE_UNIT_TEST_SCAN_V00 || \
  ENABLE_UNIT_TEST_SCAN_V01) || (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
  cl_bool                             block,
  cl::CommandQueue                    &queue,
#endif
  cl::Context                         &context,
  cl::Device                          &device,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
#if ENABLE_OPERATOR_SCAN
  this->operatorScan = new Operator_Scan
  (
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || ENABLE_UNIT_TEST_SCAN_V00 || \
  ENABLE_UNIT_TEST_SCAN_V01)
    block,
    queue,
#endif
    context,
    device,
    kernelStats,
    this->dataToSimulationLogFile,
    this->dataToReportLogFile
  );
#endif

#if ENABLE_OPERATOR_GROUP
  this->operatorGroup = new Operator_Group
  (
#if (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
    block,
    queue,
#endif
    context,
    device,
    kernelStats,
    this->dataToSimulationLogFile,
    this->dataToReportLogFile
  );
#endif
}
/**************************************************************************************************/



#if (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
void
Operator_Sort::sortByTimeFirstPass
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
  cl_uint                   currentTimeSlot,
#endif
  cl_uint                   currentTimeStep,
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt
)
/**************************************************************************************************/
{
#if (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00)
  {
#if ENABLE_UNIT_TEST_SCAN_V00

    (*this->operatorScan).scanUnitTest_v00
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt,
      currentTimeStep
    );
    
#else

#if SCAN_VERIFY_ENABLE
    (*this->synapticEvents).refreshUnsortedHistogram(queue, CL_FALSE);
#endif

    (*this->operatorScan).scan_v00
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt,
      currentTimeStep,
      SAFE_GET((*this->synapticEvents).timeSlots),
      (*this->synapticEvents).dataUnsortedEventsHistogramBuffer
    );

#if DEVICE_HOST_DATA_COHERENCE
    (*this->synapticEvents).invalidateUnsortedHistogram();
#endif

#if SCAN_VERIFY_ENABLE
    cl_uint scanTimeStep = currentTimeStep % SAFE_GET((*this->synapticEvents).timeSlots);
    
    (*this->operatorScan).verifyScan_v00
    (
      queue,
      (void*)this->synapticEvents, 
      Data_SynapticEvents::getPreviousUnsortedHistogramItem, 
      Data_SynapticEvents::getCurrentUnsortedHistogramItem,
      SAFE_GET((*this->synapticEvents).timeSlots),
      SAFE_GET((*this->synapticEvents).histogramBinCount),
      SAFE_GET((*this->synapticEvents).histogramBinSize),
      scanTimeStep
    );
#endif
#endif
  }
#endif

#if (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
  {
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V00
    (*this->operatorGroup).groupUnitTest_v00
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt,
      currentTimeStep,
      (int)(this->groupTestMode),
      (*this->synapticEvents)
    );

#else

    (*this->operatorGroup).group_v00
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt,
      currentTimeStep,
      (*this->synapticEvents)
    );
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    (*this->operatorGroup).verifyGroup_v00
    (
      queue,
      0,
      currentTimeSlot,
      (*this->synapticEvents)
    );
#endif
  }
#endif
}
/**************************************************************************************************/
#endif



#if ((ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00) || \
    (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V02))
void
Operator_Sort::sortByTime
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
  cl_uint                   currentTimeSlot,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE || \
(ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt
)
/**************************************************************************************************/
{
#if (ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
  this->sortByTimeFirstPass
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    currentTimeSlot,
#endif
    currentTimeStep,
    queue,
    ndrEvt
  );
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02
  /*
    Sort positive floats (event time) in the rage 0.0f - 1.0f 0x(00 00 00 00 - 3F 80 00 00).
    No need to flip - unflip: 
    #define FLOAT_UNFLIP(f)       (f ^ (((f >> 31) - 1) | 0x80000000))
    #define FLOAT_FLIP(f)         (f ^ (-int(f >> 31) | 0x80000000))
  */
  for
  (
    cl_uint sortStep = 1; 
    sortStep < (this->sortByTimeIterationCount); 
    sortStep++
  )
#endif
  {
#if ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01
    {
#if !(GROUP_EVENTS_ENABLE_V00)

    (*this->operatorScan).scanUnitTest_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt
    );
    
#else

#if SCAN_VERIFY_ENABLE
    (*this->synapticEvents).refreshSortedEventsHistogram(queue, CL_TRUE);
#endif

    (*this->operatorScan).scan_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt,
      (*this->synapticEvents).dataSortedEventsHistogramReadBuffer
    );
    
#if DEVICE_HOST_DATA_COHERENCE
    (*this->synapticEvents).invalidateSortedHistogram();
#endif

#if SCAN_VERIFY_ENABLE
    (*this->operatorScan).verifyScan_v01
    (
      queue,
      (void*)this->synapticEvents,
      Data_SynapticEvents::getPreviousSortedHistogramItem,
      Data_SynapticEvents::getCurrentSortedHistogramItem
    );
#endif
#endif
    }
#endif

#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02
    if(sortStep < (this->sortByTimeIterationCount)-1)
#endif
    {
#if ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01
    {
#if !(GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01)

    (*this->operatorGroup).groupUnitTest_v01
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (this->groupTestMode),
      (*this->synapticEvents)
    );

#else

    (*this->operatorGroup).group_v01
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (*this->synapticEvents)
    );
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    (*this->operatorGroup).verifyGroup_v01
    (
      queue,
      currentTimeSlot,
      sortStep,
      (*this->synapticEvents)
    );
#endif
    }
#endif
    } /*if*/
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02
    else
#endif
    {
#if ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V02
    {
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V02

    (*this->operatorGroup).groupUnitTest_v02
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (this->groupTestMode),
      (*this->synapticEvents)
    );

#else

    (*this->operatorGroup).group_v02
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (*this->synapticEvents)
    );
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    (*this->operatorGroup).verifyGroup_v02
    (
      queue,
      currentTimeSlot,
      sortStep,
      (*this->synapticEvents)
    );
#endif
    }
#endif
    } /*else*/
  } /*for*/
}
/**************************************************************************************************/
#endif



#if ((ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01) || \
    (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V03))
void
Operator_Sort::sortByNeuron
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics   &kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
  cl_uint                   currentTimeSlot,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE
  cl_uint                   currentTimeStep,
#endif
  cl::CommandQueue          &queue,
  cl::Event                 &ndrEvt
)
/**************************************************************************************************/
{
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V03
  /*
    Sort positive integers (neuron IDs).
  */
  for
  (
    cl_uint sortStep = 0; 
    sortStep < (this->sortByNeuronIterationCount); 
    sortStep++
  )
#endif
  {
#if ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01
    {
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02)

    (*this->operatorScan).scanUnitTest_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt
    );
    
#else

#if SCAN_VERIFY_ENABLE
    (*this->synapticEvents).refreshSortedEventsHistogram(queue, CL_TRUE);
#endif

    (*this->operatorScan).scan_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      queue,
      ndrEvt,
      (*this->synapticEvents).dataSortedEventsHistogramReadBuffer
    );
    
#if DEVICE_HOST_DATA_COHERENCE
    (*this->synapticEvents).invalidateSortedHistogram();
#endif

#if SCAN_VERIFY_ENABLE
    (*this->operatorScan).verifyScan_v01
    (
      queue,
      (void*)this->synapticEvents,
      Data_SynapticEvents::getPreviousSortedHistogramItem,
      Data_SynapticEvents::getCurrentSortedHistogramItem
    );
#endif
#endif
    }
#endif
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V03
    if(sortStep < ((this->sortByNeuronIterationCount)-1))
#endif
    {
#if ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V01
    {
#if !(GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02)
    (*this->operatorGroup).groupUnitTest_v01
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (this->groupTestMode),
      (*this->synapticEvents)
    );

#else

    (*this->operatorGroup).group_v01
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (*this->synapticEvents)
    );
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    (*this->operatorGroup).verifyGroup_v01
    (
      queue,
      currentTimeSlot,
      sortStep,
      (*this->synapticEvents)
    );
#endif
    }
#endif
    }
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V03
    else
#endif
    {
#if ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V03
    {
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V03

    (*this->operatorGroup).groupUnitTest_v03
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (this->groupTestMode),
      (*this->synapticEvents)
    );

#else

    (*this->operatorGroup).group_v03
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if KERNEL_LEVEL_PROFILING || GROUP_EVENTS_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
      queue,
      ndrEvt,
      sortStep,
      (*this->synapticEvents)
    );
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    (*this->operatorGroup).verifyGroup_v03
    (
      queue,
      currentTimeSlot,
      sortStep,
      (*this->synapticEvents)
    );
#endif
    }
#endif
    } /*if*/
  } /*for*/
}
/**************************************************************************************************/
#endif

#endif  /*ENABLE_OPERATOR_SORT*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
