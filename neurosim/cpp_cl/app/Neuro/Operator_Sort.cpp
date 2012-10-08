
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



void 
Operator_Sort::clearData
(
  cl::CommandQueue  &queue,
  cl_bool           block,
  cl_uint           clearBitMask
)
/**************************************************************************************************/
{
  if(clearBitMask & 0x1)
  {
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
    swap2(dataHistogramGroupEventsTik, dataPastHistogramGroupEventsTik);
#endif

    memset(this->dataHistogramGroupEventsTik, 0, this->dataHistogramGroupEventsTikSizeBytes);
    
    this->storeData(queue, block, OPERATOR_SORT_VALID_ALL);
  }
  
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
  if(clearBitMask & 0x2)
  {
    memset(this->dataPastHistogramGroupEventsTik, 0, this->dataHistogramGroupEventsTikSizeBytes);
  }
#endif
}
/**************************************************************************************************/



cl_uint
Operator_Sort::getPreviousHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Operator_Sort* mySelf = (Operator_Sort*)objectToCallThisFunction;
  return mySelf->getHistogramItem(queue, backetID, binID, itemID, Operator_Sort::PREVIOUS);
}
/**************************************************************************************************/



cl_uint
Operator_Sort::getCurrentHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Operator_Sort* mySelf = (Operator_Sort*)objectToCallThisFunction;
  return mySelf->getHistogramItem(queue, backetID, binID, itemID, Operator_Sort::RECENT);
}
/**************************************************************************************************/



cl_uint
Operator_Sort::getHistogramItem
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  const int         type
)
/**************************************************************************************************/
{
#if CLASS_VALIDATION_ENABLE
  if(backetID >= this->histogramBacketCount)
  {
    throw SimException("Operator_Sort::getHistogramItem: backetID exceeds total backet count");
  }
  
  if(binID >= this->histogramBinCount)
  {
    throw SimException("Operator_Sort::getHistogramItem: binID exceeds bin count");
  }
  
  if(itemID >= this->histogramBinSize)
  {
    throw SimException("Operator_Sort::getHistogramItem: itemID exceeds bin size");
  }
#endif

  this->getData(queue, CL_TRUE, OPERATOR_SORT_VALID_HISTOGRAM_TIK);

  cl_uint ptr = 
    /*WG offset*/
    backetID*((this->histogramBinSize)*(this->histogramBinCount)) + 
    /*WG offset for histogram out*/
    itemID +
    /*bin*/
    binID*(this->histogramBinSize);
  
  if(type == Operator_Sort::RECENT)
  {
    return this->dataHistogramGroupEventsTik[ptr];
  }
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
  else if(type == Operator_Sort::PREVIOUS)
  {
    return this->dataPastHistogramGroupEventsTik[ptr];
  }
#endif
  else
  {
    throw SimException("Operator_Sort::getHistogramItem: incorrect item type");
  }
}
/**************************************************************************************************/



void 
Operator_Sort::invalidateHistogram
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_SORT_VALID_HISTOGRAM_TIK ^ 0xFFFFFFFF);
}
/**************************************************************************************************/



void 
Operator_Sort::refreshHistogram
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  this->getData
  (
    queue, 
    block, 
    OPERATOR_SORT_VALID_HISTOGRAM_TIK
  );
}
/**************************************************************************************************/



void
Operator_Sort::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  IF_HIT_READ(selectBitMask, OPERATOR_SORT_VALID_HISTOGRAM_TIK, this->dataValid)
  {
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
    swap2(dataHistogramGroupEventsTik, dataPastHistogramGroupEventsTik);
#endif
    
    ENQUEUE_READ_BUFFER(block, queue, this->dataHistogramGroupEventsTikBuffer, 
      this->dataHistogramGroupEventsTikSizeBytes, this->dataHistogramGroupEventsTik);
      
    this->dataValid |= OPERATOR_SORT_VALID_HISTOGRAM_TIK;
  }
}
/**************************************************************************************************/



void
Operator_Sort::storeData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
  if(selectBitMask & OPERATOR_SORT_VALID_HISTOGRAM_TIK)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataHistogramGroupEventsTikBuffer, 
      this->dataHistogramGroupEventsTikSizeBytes, this->dataHistogramGroupEventsTik);
      
    this->dataValid |= OPERATOR_SORT_VALID_HISTOGRAM_TIK;
  }
}
/**************************************************************************************************/



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Operator_Sort::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
  size_t size = (this->histogramBacketCount)*(this->histogramBinSize)*(this->histogramBinCount);
  
  CALLOC(dataHistogramGroupEventsTik, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataHistogramGroupEventsTik, kernelStats);
#if OPERATOR_SORT_ENABLE_PAST_EVENTS
  this->dataPastHistogramGroupEventsTik = (cl_uint *)calloc(size, sizeof(cl_uint));
#endif
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataHistogramGroupEventsTikBuffer, 
    this->dataHistogramGroupEventsTikSizeBytes);
  
  this->storeData(queue, block, OPERATOR_SORT_VALID_ALL);
}
/**************************************************************************************************/

#endif  /*ENABLE_OPERATOR_SORT*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
