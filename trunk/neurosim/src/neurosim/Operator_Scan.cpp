
/* ===============================================================================================



  =============================================================================================== */

  
  
/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Operator_Scan.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



#if ENABLE_OPERATOR_SCAN



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



#if SCAN_DEBUG_ENABLE
void 
Operator_Scan::invalidateDebug
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_SCAN_VALID_DEBUG ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if SCAN_ERROR_TRACK_ENABLE
void 
Operator_Scan::invalidateError
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_SCAN_VALID_ERROR ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V00)
void 
Operator_Scan::invalidateUnitTestData_v00
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V01)
void 
Operator_Scan::invalidateUnitTestData_v01
()
/**************************************************************************************************/
{
  this->dataValid &= (OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if SCAN_ENABLE_V00
void 
Operator_Scan::scan_v00
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt,
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  cl::Buffer                          &dataScanBuffer
)
/**************************************************************************************************/
{
  cl_uint relativeTimeStep = currentTimeStep % timeSlotCount;

#if (SCAN_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v00)
  {
    this->setKernelArguments_v00 = false;
    
#if (SCAN_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelScanHistogramV00, this->dataScanDebugHostBuffer, 
      this->argNumScan00++);
    SET_KERNEL_ARG(this->kernelScanHistogramV00, this->dataScanDebugDeviceBuffer, 
      this->argNumScan00++);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelScanHistogramV00, this->dataScanErrorBuffer, 
      this->argNumScan00++);
#endif

    SET_KERNEL_ARG(this->kernelScanHistogramV00, dataScanBuffer, this->argNumScan00++);
  }
  
  SET_KERNEL_ARG(this->kernelScanHistogramV00, relativeTimeStep, this->argNumScan00);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelScanHistogramV00, 
    *(this->globalThreadsScanV00), *(this->localThreadsScanV00), NULL, ndrEvt, 
    "kernelScanHistogramV00");

#if (SCAN_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Scan::scan_v00:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_ERROR);}, 
    this->dataScanError
  );
#endif
}
/**************************************************************************************************/
#endif



#if SCAN_ENABLE_V01
void 
Operator_Scan::scan_v01
(
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
  cl_uint                             currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt,
  cl::Buffer                          &dataScanBuffer
)
/**************************************************************************************************/
{
#if (SCAN_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v01)
  {
    this->setKernelArguments_v01 = false;
    
#if (SCAN_DEBUG_ENABLE)
    SET_KERNEL_ARG(this->kernelScanHistogramV01, this->dataScanDebugHostBuffer, 
      this->argNumScan01++);
    SET_KERNEL_ARG(this->kernelScanHistogramV01, this->dataScanDebugDeviceBuffer, 
      this->argNumScan01++);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG(this->kernelScanHistogramV01, this->dataScanErrorBuffer, 
      this->argNumScan01++);
#endif
  }
  
  SET_KERNEL_ARG(this->kernelScanHistogramV01, dataScanBuffer, this->argNumScan01);
  
  ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, queue, this->kernelScanHistogramV01, 
    *(this->globalThreadsScanV01), *(this->localThreadsScanV01), NULL, ndrEvt, 
    "kernelScanHistogramV01");

#if (SCAN_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  THROW_ERROR_TRACK
  (
    "Operator_Scan::scan_v01:",
    currentTimeStep, 
    {this->invalidateError(); this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_ERROR);}, 
    this->dataScanError
  );
#endif
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V00)
void 
Operator_Scan::scanUnitTest_v00
(
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt,
  cl_uint                             currentTimeStep
)
/**************************************************************************************************/
{
  LOG_SIM("Operator_Scan::scanUnitTest_v00: performing unit test");
  
  cl_uint relativeTimeStep = currentTimeStep % (this->timeSlots);
  
  this->setUnitTestData_v00(queue, CL_FALSE, relativeTimeStep);

  this->scan_v00
  (
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
    queue,
    ndrEvt,
    currentTimeStep,
    this->timeSlots,
    this->dataHistogramV00Buffer
  );
  
  this->invalidateUnitTestData_v00();

#if (SCAN_ENABLE_V00 && SCAN_VERIFY_ENABLE)
  this->verifyScan_v00
  (
    queue,
    (void*)this, 
    Operator_Scan::getPreviousHistogramItem_v00, 
    Operator_Scan::getCurrentHistogramItem_v00,
    this->timeSlots,
    this->histogramBinCount_v00,
    this->histogramBinSize_v00,
    relativeTimeStep
  );
#endif
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V01)
void 
Operator_Scan::scanUnitTest_v01
(
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
  cl_uint                             currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
  struct kernelStatistics             &kernelStats,
#endif
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
  LOG_SIM("Operator_Scan::scanUnitTest_v01: performing unit test");
  
  this->setUnitTestData_v01(queue, CL_TRUE);

  this->scan_v01
  (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
    currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
    kernelStats,
#endif
    queue,
    ndrEvt,
    this->dataHistogramV01Buffer
  );
  
  this->invalidateUnitTestData_v01();

#if (SCAN_ENABLE_V01 && SCAN_VERIFY_ENABLE)
  this->verifyScan_v01
  (
    queue,
    (void*)this, 
    Operator_Scan::getPreviousHistogramItem_v01, 
    Operator_Scan::getCurrentHistogramItem_v01
  );
#endif
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_V00 && SCAN_VERIFY_ENABLE)
void 
Operator_Scan::verifyScan_v00
(
  cl::CommandQueue  &queue,
  void              *objectToCall,
  cl_uint           (*getPreviousItem)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
  cl_uint           (*getResentItem)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
  cl_uint           timeSlots,
  cl_uint           histogramBinCount,
  cl_uint           histogramBinSize,
  cl_uint           timeSlot
)
/**************************************************************************************************/
{
  cl_uint runningSum = 0, element = 0;
  
  cl_uint size = (timeSlots)*((histogramBinCount)*(histogramBinSize) + 1);
  cl_uint *dataHistogramTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  for(cl_uint binID = 0; binID < histogramBinCount; binID++)
  {
    /*Pointer is based on time slot, bin, WG*/
    cl_uint offset = 
    /*Time slot space = (time slot #) x (bins per WG) x WGs*/
    timeSlot*((histogramBinCount)*(histogramBinSize) + 1) + 
    /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
    binID*(histogramBinSize) /*+ bufferID*/;
    
    for(cl_uint bufferID = 0; bufferID < histogramBinSize; bufferID++)
    {

      element = getPreviousItem(queue, timeSlot, binID, bufferID, objectToCall);

      if(binID == 0 && bufferID == 0)
      {
        runningSum = element;
        dataHistogramTemp[offset + bufferID] = 0;
      }
      else
      {
        dataHistogramTemp[offset + bufferID] = runningSum;
        runningSum += element;
      }
    }
  }
  
  dataHistogramTemp[timeSlot*((histogramBinCount)*(histogramBinSize) + 1) + 
    (histogramBinCount)*(histogramBinSize)] = runningSum;

  cl_uint errorEventTotals = 0; 
  cl_uint offsetTimeSlot = timeSlot*((histogramBinCount)*(histogramBinSize) + 1);

  for(cl_uint binID = 0; binID < histogramBinCount; binID++)
  {
    for(cl_uint bufferID = 0; bufferID < histogramBinSize; bufferID++)
    {
      /*Pointer is based on time slot, bin, WG*/
      cl_uint offset = 
      /*Time slot space = (time slot #) x (bins per WG) x WGs*/
      offsetTimeSlot + 
      /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
      binID*(histogramBinSize) + bufferID;

      element = getResentItem(queue, timeSlot, binID, bufferID, objectToCall);
      
      if(dataHistogramTemp[offset] != element)
      {
        errorEventTotals++;
        /*
        std::cout << binID << "->(" << dataHistogramTemp[offset] << "," << element << "), ";
        */
      }
    }
  }
  
  if(dataHistogramTemp)
    free(dataHistogramTemp);

  if(errorEventTotals)
  {
    THROW_SIMEX("verifyScan_v00: failed to match verification data " << errorEventTotals 
      << " times in time slot " << timeSlot);
  }
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_V01 && SCAN_VERIFY_ENABLE)
void 
Operator_Scan::verifyScan_v01
(
  cl::CommandQueue  &queue,
  void              *objectToCall,
  cl_uint           (*getPreviousItem)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
  cl_uint           (*getResentItem)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*)
)
/**************************************************************************************************/
{
#if (ENABLE_UNIT_TEST_SCAN_V01)
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);
#endif

  /* allocate memory for histogram */
  cl_uint *temp = (cl_uint *)calloc(((this->histogramBinSize_v01)*(this->histogramBinCount_v01)+1), 
    sizeof(cl_uint));
  
  /*Compute histogram*/
  for(cl_uint j = 0; j < (this->histogramBinSize_v01); j++)
  {
    for(cl_uint b = 0; b < (this->histogramBinCount_v01); b++)
    {
      cl_uint sum = 0;
      
      for(cl_uint w = 0; w < (this->histogramBinBackets); w++)
      {
        sum += getPreviousItem(queue, w, b, j, objectToCall);
      }
      temp[j + b*(this->histogramBinSize_v01)] = sum;
    }
  }
  
  cl_uint runningSum = temp[0];
  temp[0] = 0;
  
  /*Compute offsets*/
  for(cl_uint j = 1; j < ((this->histogramBinCount_v01)*(this->histogramBinSize_v01) + 1); j++)
  {
    cl_uint d = temp[j];
    temp[j] = runningSum;
    runningSum += d;
  }
  
  unsigned long long error_count_histogram_out = 0;

  for(cl_uint i = 1; i < ((this->histogramBinCount_v01)*(this->histogramBinSize_v01) + 1); i++)
  {
    cl_uint w = (i/((this->histogramBinCount_v01)*(this->histogramBinSize_v01)))%
      (this->histogramBinBackets);
    cl_uint b = (i/(this->histogramBinSize_v01))%(this->histogramBinCount_v01);
    cl_uint j = i%(this->histogramBinSize_v01);
    
    error_count_histogram_out += (getResentItem(queue, w, b, j, objectToCall) != temp[i]);
    /*
    std::cout << j << "->(" << dataHistogram[j] << "," << temp[j] << "), ";
    */
  }
  
  if(temp)
    free(temp);

  /*Verify*/
  if(error_count_histogram_out)
  {
    THROW_SIMEX("verifyScan_v01: failed to match 01 histogram " << error_count_histogram_out 
      << " times");
  }
}
/**************************************************************************************************/
#endif



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



void
Operator_Scan::initialize
(
#if (SCAN_DEBUG_ENABLE)
  size_t                              debugBufferSizeWords,
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  size_t                              errorBufferSizeWords,
#endif
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE)
  cl::CommandQueue                    &queue,
  cl_bool                             block,
#endif
  cl::Context                         &context,
  cl::Device                          &device,
  struct kernelStatistics             &kernelStats,
  size_t                              cacheSizeWords
)
/**************************************************************************************************/
{
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
    QueryPerformanceFrequency((LARGE_INTEGER *)&(this->performanceFrequency));
#endif

  /* register device local memory buffer for stats */
  size_t lmScanData = sizeof(cl_uint)*cacheSizeWords; 
  size_t lmScanDataSizeBytes = lmScanData;
  REGISTER_MEMORY_O(device, SCAN_KERNEL_NAME, MEM_LOCAL, lmScanData, kernelStats);

#if (SCAN_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataScanDebugHost, cl_uint, debugBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataScanDebugHost, kernelStats);
  /* allocate memory for debug device buffer */
  CALLOC(dataScanDebugDevice, cl_uint, debugBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataScanDebugDevice, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataScanDebugHostBuffer, 
    this->dataScanDebugHostSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataScanDebugDeviceBuffer, 
    this->dataScanDebugDeviceSizeBytes);
    
  this->storeData(queue, block, OPERATOR_SCAN_VALID_DEBUG);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataScanError, cl_uint, errorBufferSizeWords);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataScanError, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataScanErrorBuffer, 
    this->dataScanErrorSizeBytes);
    
  this->storeData(queue, block, OPERATOR_SCAN_VALID_ERROR);
#endif

#if SCAN_ENABLE_V00
  createKernel
  (
#if LOG_SIMULATION
    this->dataToSimulationLogFile,
#endif
    context,
    device,
    this->kernelScanHistogramV00,
    SCAN_KERNEL_FILE_NAME,
    SCAN_KERNEL_NAME,
    "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D SCAN_DEVICE_V00=1",
    this->blockSizeX_kernelScanHistogramV00,
    this->blockSizeY_kernelScanHistogramV00
  );
#endif

#if SCAN_ENABLE_V01
  createKernel
  (
#if LOG_SIMULATION
    this->dataToSimulationLogFile,
#endif
    context,
    device,
    this->kernelScanHistogramV01,
    SCAN_KERNEL_FILE_NAME,
    SCAN_KERNEL_NAME,
    "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D SCAN_DEVICE_V01=1",
    this->blockSizeX_kernelScanHistogramV01,
    this->blockSizeY_kernelScanHistogramV01
  );
#endif
}
/**************************************************************************************************/



#if (ENABLE_UNIT_TEST_SCAN_V00)
void
Operator_Scan::setupUnitTest_v00
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
  size_t size = (this->timeSlots)*((this->histogramBinCount_v00)*(this->histogramBinSize_v00) + 1);
  
  CALLOC(dataHistogramV00, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataHistogramV00, kernelStats);

  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataHistogramV00Buffer, 
    this->dataHistogramV00SizeBytes);
    
  this->dataPastHistogramV00 = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM);
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V01)
void
Operator_Scan::setupUnitTest_v01
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             &kernelStats
)
/**************************************************************************************************/
{
  size_t size = 
    (this->histogramBinBackets)*(this->histogramBinCount_v01)*(this->histogramBinSize_v01);
  
  CALLOC(dataHistogramV01, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataHistogramV01, kernelStats);

  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataHistogramV01Buffer, 
    this->dataHistogramV01SizeBytes);
    
  this->dataPastHistogramV01 = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V00)
void
Operator_Scan::setUnitTestData_v00
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             timeSlot
)
/**************************************************************************************************/
{
  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG_SIM("Operator_Scan::setUnitTestData_v00: set srand seed to " << this->srandSeed);
  
  cl_uint maxElementSize = (0xFFFFFFFF/((this->histogramBinCount_v00)*(this->histogramBinSize_v00)));
  
  for(cl_uint j = 0; j < (this->histogramBinCount_v00); j++)
  {
    cl_uint offset = timeSlot*((this->histogramBinCount_v00)*(this->histogramBinSize_v00) + 1) + 
      j*(this->histogramBinSize_v00);
    
    for(cl_uint k = 0; k < (this->histogramBinSize_v00); k++)
    {
      this->dataHistogramV00[offset + k] = 
        cl_uint(abs(maxElementSize*((double)rand()/((double)RAND_MAX))));
    }
  }
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM);
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V01)
void
Operator_Scan::setUnitTestData_v01
(
  cl::CommandQueue    &queue,
  cl_bool             block
)
/**************************************************************************************************/
{
  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG_SIM("Operator_Scan::setUnitTestData_v01: set srand seed to " << this->srandSeed);
  
  memset(this->dataHistogramV01, 0, this->dataHistogramV01SizeBytes);
  
  for(cl_uint j = 0; j < (this->histogramBinSize_v01); j++)
  {
    for(cl_uint b = 0; b < (this->histogramBinCount_v01); b++)
    {
      for(cl_uint w = 0; w < (this->histogramBinBackets); w++)
      {
        cl_uint p = 
          /*WG offset*/
          w*((this->histogramBinSize_v01)*(this->histogramBinCount_v01)) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*(this->histogramBinSize_v01);
          
        this->dataHistogramV01[p] = cl_uint(abs((this->histogramMaxCount)*((double)rand()/
          ((double)RAND_MAX))));
      }
    }
  }
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V00)
cl_uint
Operator_Scan::getPreviousHistogramItem_v00
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Operator_Scan* mySelf = (Operator_Scan*)objectToCallThisFunction;
  
  mySelf->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM);

  /*Pointer is based on time slot, bin, WG*/
  cl_uint ptr = 
  /*Time slot space = (time slot #) x (bins per WG) x WGs*/
  timeSlot*((mySelf->histogramBinCount_v00)*(mySelf->histogramBinSize_v00) + 1) + 
  /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
  binID*(mySelf->histogramBinSize_v00) + bufferID;
  
  return mySelf->dataPastHistogramV00[ptr];
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V00)
cl_uint
Operator_Scan::getCurrentHistogramItem_v00
(
  cl::CommandQueue  &queue,
  cl_uint           timeSlot,
  cl_uint           binID,
  cl_uint           bufferID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Operator_Scan* mySelf = (Operator_Scan*)objectToCallThisFunction;
  
  mySelf->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM);

  /*Pointer is based on time slot, bin, WG*/
  cl_uint ptr = 
  /*Time slot space = (time slot #) x (bins per WG) x WGs*/
  timeSlot*((mySelf->histogramBinCount_v00)*(mySelf->histogramBinSize_v00) + 1) + 
  /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
  binID*(mySelf->histogramBinSize_v00) + bufferID;
  
  return mySelf->dataHistogramV00[ptr];
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V01)
cl_uint
Operator_Scan::getPreviousHistogramItem_v01
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Operator_Scan* mySelf = (Operator_Scan*)objectToCallThisFunction;
  
  mySelf->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);

  cl_uint ptr = 
    /*WG offset*/
    backetID*((mySelf->histogramBinSize_v01)*(mySelf->histogramBinCount_v01)) + 
    /*WG offset for histogram out*/
    itemID +
    /*bin*/
    binID*(mySelf->histogramBinSize_v01);

  return mySelf->dataPastHistogramV01[ptr];
}
/**************************************************************************************************/
#endif



#if (ENABLE_UNIT_TEST_SCAN_V01)
cl_uint
Operator_Scan::getCurrentHistogramItem_v01
(
  cl::CommandQueue  &queue,
  cl_uint           backetID,
  cl_uint           binID,
  cl_uint           itemID,
  void              *objectToCallThisFunction
)
/**************************************************************************************************/
{
  Operator_Scan* mySelf = (Operator_Scan*)objectToCallThisFunction;
  
  mySelf->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);

  cl_uint ptr = 
    /*WG offset*/
    backetID*((mySelf->histogramBinSize_v01)*(mySelf->histogramBinCount_v01)) + 
    /*WG offset for histogram out*/
    itemID +
    /*bin*/
    binID*(mySelf->histogramBinSize_v01);
  
  return mySelf->dataHistogramV01[ptr];
}
/**************************************************************************************************/
#endif



#if SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || \
    ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01
void
Operator_Scan::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
#if SCAN_DEBUG_ENABLE
  IF_HIT_READ(selectBitMask, OPERATOR_SCAN_VALID_DEBUG, this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataScanDebugHostBuffer, 
      this->dataScanDebugHostSizeBytes, this->dataScanDebugHost);
      
    ENQUEUE_READ_BUFFER(block, queue, this->dataScanDebugDeviceBuffer, 
      this->dataScanDebugDeviceSizeBytes, this->dataScanDebugDevice);
      
    this->dataValid |= OPERATOR_SCAN_VALID_DEBUG;
  }
#endif

#if SCAN_ERROR_TRACK_ENABLE
  IF_HIT_READ(selectBitMask, OPERATOR_SCAN_VALID_ERROR, this->dataValid)
  {
    ENQUEUE_READ_BUFFER(block, queue, this->dataScanErrorBuffer, 
      this->dataScanErrorSizeBytes, this->dataScanError);
      
    this->dataValid |= OPERATOR_SCAN_VALID_ERROR;
  }
#endif

#if ENABLE_UNIT_TEST_SCAN_V00
  IF_HIT_READ(selectBitMask, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM, this->dataValid)
  {
    swap2(this->dataHistogramV00, this->dataPastHistogramV00);

    ENQUEUE_READ_BUFFER(block, queue, this->dataHistogramV00Buffer, 
      this->dataHistogramV00SizeBytes, this->dataHistogramV00);
      
    this->dataValid |= OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM;
  }
#endif

#if ENABLE_UNIT_TEST_SCAN_V01
  IF_HIT_READ(selectBitMask, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM, this->dataValid)
  {
    swap2(this->dataHistogramV01, this->dataPastHistogramV01);

    ENQUEUE_READ_BUFFER(block, queue, this->dataHistogramV01Buffer, 
      this->dataHistogramV01SizeBytes, this->dataHistogramV01);
      
    this->dataValid |= OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM;
  }
#endif
}
/**************************************************************************************************/
#endif



#if SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || \
    ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01
void
Operator_Scan::storeData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
#if SCAN_DEBUG_ENABLE
  if(selectBitMask & OPERATOR_SCAN_VALID_DEBUG)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataScanDebugHostBuffer, 
      this->dataScanDebugHostSizeBytes, this->dataScanDebugHost);
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataScanDebugDeviceBuffer, 
      this->dataScanDebugDeviceSizeBytes, this->dataScanDebugDevice);
    
    this->dataValid |= OPERATOR_SCAN_VALID_DEBUG;
  }
#endif

#if SCAN_ERROR_TRACK_ENABLE
  if(selectBitMask & OPERATOR_SCAN_VALID_ERROR)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataScanErrorBuffer, 
      this->dataScanErrorSizeBytes, this->dataScanError);
    
    this->dataValid |= OPERATOR_SCAN_VALID_ERROR;
  }
#endif

#if ENABLE_UNIT_TEST_SCAN_V00
  if(selectBitMask & OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataHistogramV00Buffer, 
      this->dataHistogramV00SizeBytes, this->dataHistogramV00);

    this->dataValid |= OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM;
  }
#endif

#if ENABLE_UNIT_TEST_SCAN_V01
  if(selectBitMask & OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM)
  {
    ENQUEUE_WRITE_BUFFER(block, queue, this->dataHistogramV01Buffer, 
      this->dataHistogramV01SizeBytes, this->dataHistogramV01);

    this->dataValid |= OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM;
  }
#endif
}
/**************************************************************************************************/
#endif



#endif  /*ENABLE_OPERATOR_SCAN*/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
