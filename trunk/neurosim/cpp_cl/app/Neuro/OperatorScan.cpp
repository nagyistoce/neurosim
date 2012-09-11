
/* ===============================================================================================



  =============================================================================================== */


  
/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "OperatorScan.hpp"

/**************************************************************************************************/



/**************************************************************************************************/
/*  public:                                                                                       */
/**************************************************************************************************/



OperatorScan::OperatorScan()
/**************************************************************************************************/
{
  this->reset(false);
}
/**************************************************************************************************/



OperatorScan::~OperatorScan()
/**************************************************************************************************/
{
  this->reset(true);
}
/**************************************************************************************************/



void
OperatorScan::initialize
(
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             *kernelStats,
  std::stringstream                   *dataToSimulationLogFile,
  std::stringstream                   *dataToReportLogFile
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  if(!this->resetObject)
  {
    throw SimException("OperatorScan::initialize: attemp to initialize without reset");
  }
#endif
  
  this->resetObject = false;
  this->dataValid = 0;
  this->dataToSimulationLogFile = dataToSimulationLogFile;
  this->dataToReportLogFile = dataToReportLogFile;
  
#if (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)
  
  /* register device local memory buffer for stats */
  //TODO: verify if it's correct
  size_t lmScanData = sizeof(cl_uint)*SCAN_WG_SIZE_WI*2; 
  lmScanDataSizeBytes = lmScanData;
  REGISTER_MEMORY_O(device, SCAN_KERNEL_NAME, MEM_LOCAL, lmScanData, kernelStats);

#if (SCAN_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC_O(dataScanDebugHost, cl_uint, SCAN_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataScanDebugHost, kernelStats);
  /* allocate memory for debug device buffer */
  CALLOC_O(dataScanDebugDevice, cl_uint, SCAN_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataScanDebugDevice, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataScanDebugHostBuffer, 
    this->dataScanDebugHostSizeBytes);
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataScanDebugDeviceBuffer, 
    this->dataScanDebugDeviceSizeBytes);
    
  this->storeData(queue, block, OPERATOR_SCAN_VALID_DEBUG);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC_O(dataScanError, cl_uint, SCAN_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataScanError, kernelStats);
  
  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataScanErrorBuffer, 
    this->dataScanErrorSizeBytes);
    
  this->storeData(queue, block, OPERATOR_SCAN_VALID_ERROR);
#endif

#if SCAN_ENABLE_V00
  this->blockSizeX_kernelScanHistogramV00 = SCAN_WG_SIZE_WI;
  this->blockSizeY_kernelScanHistogramV00 = 1;
  this->argNumScan00 = 0;
  this->setKernelArguments_v00 = true;
  
  this->globalThreadsScanV00 = 
    new cl::NDRange(SCAN_WG_SIZE_WI*SCAN_GRID_SIZE_WG);
  this->localThreadsScanV00 = 
    new cl::NDRange(this->blockSizeX_kernelScanHistogramV00, 
    this->blockSizeY_kernelScanHistogramV00);
    
  createKernel
  (
    context,
    device,
    this->kernelScanHistogramV00,
    SCAN_KERNEL_FILE_NAME,
    SCAN_KERNEL_NAME,
    "-D SCAN_DEVICE_V00",
    this->blockSizeX_kernelScanHistogramV00,
    this->blockSizeY_kernelScanHistogramV00
  );
#endif

#if SCAN_ENABLE_V01
  this->blockSizeX_kernelScanHistogramV01 = SCAN_WG_SIZE_WI;
  this->blockSizeY_kernelScanHistogramV01 = 1;
  this->argNumScan01 = 0;
  this->setKernelArguments_v01 = true;
  
  this->globalThreadsScanV01 = 
    new cl::NDRange(SCAN_WG_SIZE_WI*SCAN_GRID_SIZE_WG);
  this->localThreadsScanV01 = 
    new cl::NDRange(this->blockSizeX_kernelScanHistogramV01, 
    this->blockSizeY_kernelScanHistogramV01);
    
  createKernel
  (
    context,
    device,
    this->kernelScanHistogramV01,
    SCAN_KERNEL_FILE_NAME,
    SCAN_KERNEL_NAME,
    "-D SCAN_DEVICE_V01",
    this->blockSizeX_kernelScanHistogramV01,
    this->blockSizeY_kernelScanHistogramV01
  );
#endif
#endif
}
/**************************************************************************************************/



#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_DEBUG_ENABLE)
void 
OperatorScan::invalidateDebug
()
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= (OPERATOR_SCAN_VALID_DEBUG ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_ERROR_TRACK_ENABLE)
void 
OperatorScan::invalidateError
()
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

  this->dataValid &= (OPERATOR_SCAN_VALID_ERROR ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_UNIT_TEST_V00)
void 
OperatorScan::invalidateUnitTestData_v00
()
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(this->initDataUnitTest_v00)
  {
    throw SimException("OperatorScan::invalidateUnitTestData_v00: attemp to ivalidate "
      "unit test data without setup");
  }
#endif

  this->dataValid &= (OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_UNIT_TEST_V01)
void 
OperatorScan::invalidateUnitTestData_v01
()
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(this->initDataUnitTest_v01)
  {
    throw SimException("OperatorScan::invalidateUnitTestData_v01: attemp to ivalidate "
      "unit test data without setup");
  }
#endif

  this->dataValid &= (OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM ^ 0xFFFFFFFF);
}
/**************************************************************************************************/
#endif



#if SCAN_ENABLE_V00
void 
OperatorScan::scan_v00
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
  cl_uint                             currentTimeStep,
  cl_uint                             timeSlotCount,
  cl::Buffer                          &dataScanBuffer,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

  cl_uint relativeTimeStep = currentTimeStep % timeSlotCount;

#if (SCAN_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v00)
  {
    this->setKernelArguments_v00 = false;
    
#if (SCAN_DEBUG_ENABLE)
    SET_KERNEL_ARG_O(this->kernelScanHistogramV00, this->dataScanDebugHostBuffer, 
      this->argNumScan00++);
    SET_KERNEL_ARG_O(this->kernelScanHistogramV00, this->dataScanDebugDeviceBuffer, 
      this->argNumScan00++);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG_O(this->kernelScanHistogramV00, this->dataScanErrorBuffer, 
      this->argNumScan00++);
#endif

    SET_KERNEL_ARG_O(this->kernelScanHistogramV00, dataScanBuffer, this->argNumScan00++);
  }
  
  SET_KERNEL_ARG_O(this->kernelScanHistogramV00, relativeTimeStep, this->argNumScan00);
  
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  double startAppTime = 0, endAppTime = 0;
  if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif

  ENQUEUE_KERNEL(this->kernelScanHistogramV00, *(this->globalThreadsScanV00), 
    *(this->localThreadsScanV00), NULL, ndrEvt);
  
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  if(currentTimeStep >= START_PROFILING_AT_STEP)
  {
    cl_uint status = ndrEvt.wait();
    endAppTime = timeStampNs();
    
    if(status != CL_SUCCESS)
    {
      std::stringstream ss;
      ss << "OperatorScan::scan_v00: Failed cl:Event.wait() due to error code " << status << "\n";
      throw SimException(ss.str());
    }
    REGISTER_TIME(kernelScanHistogramV00, (endAppTime-startAppTime), 1.0)
  }
#endif

#if (SCAN_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  this->invalidateError();

  if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
  {
    this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_ERROR);
    
    if(this->dataScanError[0] != 0)
    {
      std::stringstream ss;
      ss << "OperatorScan::scan_v00: received error code from the device: " 
        << this->dataScanError[0] << std::endl;
      throw SimException(ss.str());
    }
  }
#endif
}
/**************************************************************************************************/
#endif



#if SCAN_ENABLE_V01
void 
OperatorScan::scan_v01
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
  cl_uint                             currentTimeStep,
  cl::Buffer                          &dataScanBuffer,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

#if (SCAN_DEBUG_ENABLE)
  this->storeData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

  if(this->setKernelArguments_v01)
  {
    this->setKernelArguments_v01 = false;
    
#if (SCAN_DEBUG_ENABLE)
    SET_KERNEL_ARG_O(this->kernelScanHistogramV01, this->dataScanDebugHostBuffer, 
      this->argNumScan01++);
    SET_KERNEL_ARG_O(this->kernelScanHistogramV01, this->dataScanDebugDeviceBuffer, 
      this->argNumScan01++);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
    SET_KERNEL_ARG_O(this->kernelScanHistogramV01, this->dataScanErrorBuffer, 
      this->argNumScan01++);
#endif
  }
  
  SET_KERNEL_ARG_O(this->kernelScanHistogramV01, dataScanBuffer, this->argNumScan01);
  /*TODO: get rid of unnecessary argument with preprocessor*/
  SET_KERNEL_ARG_O(this->kernelScanHistogramV01, 0, this->argNumScan01+1);
  
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  double startAppTime = 0, endAppTime = 0;
  if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif

  ENQUEUE_KERNEL(this->kernelScanHistogramV01, *(this->globalThreadsScanV01), 
    *(this->localThreadsScanV01), NULL, ndrEvt);
  
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  if(currentTimeStep >= START_PROFILING_AT_STEP)
  {
    cl_uint status = ndrEvt.wait();
    endAppTime = timeStampNs();
    
    if(status != CL_SUCCESS)
    {
      std::stringstream ss;
      ss << "OperatorScan::scan_v01: Failed cl:Event.wait() due to error code " << status << "\n";
      throw SimException(ss.str());
    }
    REGISTER_TIME(kernelScanHistogramV01, (endAppTime-startAppTime), 1.0)
  }
#endif

#if (SCAN_DEBUG_ENABLE)
  this->invalidateDebug();
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_DEBUG);
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
  this->invalidateError();

  if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
  {
    this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_ERROR);
    
    if(this->dataScanError[0] != 0)
    {
      std::stringstream ss;
      ss << "OperatorScan::scan_v01: received error code from the device: " 
        << this->dataScanError[0] << std::endl;
      throw SimException(ss.str());
    }
  }
#endif
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_UNIT_TEST_V00)
void 
OperatorScan::scanUnitTest_v00
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
  cl_uint                             currentTimeStep,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(this->initDataUnitTest_v00)
  {
    throw SimException("OperatorScan::scanUnitTest_v00: attemp to run unit test without setup");
  }
#endif
  LOG("OperatorScan::scanUnitTest_v00: performing unit test", 0);
  
  cl_uint relativeTimeStep = currentTimeStep % (this->timeSlots);
  
  this->setUnitTestData_v00(queue, CL_FALSE, relativeTimeStep);

  this->scan_v00
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    kernelStats,
#endif
    currentTimeStep,
    this->timeSlots,
    this->dataHistogramV00Buffer,
    queue,
    ndrEvt
  );
  
  this->invalidateUnitTestData_v00();

#if SCAN_VERIFY_ENABLE
  this->verifyScan_v00
  (
    this->timeSlots,
    this->histogramBinCount_v00,
    this->histogramBinSize_v00,
    queue,
    relativeTimeStep
  );
#endif
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_UNIT_TEST_V01)
void 
OperatorScan::scanUnitTest_v01
(
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  struct kernelStatistics             *kernelStats,
#endif
  cl_uint                             currentTimeStep,
  cl_uint                             maxCount,
  cl::CommandQueue                    &queue,
  cl::Event                           &ndrEvt
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(this->initDataUnitTest_v01)
  {
    throw SimException("OperatorScan::scanUnitTest_v01: attemp to run unit test without setup");
  }
#endif
  
  LOG("OperatorScan::scanUnitTest_v01: performing unit test", 0);
  
  this->setUnitTestData_v01(queue, CL_TRUE, maxCount);

  this->scan_v01
  (
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
    kernelStats,
#endif
    currentTimeStep,
    this->dataHistogramV01Buffer,
    queue,
    ndrEvt
  );
  
  this->invalidateUnitTestData_v01();

#if SCAN_VERIFY_ENABLE
  this->verifyScan_v01
  (
    this->histogramBinBackets,
    this->histogramBinCount_v01,
    this->histogramBinSize_v01,
    queue
  );
#endif
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_UNIT_TEST_V00)
void
OperatorScan::setupUnitTest_v00
(
  cl_uint                             timeSlots,
  cl_uint                             histogramBinCount,
  cl_uint                             histogramBinSize,
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             *kernelStats
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(!this->initDataUnitTest_v00)
  {
    throw SimException("OperatorScan::setupUnitTest_v00: attemp to setup unit test without reset");
  }
#endif
  
  this->initDataUnitTest_v00 = false;
  this->timeSlots = timeSlots;
  this->histogramBinCount_v00 = histogramBinCount;
  this->histogramBinSize_v00 = histogramBinSize;
  size_t size = (this->timeSlots)*((this->histogramBinCount_v00)*(this->histogramBinSize_v00) + 1);
  
  CALLOC_O(dataHistogramV00, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataHistogramV00, kernelStats);

  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataHistogramV00Buffer, 
    this->dataHistogramV00SizeBytes);
    
  this->dataPastHistogramV00 = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM);
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_UNIT_TEST_V01)
void
OperatorScan::setupUnitTest_v01
(
  cl_uint                             histogramBinBackets,
  cl_uint                             histogramBinCount,
  cl_uint                             histogramBinSize,
  cl::Context                         &context,
  cl::Device                          &device,
  cl::CommandQueue                    &queue,
  cl_bool                             block,
  struct kernelStatistics             *kernelStats
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(!this->initDataUnitTest_v01)
  {
    throw SimException("OperatorScan::setupUnitTest_v01: attemp to setup unit test without reset");
  }
#endif
  
  this->initDataUnitTest_v01 = false;

  this->histogramBinBackets = histogramBinBackets;
  this->histogramBinCount_v01 = histogramBinCount;
  this->histogramBinSize_v01 = histogramBinSize;
  size_t size = 
    (this->histogramBinBackets)*(this->histogramBinCount_v01)*(this->histogramBinSize_v01);
  
  CALLOC_O(dataHistogramV01, cl_uint, size);
  REGISTER_MEMORY_O(device, KERNEL_ALL, MEM_GLOBAL, dataHistogramV01, kernelStats);

  CREATE_BUFFER_O(context, CL_MEM_READ_WRITE, this->dataHistogramV01Buffer, 
    this->dataHistogramV01SizeBytes);
    
  this->dataPastHistogramV01 = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_V00 && SCAN_VERIFY_ENABLE)
void 
OperatorScan::verifyScan_v00
(
#if !(SCAN_ENABLE_UNIT_TEST_V00)
  void              *objectToCall,
  cl_uint           (*getPreviousItem)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
  cl_uint           (*getResentItem)(cl::CommandQueue &, cl_uint, cl_uint, cl_uint, void*),
#endif
  cl_uint           timeSlots,
  cl_uint           histogramBinCount,
  cl_uint           histogramBinSize,
  cl::CommandQueue  &queue,
  cl_uint           timeSlot
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

#if (SCAN_ENABLE_UNIT_TEST_V00)
  if(this->initDataUnitTest_v00)
  {
    throw SimException("OperatorScan::verifyScan_v00: attemp to verify "
      "unit test without unit test setup");
  }
  
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM);
#endif

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
#if (SCAN_ENABLE_UNIT_TEST_V00)
      element = this->dataPastHistogramV00[offset + bufferID];
#else
      element = getPreviousItem(queue, timeSlot, binID, bufferID, objectToCall);
#endif
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
#if (SCAN_ENABLE_UNIT_TEST_V00)
      element = this->dataHistogramV00[offset];
#else
      element = getResentItem(queue, timeSlot, binID, bufferID, objectToCall);
#endif
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
    std::stringstream ss;
    ss << "verifyScan_v00: failed to match verification data " 
      << errorEventTotals << " times in time slot " << timeSlot << std::endl;
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/
#endif



#if (SCAN_ENABLE_V01 && SCAN_VERIFY_ENABLE)
void 
OperatorScan::verifyScan_v01
(
#if !(SCAN_ENABLE_UNIT_TEST_V01)
  cl_uint           *dataHistogram,
  cl_uint           *dataPastHistogram,
#endif
  cl_uint           histogramBinBackets,
  cl_uint           histogramBinCount,
  cl_uint           histogramBinSize,
  cl::CommandQueue  &queue
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

#if (SCAN_ENABLE_UNIT_TEST_V01)
  if(this->initDataUnitTest_v01)
  {
    throw SimException("OperatorScan::verifyScan_v01: attemp to verify "
      "unit test without unit test setup");
  }
  
  this->getData(queue, CL_TRUE, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);
#endif

  /* allocate memory for histogram */
  cl_uint *temp = (cl_uint *)calloc((histogramBinSize*histogramBinCount + 1), 
    sizeof(cl_uint));
  
  /*Compute histogram*/
  for(cl_uint j = 0; j < histogramBinSize; j++)
  {
    for(cl_uint b = 0; b < histogramBinCount; b++)
    {
      cl_uint sum = 0;
      
      for(cl_uint w = 0; w < histogramBinBackets; w++)
      {
        cl_uint p = 
          /*WG offset*/
          w*(histogramBinSize*histogramBinCount) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*histogramBinSize;
#if (SCAN_ENABLE_UNIT_TEST_V01)
        sum += this->dataPastHistogramV01[p];
#else
        sum += dataPastHistogram[p];
#endif
      }
      temp[j + b*histogramBinSize] = sum;
    }
  }
  
  cl_uint runningSum = temp[0];
  temp[0] = 0;
  
  /*Compute offsets*/
  for(cl_uint j = 1; j < (histogramBinCount*histogramBinSize + 1); j++)
  {
    cl_uint d = temp[j];
    temp[j] = runningSum;
    runningSum += d;
  }
  
  unsigned long long error_count_histogram_out = 0;

  for(cl_uint j = 1; j < (histogramBinCount*histogramBinSize + 1); j++)
  {
#if (SCAN_ENABLE_UNIT_TEST_V01)
    error_count_histogram_out += (this->dataHistogramV01[j] != temp[j]);
#else
    error_count_histogram_out += (dataHistogram[j] != temp[j]);
#endif
    /*
    std::cout << j << "->(" << dataHistogram[j] << "," << temp[j] << "), ";
    */
  }
  
  if(temp)
    free(temp);

  /*Verify*/
  if(error_count_histogram_out)
  {
    std::stringstream ss;
    ss << "verifyScan_v01: failed to match 01 histogram " << error_count_histogram_out << " times"
    << std::endl;
    throw SimException(ss.str());
  }
}
/**************************************************************************************************/
#endif



/**************************************************************************************************/
/*  private:                                                                                      */
/**************************************************************************************************/



#if (SCAN_ENABLE_UNIT_TEST_V00)
void
OperatorScan::setUnitTestData_v00
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             timeSlot
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(this->initDataUnitTest_v00)
  {
    throw SimException("OperatorScan::setUnitTestData_v00: attemp to set "
      "unit test data without unit test setup");
  }
#endif

  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG("OperatorScan::setUnitTestData_v00: set srand seed to " << this->srandSeed, 0);
  
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



#if (SCAN_ENABLE_UNIT_TEST_V01)
void
OperatorScan::setUnitTestData_v01
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             maxCount
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
  
  if(this->initDataUnitTest_v01)
  {
    throw SimException("OperatorScan::setUnitTestData_v01: attemp to set "
      "unit test data without unit test setup");
  }
#endif

  SET_RANDOM_SEED(this->srandSeed, this->srandCounter);
  LOG("OperatorScan::setUnitTestData_v01: set srand seed to " << this->srandSeed, 0);
  
  memset(this->dataHistogramV01, 0, this->dataHistogramV01SizeBytes);
  
  for(cl_uint j = 0; j < (this->histogramBinSize_v01); j++)
  {
    for(cl_uint b = 0; b < (this->histogramBinCount_v01); b++)
    {
      cl_uint sum = 0;
      
      for(cl_uint w = 0; w < (this->histogramBinBackets); w++)
      {
        cl_uint p = 
          /*WG offset*/
          w*((this->histogramBinSize_v01)*(this->histogramBinCount_v01)) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*(this->histogramBinSize_v01);
          
        this->dataHistogramV01[p] = cl_uint(abs(maxCount*((double)rand()/((double)RAND_MAX))));
      }
    }
  }
  
  this->storeData(queue, block, OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM);
}
/**************************************************************************************************/
#endif



void
OperatorScan::storeData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_DEBUG_ENABLE)
  if(selectBitMask & OPERATOR_SCAN_VALID_DEBUG)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataScanDebugHostBuffer, 
      this->dataScanDebugHostSizeBytes, this->dataScanDebugHost);
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataScanDebugDeviceBuffer, 
      this->dataScanDebugDeviceSizeBytes, this->dataScanDebugDevice);
    
    this->dataValid |= OPERATOR_SCAN_VALID_DEBUG;
  }
#endif

#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_ERROR_TRACK_ENABLE)
  if(selectBitMask & OPERATOR_SCAN_VALID_ERROR)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataScanErrorBuffer, 
      this->dataScanErrorSizeBytes, this->dataScanError);
    
    this->dataValid |= OPERATOR_SCAN_VALID_ERROR;
  }
#endif

#if SCAN_ENABLE_UNIT_TEST_V00
  if(selectBitMask & OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataHistogramV00Buffer, 
      this->dataHistogramV00SizeBytes, this->dataHistogramV00);

    this->dataValid |= OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM;
  }
#endif

#if SCAN_ENABLE_UNIT_TEST_V01
  if(selectBitMask & OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM)
  {
    ENQUEUE_WRITE_BUFFER_O(block, queue, this->dataHistogramV01Buffer, 
      this->dataHistogramV01SizeBytes, this->dataHistogramV01);

    this->dataValid |= OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM;
  }
#endif
}
/**************************************************************************************************/



void
OperatorScan::getData
(
  cl::CommandQueue    &queue,
  cl_bool             block,
  cl_uint             selectBitMask
)
/**************************************************************************************************/
{
#if OPERATOR_SCAN_VALIDATION_ENABLE
  this->isInitialized();
#endif

#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_DEBUG_ENABLE)
  if
  (
    (selectBitMask & OPERATOR_SCAN_VALID_DEBUG) && 
    !((this->dataValid) & OPERATOR_SCAN_VALID_DEBUG)
  )
  {
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataScanDebugHostBuffer, 
      this->dataScanDebugHostSizeBytes, this->dataScanDebugHost);
      
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataScanDebugDeviceBuffer, 
      this->dataScanDebugDeviceSizeBytes, this->dataScanDebugDevice);
      
    this->dataValid |= OPERATOR_SCAN_VALID_DEBUG;
  }
#endif

#if ((SCAN_ENABLE_V00 || SCAN_ENABLE_V01) && SCAN_ERROR_TRACK_ENABLE)
  if
  (
    (selectBitMask & OPERATOR_SCAN_VALID_ERROR) && 
    !((this->dataValid) & OPERATOR_SCAN_VALID_ERROR)
  )
  {
    ENQUEUE_READ_BUFFER_O(block, queue, this->dataScanErrorBuffer, 
      this->dataScanErrorSizeBytes, this->dataScanError);
      
    this->dataValid |= OPERATOR_SCAN_VALID_ERROR;
  }
#endif

#if SCAN_ENABLE_UNIT_TEST_V00
  if
  (
    (selectBitMask & OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM) && 
    !((this->dataValid) & OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM)
  )
  {
    swap2(this->dataHistogramV00, this->dataPastHistogramV00);

    ENQUEUE_READ_BUFFER_O(block, queue, this->dataHistogramV00Buffer, 
      this->dataHistogramV00SizeBytes, this->dataHistogramV00);
      
    this->dataValid |= OPERATOR_SCAN_VALID_UT_V00_HISTOGRAM;
  }
#endif

#if SCAN_ENABLE_UNIT_TEST_V01
  if
  (
    (selectBitMask & OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM) && 
    !((this->dataValid) & OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM)
  )
  {
    swap2(this->dataHistogramV01, this->dataPastHistogramV01);

    ENQUEUE_READ_BUFFER_O(block, queue, this->dataHistogramV01Buffer, 
      this->dataHistogramV01SizeBytes, this->dataHistogramV01);
      
    this->dataValid |= OPERATOR_SCAN_VALID_UT_V01_HISTOGRAM;
  }
#endif
}
/**************************************************************************************************/



void
OperatorScan::reset
(
  bool checkForNull
)
/**************************************************************************************************/
{
  try
  {
    this->resetObject = true;
    this->dataValid = 0;
    this->dataToSimulationLogFile = NULL;
    this->dataToReportLogFile = NULL;
    
#if (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)

#if SCAN_ENABLE_V00
    this->setKernelArguments_v00 = true;
    this->blockSizeX_kernelScanHistogramV00 = 0;
    this->blockSizeY_kernelScanHistogramV00 = 0;
    this->argNumScan00 = 0;
#endif

#if (SCAN_ENABLE_UNIT_TEST_V00)
    this->initDataUnitTest_v00 = true;
    this->timeSlots = 0;
    this->histogramBinCount_v00 = 0;
    this->histogramBinSize_v00 = 0;
    this->dataHistogramV00Size = 0;
    this->dataHistogramV00SizeBytes = 0;
#endif

#if SCAN_ENABLE_V01
    this->setKernelArguments_v01 = true;
    this->blockSizeX_kernelScanHistogramV01 = 0;
    this->blockSizeY_kernelScanHistogramV01 = 0;
    this->argNumScan01 = 0;
#endif

#if (SCAN_ENABLE_UNIT_TEST_V01)
    this->initDataUnitTest_v01 = true;
    this->histogramBinBackets = 0;
    this->histogramBinCount_v01 = 0;
    this->histogramBinSize_v01 = 0;
    this->dataHistogramV01Size = 0;
    this->dataHistogramV01SizeBytes = 0;
#endif

#if (SCAN_ENABLE_UNIT_TEST_V00 || SCAN_ENABLE_UNIT_TEST_V01)
    this->srandSeed = 0;
    this->srandCounter = 0;
#endif

    if(checkForNull)
    {
#if (SCAN_ENABLE_UNIT_TEST_V00)
      if(this->dataHistogramV00)
      {
        free(this->dataHistogramV00);
        this->dataHistogramV00 = NULL;
      }
      if(this->dataPastHistogramV00)
      {
        free(this->dataPastHistogramV00);
        this->dataPastHistogramV00 = NULL;
      }
#endif

#if (SCAN_ENABLE_UNIT_TEST_V01)
      if(this->dataHistogramV01)
      {
        free(this->dataHistogramV01);
        this->dataHistogramV01 = NULL;
      }
      if(this->dataPastHistogramV01)
      {
        free(this->dataPastHistogramV01);
        this->dataPastHistogramV01 = NULL;
      }
#endif
    
#if (SCAN_DEBUG_ENABLE)
      if(this->dataScanDebugHost)
      {
        free(this->dataScanDebugHost);
        this->dataScanDebugHost = NULL;
      }
      if(this->dataScanDebugDevice)
      {
        free(this->dataScanDebugDevice);
        this->dataScanDebugDevice = NULL;
      }
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
      if(this->dataScanError)
      {
        free(this->dataScanError);
        this->dataScanError = NULL;
      }
#endif 
    }
    else
    {
#if (SCAN_ENABLE_UNIT_TEST_V00)
      this->dataHistogramV00 = NULL;
      this->dataPastHistogramV00 = NULL;
#endif
    
#if (SCAN_ENABLE_UNIT_TEST_V01)
      this->dataHistogramV01 = NULL;
      this->dataPastHistogramV01 = NULL;
#endif
    
#if (SCAN_DEBUG_ENABLE)
      this->dataScanDebugHost = NULL;
      this->dataScanDebugDevice = NULL;
#endif

#if (SCAN_ERROR_TRACK_ENABLE)
      this->dataScanError = NULL;
#endif
    }
#endif
  }
  CATCH(std::cerr, OperatorScan::reset, throw SimException("OperatorScan::reset: failed.");)
}
/**************************************************************************************************/



void
OperatorScan::isInitialized
()
/**************************************************************************************************/
{
  if(this->resetObject)
  {
    throw SimException("OperatorScan::isInitialized: the object was not initialized");
  }
}
/**************************************************************************************************/
