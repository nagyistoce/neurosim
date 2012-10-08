
/* ===============================================================================================

                                              :-)
  
  =============================================================================================== */



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Neurosim.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
void
Neurosim::allocateHostData
(
  cl::Device  &device
)
/**************************************************************************************************/
{
  cl_uint size = 0;
  
#if ((GROUP_EVENTS_ENABLE_V00) ||\
     (GROUP_EVENTS_ENABLE_V01) ||\
     (GROUP_EVENTS_ENABLE_V02) ||\
     (GROUP_EVENTS_ENABLE_V03))
  {
  size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;

  CALLOC(dataHistogramGroupEventsVerify, cl_uint, size);
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
#if UPDATE_NEURONS_ENABLE_V00
  size = UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE*UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS;
#endif
#if MAKE_EVENT_PTRS_ENABLE
  size = MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE*MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
#endif
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  size = GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
#endif
  
  /* allocate memory for destination event data and verification space*/
  CALLOC(dataGroupEventsTik, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataGroupEventsTik);
  
  CALLOC(dataGroupEventsTikVerify, cl_uint, size);
  
  /*Size alignment between kernles is not enforced since the size is querried from the data.
    Only pitck is enforced*/
#if (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03) && MAKE_EVENT_PTRS_ENABLE
#if MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS
  #error (MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
#if GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE
  #error (GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE)
#endif
#endif
#if (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03) && UPDATE_NEURONS_ENABLE_V00
#if UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS
  #error (UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
#if UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE
  #error (UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE)
#endif
#endif
#if MAKE_EVENT_PTRS_ENABLE && UPDATE_NEURONS_ENABLE_V00
#if UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS
  #error (UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS)
#endif
#if UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE
  #error (UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE)
#endif
#endif
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataDebugHostGroupEvents, cl_uint, GROUP_EVENTS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataDebugHostGroupEvents);

  /* allocate memory for debug device buffer */
  CALLOC(dataDebugDeviceGroupEvents, cl_uint, GROUP_EVENTS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataDebugDeviceGroupEvents);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  /* allocate memory for error tracking */
  CALLOC(dataErrorGroupEvents, cl_uint, GROUP_EVENTS_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataErrorGroupEvents);
#endif
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
  {
#if MAKE_EVENT_PTRS_ENABLE
  size = MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET+1;
#endif
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
#endif

  /* allocate memory for histogram and its verification*/
  CALLOC(dataHistogramGroupEventsTok, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataHistogramGroupEventsTok);
  }
#if (GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03) &&\
    MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET != \
    (GROUP_EVENTS_HISTOGRAM_BIN_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS))
  #error (MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET != \
         (GROUP_EVENTS_HISTOGRAM_BIN_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS))
#endif
#if (MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
  #error (MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
#endif
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  {
  /* allocate memory for destination event data and verification space*/
  size = GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  CALLOC(dataGroupEventsTok, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataGroupEventsTok);
  
  CALLOC(dataGroupEventsTokVerify, cl_uint, size);
  }
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
  {
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataMakeEventPtrsDebugHost, cl_uint, MAKE_EVENT_PTRS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataMakeEventPtrsDebugHost);

  /* allocate memory for debug device buffer */
  CALLOC(dataMakeEventPtrsDebugDevice, cl_uint, MAKE_EVENT_PTRS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataMakeEventPtrsDebugDevice);
#endif

#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  /* allocate memory for error tracking */
  CALLOC(dataMakeEventPtrsError, cl_uint, MAKE_EVENT_PTRS_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataMakeEventPtrsError);
#endif
  }
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
  /* allocate memory for event pointer struct*/
#if MAKE_EVENT_PTRS_ENABLE
  size = 
    /*structs*/
    MAKE_EVENT_PTRS_STRUCTS*MAKE_EVENT_PTRS_STRUCT_SIZE;
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1
  /*last small struct is for storing last limiting address*/
  size += (1+MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE);
#endif
#endif
#if UPDATE_NEURONS_ENABLE_V00
  cl_uint size1 = 
    /*structs*/
    UPDATE_NEURONS_STRUCTS_V00*UPDATE_NEURONS_STRUCT_SIZE_V00;
    
#if MAKE_EVENT_PTRS_ENABLE && UPDATE_NEURONS_ENABLE_V00
  if(size != size1)
  {
    throw SimException("Neurosim::allocateHostData: MAKE_EVENT_PTRS_ENABLE and "
      "UPDATE_NEURONS_ENABLE_V00 pointer struct size mismatch");
  }
#endif
  size = size1;
#endif

  CALLOC(dataMakeEventPtrsStruct, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataMakeEventPtrsStruct);
  
#if MAKE_EVENT_PTRS_ENABLE && UPDATE_NEURONS_ENABLE_V00
#if MAKE_EVENT_PTRS_STRUCTS != UPDATE_NEURONS_STRUCTS_V00
  #error(MAKE_EVENT_PTRS_STRUCTS != UPDATE_NEURONS_STRUCTS_V00)
#endif
#if MAKE_EVENT_PTRS_STRUCT_SIZE != UPDATE_NEURONS_STRUCT_SIZE_V00
  #error(MAKE_EVENT_PTRS_STRUCT_SIZE != UPDATE_NEURONS_STRUCT_SIZE_V00)
#endif
#if MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE != UPDATE_NEURONS_STRUCT_ELEMENT_SIZE
  #error(MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE != UPDATE_NEURONS_STRUCT_ELEMENT_SIZE)
#endif
#endif
  }
#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
#if (UPDATE_NEURONS_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataUpdateNeuronsDebugHost, cl_uint, UPDATE_NEURONS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUpdateNeuronsDebugHost);
  /* allocate memory for debug device buffer */
  CALLOC(dataUpdateNeuronsDebugDevice, cl_uint, UPDATE_NEURONS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUpdateNeuronsDebugDevice);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataUpdateNeuronsError, cl_uint, UPDATE_NEURONS_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUpdateNeuronsError);
#endif
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
  /* allocate memory for toleraces */
  CALLOC(psTolerance, CL_DATA_TYPE, UPDATE_NEURONS_TOLERANCE_CHUNKS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, psTolerance);
#endif

  /*Parameters for each neuron in the network: */
  CALLOC(modelParameters, cl_float, UPDATE_NEURONS_TOTAL_NEURONS*
    UPDATE_NEURONS_MODEL_PARAMETERS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, modelParameters);
  
  /*Variables for each neuron in the network: */
  CALLOC(modelVariables, cl_float, UPDATE_NEURONS_TOTAL_NEURONS*
    UPDATE_NEURONS_MODEL_VARIABLES);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, modelVariables);
  
  /*Variables for each neuron in the network: */
  CALLOC(constantCoefficients, cl_float, CONST_SIZE);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, constantCoefficients);
  
  /*Neuron variables and parameter for verification*/
  ne = (int *)malloc(UPDATE_NEURONS_TOTAL_NEURONS*sizeof(int));
  te_ps = (DATA_TYPE *)malloc(UPDATE_NEURONS_TOTAL_NEURONS*sizeof(DATA_TYPE));
  nrn_ps = (neuron_iz_ps *)malloc(UPDATE_NEURONS_TOTAL_NEURONS*sizeof(neuron_iz_ps));
  co = (DATA_TYPE **)malloc(4*sizeof(DATA_TYPE *));
	for(int i=0; i<4; i++)
    {co[i] = (DATA_TYPE *)malloc((UPDATE_NEURONS_PS_ORDER_LIMIT+1)*sizeof(DATA_TYPE));}
#endif
}
/**************************************************************************************************/
#endif



int
Neurosim::initializeUnsortedEvents
/**************************************************************************************************/
(
  cl_uint totalBuffers,
  cl_uint bufferSize,
  cl_uint totalNeurons,
  cl_uint maxDelay,
  cl_float minDelay,
  cl_uint histogramBinSize,
  cl_uint histogramBitMask,
  cl_uint histogramBitShift,
  cl_uint keyOffset,
  cl_float percentInhibitory,
  double   percentBufferSizeDeviation,
  cl_uint *dataUnsortedEventCounts,
  cl_uint *dataUnsortedEventTargets,
  cl_uint *dataUnsortedEventDelays,
  cl_uint *dataUnsortedEventWeights,
  cl_uint *dataHistogram
){
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < totalBuffers; b++)
  {
    cl_uint event_total = cl_uint(bufferSize*(1 - percentBufferSizeDeviation/100.0) +
        abs(bufferSize*(percentBufferSizeDeviation/100.0)*(double)rand()/((double)RAND_MAX)));

    if(event_total >= bufferSize){event_total = bufferSize;}

    dataUnsortedEventCounts[b] = event_total;
    
    cl_uint count_inhibitory = cl_uint(event_total*(percentInhibitory/100.0));
    
    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Compute pointer to event data*/
      cl_uint ptr = 
        /*Event data buffers*/
        b * 
        (bufferSize) +
        /*Current event*/
        e;

      /*Compute event data*/
      /*target neuron*/
      cl_uint target_neuron = 
        cl_uint(abs((totalNeurons-1)*((double)rand()/((double)RAND_MAX))));
      /*weight*/
      cl_float weight = 6.0f/1.4f;
      if (count_inhibitory != 0)
      {
        count_inhibitory--;
        weight = -67.0f/1.4f;
      }
      cl_float ev_weight = cl_float(weight*((double)rand()/((double)RAND_MAX)));
      /*delay*/
      cl_float ev_delay = 
        cl_float(minDelay+abs((maxDelay-minDelay)*((double)rand()/((double)RAND_MAX))));

      /*Store event*/
      dataUnsortedEventTargets[ptr] = target_neuron;
      *((cl_float *)(&dataUnsortedEventDelays[ptr])) = ev_delay;
      *((cl_float *)(&dataUnsortedEventWeights[ptr])) = ev_weight;
      
      /*Compute histogram key for target neuron based on MSBs*/
      cl_uint bin = 0;
      if(keyOffset == 0)
      {
        bin = (dataUnsortedEventTargets[ptr]>>histogramBitShift) & histogramBitMask;
      }
      else if(keyOffset == 1)
      {
        bin = (dataUnsortedEventDelays[ptr]>>histogramBitShift) & histogramBitMask;
      }
      else if(keyOffset == 2)
      {
        bin = (dataUnsortedEventWeights[ptr]>>histogramBitShift) & histogramBitMask;
      }
      
      /*Offset is based on time slot, bin, WG*/
      cl_uint offset = 
      /*WG offset*/
      b +
      /*time slot + bin with histogramBinSize as a pitch*/
      bin*histogramBinSize;

      /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
      dataHistogram[offset]++;
    }
  }
  return 1;
}
/**************************************************************************************************/



#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
int
Neurosim::initializeGrouppedEvents
/**************************************************************************************************/
(
  cl_uint enableValues,
  cl_uint totalBuffers,
  cl_uint bufferSize,
  cl_uint destinationBufferSize,
  cl_uint histogramBinSize,
  cl_uint histogramTotalBins,
  cl_uint histogramBitMask,
  cl_uint histogramBitShift,
  cl_uint histogramOutBitShift,
  cl_uint histogramOutTotalGroups,
  cl_uint histogramOutBitMask,
  size_t  histogramOutSize,
  cl_uint keyOffset,
  cl_uint *dataUnsortedEventCounts,
  cl_uint *dataUnsortedEventTargets,
  cl_uint *dataUnsortedEventDelays,
  cl_uint *dataUnsortedEventWeights,
  cl_uint *dataGroupedEvents,
  cl_uint *dataHistogram,
  cl_uint *dataHistogramOut
){
  int result = 1;

  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = (histogramTotalBins*histogramBinSize + 1);
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogram, size*sizeof(cl_uint));
  
  /*Init data for verification*/
  cl_uint total_synaptic_events = 
    dataHistogram[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
  
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = total_synaptic_events/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI) < total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = total_synaptic_event_chunks/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  if(wg_chunk_size*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);

  /*Group events for verification*/
  for(cl_uint b = 0; b < totalBuffers; b++)
  {
    cl_uint event_total = dataUnsortedEventCounts[b];
    
    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Compute pointer to event data*/
      cl_uint ptr = 
        /*Event data buffers*/
        b *
        (bufferSize) +
        /*Current event*/
        e;

      /*Access event*/
      cl_uint key = 0;
      if(keyOffset == 0)
      {
        key = dataUnsortedEventTargets[ptr];
      }
      else if(keyOffset == 1)
      {
        key = dataUnsortedEventDelays[ptr];
      }
      else if(keyOffset == 2)
      {
        key = dataUnsortedEventWeights[ptr];
      }
      
      /*Compute offset key for target neuron*/
      cl_uint bin = (key>>histogramBitShift)&histogramBitMask;
        
      /*Offset is based on time slot, bin, WG*/
      cl_uint bin_offset = 
      /*WG offset*/
      b +
      /*time slot + bin with histogramBinSize as a pitch*/
      bin*histogramBinSize;

      /*Check for offset overlap*/
      if(dataOffsetGroupEventsCopy[bin_offset] >= dataHistogram[bin_offset+1])
      {
        std::cout << "initializeGrouppedEvents: Destination event bin pointer overlaps with "
          << "next pointer for bin " << 
          bin << ", buffer " << b << std::endl;
        result = 0; 
        break;
      }
      
      /*Calculate offset in the grouped data space*/
      cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
        
      /*Compute histogram key for target neuron*/
      cl_uint hist_out_ptr = 
        /*WG offset for histogram out*/
        (dest_offset/wg_chunk_size) + 
        /*bin*/
        histogramOutTotalGroups*
        ((key>>histogramOutBitShift)&histogramOutBitMask);
      /*Verify*/
      if(hist_out_ptr > histogramOutSize)
      {
        std::cout << "initializeGrouppedEvents: Pointer to an element in output histogram is "
          << "outside of its range" << std::endl;
        result = 0;
        break;
      }
      /*Increment counter for this bin.*/
      dataHistogramOut[hist_out_ptr]++;
      
      /*Store event at its group location (grouped by bins)*/
      dataGroupedEvents[dest_offset] = key;
      if(enableValues)
      {
        dataGroupedEvents[destinationBufferSize + dest_offset] = ptr;
      }
      /*Increment ptr for next data item*/
      dataOffsetGroupEventsCopy[bin_offset]++;
    }
    if(result != 1){break;}
  }
  free(dataOffsetGroupEventsCopy);
  return result;
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V01
int 
Neurosim::initializeDataForKernelGroupEventsV01
/**************************************************************************************************/
(
  int step, 
  cl_uint keyOffset,
  double percentBufferSizeDeviation
){

  int result = 1;

  cl_uint max_offset = 0;
  
  cl_uint shiftFirstStage, shiftNextStage;
  if(step < 0)
  {
    shiftFirstStage = 0, shiftNextStage = 0;
  }
  else
  {
    shiftFirstStage = GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01 * (cl_uint)step;;
    shiftNextStage = GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01 * ((cl_uint)step + 1);
  }

#if (GROUP_EVENTS_DEBUG_ENABLE)
  memset(dataDebugHostGroupEvents, 0, dataDebugHostGroupEventsSizeBytes);
  memset(dataDebugDeviceGroupEvents, 0, dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  /*Set it with non-zero since kernel should zero it out*/
  memset(dataErrorGroupEvents, 0xF, dataErrorGroupEventsSizeBytes);
#endif
  memset((*operatorSort).dataHistogramGroupEventsTik, 0, (*operatorSort).dataHistogramGroupEventsTikSizeBytes);
  memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
  memset(dataGroupEventsTok, 0, dataGroupEventsTokSizeBytes);
  memset(dataHistogramGroupEventsTok, 0, dataHistogramGroupEventsTokSizeBytes);

  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS;
  cl_uint *dataUnsortedEventCountsTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  /* allocate memory for offset copies used in verification data generation*/
  size = GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS * GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE;
  cl_uint *dataUnsortedEventTargetsTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  cl_uint *dataUnsortedEventDelaysTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  cl_uint *dataUnsortedEventWeightsTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  size = (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + 1);
  cl_uint *dataHistogramTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  /*Initialize syn events and neuron counters*/
  result = initializeUnsortedEvents
  (
    GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
    GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
    GROUP_EVENTS_TOTAL_NEURONS,
    GROUP_EVENTS_MAX_DELAY,
    GROUP_EVENTS_MIN_DELAY,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    shiftFirstStage,
    keyOffset,
    10.0f,
    percentBufferSizeDeviation,
    dataUnsortedEventCountsTemp,
    dataUnsortedEventTargetsTemp,
    dataUnsortedEventDelaysTemp,
    dataUnsortedEventWeightsTemp,
    dataHistogramTemp
  );

  cl_uint print_bins = 0;
  
  /*Compute offsets based on histogram.*/
  
  print_bins ? std::cout 
    << "initializeDataForKernelGroupEventsV01: Number of synaptic events in bins: " 
    << std::endl, true : false;
  
  cl_uint runningSum = 0;
  cl_uint runningSize = 0;
  
  for(cl_uint j = 0; j < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE; j++)
  {
    cl_uint temp = dataHistogramTemp[j];
    dataHistogramTemp[j] = runningSum;
    runningSum += temp;

    if(j%GROUP_EVENTS_HISTOGRAM_BIN_SIZE == 0 && j != 0)
    {
      print_bins ? std::cout << runningSum-runningSize << ", ", true : false;
      runningSize = runningSum;
    }
  }
  
  dataHistogramTemp[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE] = 
    runningSum;
  
  print_bins ? std::cout << runningSum-runningSize << "] = " << runningSum << 
    std::endl, true : false;
  
  if(max_offset < runningSum){max_offset = runningSum;}

  if(max_offset > GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE)
  {
    std::cout << "initializeDataForKernelGroupEventsV01: Destination event buffer overflow. " << 
      "Need to increase GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE, which is currently " << 
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE << " above " << max_offset << std::endl;
    result = 0;
  }

  if(result != 0)
  {
    result = initializeGrouppedEvents
    (
      1,
      GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
      GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
      GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK,
      shiftFirstStage,
      shiftNextStage,
      GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT,
      (*operatorSort).dataHistogramGroupEventsTikSize,
      keyOffset,
      dataUnsortedEventCountsTemp,
      dataUnsortedEventTargetsTemp,
      dataUnsortedEventDelaysTemp,
      dataUnsortedEventWeightsTemp,
      dataGroupEventsTik,
      dataHistogramTemp,
      (*operatorSort).dataHistogramGroupEventsTik
    );
  }

  if(result != 0)
  {
  runningSum = (*operatorSort).dataHistogramGroupEventsTik[0];
  (*operatorSort).dataHistogramGroupEventsTik[0] = 0;
  
  /*Compute offsets*/
  for(cl_uint j = 1; 
    j < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE + 1); j++)
  {
    
    cl_uint d = (*operatorSort).dataHistogramGroupEventsTik[j];
    (*operatorSort).dataHistogramGroupEventsTik[j] = runningSum;
    runningSum += d;
  }
  }

  free(dataHistogramTemp);
  free(dataUnsortedEventTargetsTemp);
  free(dataUnsortedEventDelaysTemp);
  free(dataUnsortedEventWeightsTemp);
  free(dataUnsortedEventCountsTemp);
  return result;
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
int 
Neurosim::initializeDataForKernelGroupEventsV02_V03
/**************************************************************************************************/
(
  int step, 
  cl_uint keyOffset,
  double percentBufferSizeDeviation
){
  int result = 1;

  cl_uint shiftFirstStage, shiftNextStage;
  if(step < 0)
  {
    shiftFirstStage = 0, shiftNextStage = 0;
  }
  else
  {
    shiftFirstStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS * (cl_uint)step;;
    shiftNextStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT * ((cl_uint)step+1);
  }

#if (GROUP_EVENTS_DEBUG_ENABLE)
  memset(dataDebugHostGroupEvents, 0, dataDebugHostGroupEventsSizeBytes);
  memset(dataDebugDeviceGroupEvents, 0, dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  memset(dataErrorGroupEvents, 0, dataErrorGroupEventsSizeBytes);
#endif
  memset((*operatorSort).dataHistogramGroupEventsTik, 0, (*operatorSort).dataHistogramGroupEventsTikSizeBytes);
  memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
  memset((*synapticEvents).dataUnsortedEventCounts, 0, (*synapticEvents).dataUnsortedEventCountsSizeBytes);
  memset((*synapticEvents).dataUnsortedEventTargets, 0, (*synapticEvents).dataUnsortedEventTargetsSizeBytes);
  memset((*synapticEvents).dataUnsortedEventDelays, 0, (*synapticEvents).dataUnsortedEventDelaysSizeBytes);
  memset((*synapticEvents).dataUnsortedEventWeights, 0, (*synapticEvents).dataUnsortedEventWeightsSizeBytes);
  memset(dataGroupEventsTok, 0, dataGroupEventsTokSizeBytes);
  memset(dataHistogramGroupEventsTok, 0, dataHistogramGroupEventsTokSizeBytes);

  cl_uint size = (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + 1);
  cl_uint *dataHistogramTemp = (cl_uint *)calloc(size, sizeof(cl_uint));
  
  /*Initialize syn events and neuron counters*/
  result = initializeUnsortedEvents
  (
    GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
    GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
    GROUP_EVENTS_TOTAL_NEURONS,
    GROUP_EVENTS_MAX_DELAY,
    GROUP_EVENTS_MIN_DELAY,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    shiftFirstStage,
    keyOffset,
    10.0f,
    percentBufferSizeDeviation,
    (*synapticEvents).dataUnsortedEventCounts,
    (*synapticEvents).dataUnsortedEventTargets,
    (*synapticEvents).dataUnsortedEventDelays,
    (*synapticEvents).dataUnsortedEventWeights,
    dataHistogramTemp
  );

  /*Compute offsets based on histogram.*/
  cl_uint print_bins = 0;
  print_bins ? std::cout 
    << "initializeDataForKernelGroupEventsV02_V03: Number of synaptic events in bins: " 
    << std::endl, true : false;
  
  cl_uint runningSum = 0;
  cl_uint runningSize = 0;
  
  for(cl_uint j = 0; j < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE; j++)
  {
    cl_uint temp = dataHistogramTemp[j];
    dataHistogramTemp[j] = runningSum;
    runningSum += temp;

    if(j%GROUP_EVENTS_HISTOGRAM_BIN_SIZE == 0 && j != 0)
    {
      print_bins ? std::cout << runningSum-runningSize << ", ", true : false;
      runningSize = runningSum;
    }
  }
  
  dataHistogramTemp[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE] = 
    runningSum;
  
  print_bins ? std::cout << runningSum-runningSize << "] = " << runningSum << 
    std::endl, true : false;

  if(runningSum > GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE)
  {
    std::cout << "initializeDataForKernelGroupEventsV02_V03: Destination event buffer overflow. " << 
      "Need to increase GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE, which is currently " << 
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE << " above " << runningSum << std::endl;
    result = 0;
  }

  if(result != 0)
  {
    result = initializeGrouppedEvents
    (
      1,
      GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
      GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
      GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK,
      shiftFirstStage,
      shiftNextStage,
      GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT,
      (*operatorSort).dataHistogramGroupEventsTikSize,
      keyOffset,
      (*synapticEvents).dataUnsortedEventCounts,
      (*synapticEvents).dataUnsortedEventTargets,
      (*synapticEvents).dataUnsortedEventDelays,
      (*synapticEvents).dataUnsortedEventWeights,
      dataGroupEventsTik,
      dataHistogramTemp,
      (*operatorSort).dataHistogramGroupEventsTik
    );
  }

  if(result != 0)
  {
    runningSum = (*operatorSort).dataHistogramGroupEventsTik[0];
    (*operatorSort).dataHistogramGroupEventsTik[0] = 0;
    
    /*Compute offsets*/
    for(cl_uint j = 1; 
      j < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE + 1); j++)
    {
      
      cl_uint d = (*operatorSort).dataHistogramGroupEventsTik[j];
      (*operatorSort).dataHistogramGroupEventsTik[j] = runningSum;
      runningSum += d;
    }
  }

  free(dataHistogramTemp);
  return result;
}
/**************************************************************************************************/
#endif



int 
Neurosim::initializeSortedEvents
(
  cl_uint mode,
  cl_uint maxEventId,
  double eventStdDev,
  double structBufferSizeMargin,
  double gabaRatio,
  cl_uint totalEvents,
  cl_uint wfWorkSize,
  cl_uint structSize,
  cl_uint structPitch,
  cl_uint eventBufferSize,
  cl_uint *sortedEvents
)
/**************************************************************************************************/
{
#define PRINT_initializeSortedEvents  0
  
  bool enableGaba = true;
  int result = 1;
  cl_uint totalBoundaries = 0, maxBoundariesPerWf = 0, currentWindow = 0, currentId = 0;
  
#if PRINT_initializeSortedEvents
  std::cout << "initializeSortedEvents: Test conditions." 
    << "\n  mode: " << mode 
    << "\n  maxEventId: " << maxEventId 
    << "\n  totalEvents: " << totalEvents 
    << "\n  eventStdDev: " << eventStdDev 
    << "\n  structBufferSizeMargin: " << structBufferSizeMargin 
    << std::endl;
#endif

  for(cl_uint j = 1; j < totalEvents; j++)
  {
    /*Reinint window for the next row of consequitive elements*/
    if(!currentWindow && (currentId != maxEventId))
    {
      /*Adjust window according to the current events per neuron*/
      double window = ((double)(totalEvents-j)/(double)(maxEventId-currentId));
      
      if(window < 2.0){window = 2.0;}
      currentWindow = cl_uint(window*(1-eventStdDev/200.0) + (abs((window*eventStdDev/100.0)*
        ((double)rand()/((double)RAND_MAX)))));

      /*neuron ID*/
      sortedEvents[(j-1)] = currentId;
      
      currentId++;
      if(currentId >= maxEventId){currentId = maxEventId;}
      
      if(mode > 0)
      {
        /*event time*/
        cl_float time = cl_float(abs(0.1*((double)rand()/((double)RAND_MAX))));
        sortedEvents[(j-1) + eventBufferSize] = *((cl_uint *)(&time));
      }
      if(mode > 1)
      {
        /*weight*/
        double weightType = abs(100.0*((double)rand()/((double)RAND_MAX)));
        enableGaba = abs(100.0*((double)rand()/((double)RAND_MAX))) < 50.0;
        cl_float weight = 6.0f/1.4f;
        if(enableGaba && weightType < gabaRatio)
        {
          weight = -67.0f/1.4f;
        }
        weight = cl_float(weight*((double)rand()/((double)RAND_MAX)));
        sortedEvents[(j-1) + 2*eventBufferSize] = *((cl_uint *)(&weight));
      }
      totalBoundaries++;
    }
    currentWindow--;
    
    sortedEvents[j] = sortedEvents[(j-1)];
      
    if(mode > 0)
    {
      /*event time*/
      cl_float time = *((cl_float *)(&sortedEvents[(j-1) + eventBufferSize]));
      time += cl_float(abs((0.99-time)*((double)rand()/((double)RAND_MAX))));
      sortedEvents[j + eventBufferSize] = *((cl_uint *)(&time));
      if(time > 1.0)
      {
        std::cerr << "ERROR, initializeSortedEvents, event time exceeds max" << std::endl;
        result = 0;
      }
    }
    if(mode > 1)
    {
      /*weight*/
      double weightType = abs(100.0*((double)rand()/((double)RAND_MAX)));
      cl_float weight = 6.0f/1.4f;
      if(enableGaba && weightType < gabaRatio)
      {
        weight = -67.0f/1.4f;
      }
      weight = cl_float(weight*((double)rand()/((double)RAND_MAX)));
      sortedEvents[j + 2*eventBufferSize] = *((cl_uint *)(&weight));
    }

    /*Detect allocation for the next WF*/
    if(wfWorkSize > 0)
    {
      if(!(j%(wfWorkSize)))
      {
        if(maxBoundariesPerWf < totalBoundaries){maxBoundariesPerWf = totalBoundaries;}
        totalBoundaries = 0;
      }
    }
  }
  
  if(wfWorkSize > 0)
  {
    /*Detect allocation for the next WF*/
    if(maxBoundariesPerWf < totalBoundaries){maxBoundariesPerWf = totalBoundaries;}

    double structSizeLimit = ((double)(structSize/structPitch-1))*
      ((100.0-structBufferSizeMargin)/100.0);
    if((double)maxBoundariesPerWf > structSizeLimit)
    {
      std::cerr << "ERROR, initializeSortedEvents, pointer struct buffer size " 
        << " is outside of valid range: " << maxBoundariesPerWf << " > " << structSizeLimit 
        << " (" << 100.0-structBufferSizeMargin << "% of " << structSize/structPitch << std::endl;
      result = 0;
    }
  }

  return result;
}
/**************************************************************************************************/



#if MAKE_EVENT_PTRS_ENABLE
int 
Neurosim::initializeDataForKernelMakeEventPtrs
(
  int mode,
  cl_uint step
)
/**************************************************************************************************/
{
/*
  TODO
  - enhence with test cases:
    - insert key change at the inter WI, inter WF boundaries
*/

  int result = 1;

  /*Reset buffers*/
  memset(dataHistogramGroupEventsTok, 0, dataHistogramGroupEventsTokSizeBytes);
  memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
  memset(dataMakeEventPtrsStruct, 0, dataMakeEventPtrsStructSizeBytes);
  
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  for
  (
    cl_uint i = 0; 
    i < 2*MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF; 
    i++
  ){
    cl_uint gm_offset = MAKE_EVENT_PTRS_TOTAL_NEURONS*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 
      i*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
      
    dataMakeEventPtrsStruct[gm_offset]    = MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE;
    dataMakeEventPtrsStruct[gm_offset+1]  = MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE;
  }
#endif
  
  if(!mode){return result;}
  
  cl_uint totalEvents = 0, maxNeuronId = 0;
  double eventsPerNeuronDeviation = 30.0;
  double gabaRatio = 5.0*(!(step%10));
  
  if(mode == 2)
  {
    /*Init total events*/
    if(step > 10)
    {
      cl_uint structSize = cl_uint(((double)MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE)*0.1);
      double minSizeFraction = abs(1.0*((double)rand()/((double)RAND_MAX)));
      totalEvents = cl_uint(minSizeFraction*(double)structSize + 
        abs((1.0-minSizeFraction)*(double)structSize*
        ((double)rand()/((double)RAND_MAX))));
      maxNeuronId = cl_uint(abs(((double)totalEvents)*((double)rand()/((double)RAND_MAX))));
      eventsPerNeuronDeviation = abs(100.0*((double)rand()/((double)RAND_MAX)));
    }
    else if(step > 0)
    {
      double minSizeFraction = 0.8;
      totalEvents = cl_uint(minSizeFraction*(double)MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE + 
        abs((1.0-minSizeFraction)*(double)(MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE-1)*
        ((double)rand()/((double)RAND_MAX))));
      maxNeuronId = totalEvents/(1<<(step+4));
    }
  }
  else if(mode == 1)
  {
    maxNeuronId = (MAKE_EVENT_PTRS_TOTAL_NEURONS-1);
    eventsPerNeuronDeviation = 20.0;
    totalEvents = MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE;
  }

  dataHistogramGroupEventsTok[MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET] = totalEvents;
  
  /*Init sorted events*/
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  if(maxNeuronId >= MAKE_EVENT_PTRS_TOTAL_NEURONS-1)
  {
    maxNeuronId = MAKE_EVENT_PTRS_TOTAL_NEURONS-1;
  }
  
  result = initializeSortedEvents
  (
    2,
    maxNeuronId,
    eventsPerNeuronDeviation,
    10.0,
    gabaRatio,
    totalEvents,
    0,
    0,
    0,
    dataGroupEventsTikSize/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS,
    dataGroupEventsTik
  );
  
#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

  /*Compute total chunks in the grid*/
  cl_uint totalEventChunks = totalEvents/
    (MAKE_EVENT_PTRS_WF_SIZE_WI*MAKE_EVENT_PTRS_ELEMENTS_PER_WI);
  if(totalEventChunks*(MAKE_EVENT_PTRS_WF_SIZE_WI*MAKE_EVENT_PTRS_ELEMENTS_PER_WI) < totalEvents)
  {
    totalEventChunks++;
  }
  
  /*Compute total chunks per WF*/
  cl_uint chunksPerWf = totalEventChunks/(MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF);
  if(chunksPerWf*(MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF)  < totalEventChunks)
  {
    chunksPerWf++;
  }
  
  cl_uint wfWorkSize = chunksPerWf*(MAKE_EVENT_PTRS_WF_SIZE_WI*MAKE_EVENT_PTRS_ELEMENTS_PER_WI);

  result = initializeSortedEvents
  (
    0,
    maxNeuronId,
    eventsPerNeuronDeviation,
    5.0,
    gabaRatio,
    totalEvents,
    wfWorkSize,
    MAKE_EVENT_PTRS_STRUCT_SIZE,
    MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE,
    dataGroupEventsTikSize/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS,
    dataGroupEventsTik
  );
#endif
  
  return result;
}
/**************************************************************************************************/
#endif



int 
Neurosim::initializeEventPointers
(
  bool    verify,
  cl_uint totalSortedEvents,
  cl_uint totalPointers,
  cl_uint pointersPitch,
  cl_uint *sortedEvents,
  cl_uint *pointersToEvents
)
/**************************************************************************************************/
{
  int result = 1;
  
  pointersToEvents[sortedEvents[0]*pointersPitch] = 0;

  cl_uint count = 1, j = 1;
  for(j = 1; j < totalSortedEvents; j++)
  {
    /*Detect boundary*/
    if(sortedEvents[j] > sortedEvents[j-1])
    {
      if(sortedEvents[j] >= totalPointers)
      {
        std::cerr << "ERROR, initializeEventPointers, pointer data structure overflow: " 
          << sortedEvents[j] << " >= " << totalPointers << std::endl;
        result = 0;
        break;
      }
      else
      {
        pointersToEvents[sortedEvents[j]*pointersPitch] = j;
        pointersToEvents[sortedEvents[j-1]*pointersPitch + 1] = count;
        count = 0;
      }
    }
    else if(verify)
    {
      if(sortedEvents[j] < sortedEvents[j-1])
      {
        std::cerr << "ERROR, initializeEventPointers, detected violation of key sort order " 
          << sortedEvents[j] << " < " 
          << sortedEvents[j-1] << std::endl;
        result = 0;
        break;
      }
      if(*((cl_float *)(&sortedEvents[j + totalSortedEvents])) < 
        *((cl_float *)(&sortedEvents[j-1 + totalSortedEvents])))
      {
        std::cerr << "ERROR, initializeEventPointers, detected violation of value sort order " 
          << *((cl_float *)(&sortedEvents[j + totalSortedEvents])) << " < " 
          << *((cl_float *)(&sortedEvents[j-1 + totalSortedEvents])) << std::endl;
        result = 0;
        break;
      }
    }
    
    count++;
  }
  
  if(totalSortedEvents)
  {
    pointersToEvents[sortedEvents[j-1]*pointersPitch + 1] = count;
  }
  
  return result;
}
/**************************************************************************************************/



#if UPDATE_NEURONS_ENABLE_V00
int
Neurosim::psInit
(
  cl_uint     totalNeurons,
  int         injectCurrentUntilStep,
  const char  *neuronVariablesSampleFile
)
/**************************************************************************************************/
{
  /*
  const double tols[16] = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,
                           1e-11,1e-12,1e-13,1e-14,1e-15,1e-16};
  */                         
  /*
  const double dt_vals[16] = {(double)1/4,(double)1/6,(double)1/8,(double)1/10,
                              (double)1/20,(double)1/40,(double)1/60,
                              (double)1/80,(double)1/100,(double)1/200,
                              (double)1/400,(double)1/600,(double)1/800,
                              (double)1/1000,(double)1/2000};
  */
  double
    dt = UPDATE_NEURONS_DT,
    C=UPDATE_NEURONS_C,
    vr=-65,
    vt=-50,
    k=1.3,
    a=UPDATE_NEURONS_a,
    b=-9.5,
    v_reset=-85,
    u_step=0, 
    v_peak=48,
    tau_ampa = UPDATE_NEURONS_TAU_AMPA,
    tau_gaba = UPDATE_NEURONS_TAU_GABA,
    E_ampa = 0,
    E_gaba = -80,
    E = 1.0/C,
    iMax = 400.0;

  dt_ps = (DATA_TYPE)dt;
  steps_ps = (int)(floor((1.0/dt_ps)+0.5));
  tol_ps = (DATA_TYPE)UPDATE_NEURONS_PS_TOLERANCE;

  E_ps = (DATA_TYPE)(E*dt_ps);
  a_ps = (DATA_TYPE)(a*dt_ps);

  /*Pointers to neuron variables*/
  neuron_iz_ps   *nrnp_ps, *nrnx_ps;
  nrnx_ps = nrn_ps + totalNeurons;
  
  /*Read a sample of variables from a file: */
  std::vector <std::vector <CL_DATA_TYPE> > neuronVariablesSample;
  if(neuronVariablesSampleFile != NULL)
  {
    std::ifstream infile(neuronVariablesSampleFile);
    
    if(!infile.is_open())
    {
      std::cerr << "ERROR, psInit: Failed to open neuronVariablesSampleFile: " << 
        neuronVariablesSampleFile << std::endl;
      return 0;
    }

    bool header = true;
    while(infile)
    {
      std::string s;
      if(!getline(infile, s)) break;
      if(header){header = false; continue;}

      std::istringstream ss(s);
      std::vector <CL_DATA_TYPE> record;

      while(ss)
      {
        std::string s;
        if(!getline(ss, s, ',')) break;
        record.push_back((CL_DATA_TYPE)atof(s.c_str()));
      }

      neuronVariablesSample.push_back( record );
    }
    
    if(!infile.eof())
    {
      std::cerr << "ERROR, psInit: Failed to close neuronVariablesSampleFile: " << 
        neuronVariablesSampleFile << std::endl;
      return 0;
    }
    
    infile.close();
  }
  
  /*Initialise neuron parameters and variables: */
  int i; 
  size_t sampleSize = neuronVariablesSample.size();
  for(nrnp_ps = nrn_ps, i=0; nrnp_ps < nrnx_ps; nrnp_ps++, i++)
  {
    nrnp_ps->vr      = (DATA_TYPE)vr; 
    nrnp_ps->k       = (DATA_TYPE)k; 
    nrnp_ps->u_step  = (DATA_TYPE)u_step;
    nrnp_ps->b       = (DATA_TYPE)b;
    
    /*Current initialization: random vs none */
    if(injectCurrentUntilStep > 0)
    {
      double iMaxPercent = 0.7;
      nrnp_ps->I = (DATA_TYPE)(iMax*iMaxPercent + 
        abs((iMax*(1.0-iMaxPercent))*((double)rand()/((double)RAND_MAX))));
    }
    else
    {
      nrnp_ps->I = 0;
    }
    
    /*Voltages are shifted relative to vr: */
    nrnp_ps->l       = (DATA_TYPE)(-k*(vt-vr))*(DATA_TYPE)abs(1.0 - 0.3*
      ((double)rand()/((double)RAND_MAX))); 
    nrnp_ps->v_reset = (DATA_TYPE)(v_reset-vr)*(DATA_TYPE)abs(1.0 - 0.3*
      ((double)rand()/((double)RAND_MAX)));
    nrnp_ps->v_peak  = (DATA_TYPE)(v_peak-vr)*(DATA_TYPE)abs(1.0 - 0.3*
      ((double)rand()/((double)RAND_MAX))); 
    nrnp_ps->E_ampa  = (DATA_TYPE)(E_ampa-vr); 
    nrnp_ps->E_gaba  = (DATA_TYPE)(E_gaba-vr);
    
    /*Scale time/rate constants such that dt=1 in the equations: */ 
    nrnp_ps->E = E_ps;
    nrnp_ps->a = a_ps;
    
    /*Variables: */
    if(neuronVariablesSampleFile != NULL)
    {
      std::vector <CL_DATA_TYPE> record = neuronVariablesSample[i%sampleSize];
      nrnp_ps->v      = record[0];
      nrnp_ps->u      = record[1];
      nrnp_ps->g_ampa = record[2];
      nrnp_ps->g_gaba = record[3];
    }
    else
    {
      nrnp_ps->v      = 0.0f;
      nrnp_ps->u      = 0.0f;
      nrnp_ps->g_ampa = 0.0f;
      nrnp_ps->g_gaba = 0.0f;
    }
      
    nrnp_ps->n_in   = 0;
  }
  
	/*Scale time constants to time step size*/
	tau_ampa_ps = ((DATA_TYPE)tau_ampa)/dt_ps;
	tau_gaba_ps = ((DATA_TYPE)tau_gaba)/dt_ps;
	co_g_ampa_ps = -1.0f/tau_ampa_ps;
  co_g_gaba_ps = -1.0f/tau_gaba_ps;

  for(int p = 1; p < UPDATE_NEURONS_PS_ORDER_LIMIT; p++)
  {
  	co[0][p] = E_ps/((DATA_TYPE)(p+1));/*assumes all E are the same*/
  	co[1][p] = a_ps/((DATA_TYPE)(p+1));/*assumes all a are the same*/
  	co[2][p] = -1.0f/(tau_ampa_ps*(DATA_TYPE)(p+1));
  	co[3][p] = -1.0f/(tau_gaba_ps*(DATA_TYPE)(p+1)); 
	}
  
  memset(te_ps, 0, UPDATE_NEURONS_TOTAL_NEURONS*sizeof(DATA_TYPE));
  
  /*Initialize delay time counter: */
  for(unsigned int i=0; i < UPDATE_NEURONS_TOTAL_NEURONS; i++){ne[i] = -1;}
  
  return 1;

}
/**************************************************************************************************/
#endif



#if UPDATE_NEURONS_ENABLE_V00
int 
Neurosim::initializeDataForKernelUpdateNeurons
(
  bool          resetEvents,
  bool          resetParameters,
  bool          resetVariables,
  double        gabaRatio,
  const char    *neuronVariablesSampleFile
)
/**************************************************************************************************/
{
  int result = 1;
  
  SET_RANDOM_SEED(srandSeed, srandCounter);
  LOG_SIM("initializeDataForKernelUpdateNeurons: set srand seed to " << srandSeed);

  if(resetVariables || resetParameters)
  {
    result = psInit(UPDATE_NEURONS_TOTAL_NEURONS, UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP, 
      neuronVariablesSampleFile);
    if(result != 1){return result;}
  }

  if(resetEvents)
  {
    memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
    result = initializeSortedEvents
    (
      2,
      (UPDATE_NEURONS_TOTAL_NEURONS-1),
      30.0,
      5.0,
      gabaRatio,
      dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
      0,
      0,
      0,
      dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
      dataGroupEventsTik
    );
    if(result != 1){return result;}
    
    memset(dataMakeEventPtrsStruct, 0, dataMakeEventPtrsStructSizeBytes);
    
    result = initializeEventPointers
    (
      true,
      dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
      UPDATE_NEURONS_TOTAL_NEURONS,
      UPDATE_NEURONS_STRUCT_ELEMENT_SIZE,
      dataGroupEventsTik,
      dataMakeEventPtrsStruct
    );
    if(result != 1){return result;}
  }

  if(resetVariables)
  {
    memset(modelVariables, 0, modelVariablesSizeBytes);
    for(cl_uint i=0; i<UPDATE_NEURONS_TOTAL_NEURONS; i++)
    {
      modelVariables[i]                                 = nrn_ps[i].v;
      modelVariables[UPDATE_NEURONS_TOTAL_NEURONS+i]    = nrn_ps[i].u;
      modelVariables[2*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].g_ampa;
      modelVariables[3*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].g_gaba;
    }
  }
  
  if(resetParameters)
  {
    memset(modelParameters, 0, modelParametersSizeBytes);
    for(cl_uint i=0; i<UPDATE_NEURONS_TOTAL_NEURONS; i++)
    {
      modelParameters[i]                                 = nrn_ps[i].I;
      modelParameters[UPDATE_NEURONS_TOTAL_NEURONS+i]    = nrn_ps[i].k;
      modelParameters[2*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].l;
      modelParameters[3*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].b;
      modelParameters[4*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].v_reset;
      modelParameters[5*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].u_step;
      modelParameters[6*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].v_peak;
      modelParameters[7*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].E_ampa;
      modelParameters[8*UPDATE_NEURONS_TOTAL_NEURONS+i]  = nrn_ps[i].E_gaba;
    }
    
    memset(constantCoefficients, 0, constantCoefficientsSizeBytes);
    CONST_CO(constantCoefficients,0,0) = E_ps;
    CONST_CO(constantCoefficients,1,0) = a_ps;
    CONST_CO(constantCoefficients,2,0) = co_g_ampa_ps;
    CONST_CO(constantCoefficients,3,0) = co_g_gaba_ps;
    for(cl_uint i = 1; i < UPDATE_NEURONS_PS_ORDER_LIMIT; i++)
    {
      CONST_CO(constantCoefficients,0,i) = co[0][i];
      CONST_CO(constantCoefficients,1,i) = co[1][i];
      CONST_CO(constantCoefficients,2,i) = co[2][i];
      CONST_CO(constantCoefficients,3,i) = co[3][i];
    }
    
    CONST_TOL(constantCoefficients) = tol_ps;
    
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    {
      #define psToleranceValueCount     17
      #define psToleranceValuesOffset   0
      const double psToleranceValues[psToleranceValueCount] = 
        {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,1e-16,0.0};
      
      for(cl_uint i = 0; i < UPDATE_NEURONS_TOLERANCE_CHUNKS; i++)
      {
        cl_uint selection = cl_uint(abs(psToleranceValuesOffset + 
          (psToleranceValueCount-1-psToleranceValuesOffset)*((double)rand()/((double)RAND_MAX))));
        psTolerance[i] = CL_DATA_TYPE(psToleranceValues[selection]);
      }
      #undef psToleranceValueCount
      #undef psToleranceValuesOffset
    }
#endif
  }
  
  return result;
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
void
Neurosim::registerLocalMemory
(
  cl::Device  &device
)
/**************************************************************************************************/
{
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  {
  cl_uint lmCacheGroupEvents = sizeof(cl_uint)*(GROUP_EVENTS_CACHE_SIZE_WORDS); 
  lmCacheGroupEventsSizeBytes = lmCacheGroupEvents;
  REGISTER_MEMORY(GROUP_EVENTS_KERNEL_NAME, MEM_LOCAL, lmCacheGroupEvents);

  cl_uint lmlocalHistogramReference = sizeof(cl_uint)*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS; 
  lmlocalHistogramReferenceSizeBytes = lmlocalHistogramReference;
  REGISTER_MEMORY(GROUP_EVENTS_KERNEL_NAME, MEM_LOCAL, lmlocalHistogramReference);

#if (GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT)
  cl_uint lmGroupEventsHistogramOut = sizeof(cl_uint)*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT);
  lmGroupEventsHistogramOutSizeBytes = lmGroupEventsHistogramOut;
  REGISTER_MEMORY(GROUP_EVENTS_KERNEL_NAME, MEM_LOCAL, lmGroupEventsHistogramOut);
#endif
  }
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
  {
  cl_uint lmGenericMakeEventPtrs = sizeof(cl_uint)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE; 
  lmGenericMakeEventPtrsSizeBytes = lmGenericMakeEventPtrs;
  REGISTER_MEMORY(MAKE_EVENT_PTRS_KERNEL_NAME, MEM_LOCAL, lmGenericMakeEventPtrs);
  
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  cl_uint lmLastElementMakeEventPtrs = sizeof(cl_uint)*
    (MAKE_EVENT_PTRS_WG_SIZE_WF*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE); 
  lmLastElementMakeEventPtrsSizeBytes = lmLastElementMakeEventPtrs;
  REGISTER_MEMORY(MAKE_EVENT_PTRS_KERNEL_NAME, MEM_LOCAL, lmLastElementMakeEventPtrs);
  
  cl_uint lmGenericGlueEventPtrs = sizeof(cl_uint)*GLUE_EVENT_WF_LM_SHARE_SIZE; 
  lmGenericGlueEventPtrsSizeBytes = lmGenericGlueEventPtrs;
  REGISTER_MEMORY(GLUE_EVENT_PTRS_KERNEL_NAME, MEM_LOCAL, lmGenericGlueEventPtrs);
#endif
  }
#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
  {
  std::string kernelTag = std::string(UPDATE_NEURONS_KERNEL_NAME) + std::string("_V00");

  cl_uint lmSpikePackets = sizeof(cl_uint)*UPDATE_NEURONS_WG_SIZE_WF_V00*
    UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS; 
  lmSpikePacketsSizeBytes = lmSpikePackets;
  REGISTER_MEMORY(kernelTag, MEM_LOCAL, lmSpikePackets);
  
  kernelTag = std::string(UPDATE_NEURONS_SPIKED_KERNEL_NAME) + std::string("_V00");

  lmSpikePackets = sizeof(cl_uint)*UPDATE_NEURONS_WG_SIZE_WF_V00*
    UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS; 
  lmSpikePacketsSizeBytes = lmSpikePackets;
  REGISTER_MEMORY(kernelTag, MEM_LOCAL, lmSpikePackets);
  }
#endif
}
/**************************************************************************************************/
#endif



void
Neurosim::getPlatformStats()
{
  cl_int err;

  /* Plaform info */
  err = cl::Platform::get(&platforms);

  if (err != CL_SUCCESS) 
  {
    std::stringstream ss;
    ss << "Neurosim::getPlatformStats: " <<  "cl::Platform::get()" << " (" << err << ")" 
      << std::endl;
    throw SimException(ss.str());
  }

  /* Iteratate over platforms */
  LOG_SIM("Number of platforms:\t\t\t\t " << platforms.size());

  for
  (
    std::vector<cl::Platform>::iterator i = platforms.begin(); 
    i != platforms.end(); 
    ++i
  ){
    LOG_SIM("  Plaform Profile:\t\t\t\t " << (*i).getInfo<CL_PLATFORM_PROFILE>().c_str());
    LOG_SIM("  Plaform Version:\t\t\t\t " << (*i).getInfo<CL_PLATFORM_VERSION>().c_str());
    LOG_SIM("  Plaform Name:\t\t\t\t\t " << (*i).getInfo<CL_PLATFORM_NAME>().c_str());
    LOG_SIM("  Plaform Vendor:\t\t\t\t " << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str());
    if ((*i).getInfo<CL_PLATFORM_EXTENSIONS>().size() > 0) 
    {
      LOG_SIM("  Plaform Extensions:\t\t\t " << (*i).getInfo<CL_PLATFORM_EXTENSIONS>().c_str());
    }
  }
  LOG_SIM(std::endl << std:: endl);

  /* Now Iteratate over each platform and its devices */
  for 
  (
    std::vector<cl::Platform>::iterator p = platforms.begin(); 
    p != platforms.end(); 
    ++p
  ){
    LOG_SIM("  Plaform Name:\t\t\t\t\t " << (*p).getInfo<CL_PLATFORM_NAME>().c_str());
       
    std::vector<cl::Device> devices;
    (*p).getDevices(CL_DEVICE_TYPE_ALL, &devices);
    
    LOG_SIM("Number of devices:\t\t\t\t " << devices.size());
  
    for (std::vector<cl::Device>::iterator i = devices.begin(); i != devices.end(); ++i) 
    {
      LOG_SIM("  Device Type:\t\t\t\t\t ");
      cl_device_type dtype = (*i).getInfo<CL_DEVICE_TYPE>();
      switch (dtype) 
      {
        case CL_DEVICE_TYPE_ACCELERATOR:
          LOG_SIM("CL_DEVICE_TYPE_ACCRLERATOR");
        break;
        case CL_DEVICE_TYPE_CPU:
          LOG_SIM("CL_DEVICE_TYPE_CPU");
        break;
        case CL_DEVICE_TYPE_DEFAULT:
          LOG_SIM("CL_DEVICE_TYPE_DEFAULT");
        break;
        case CL_DEVICE_TYPE_GPU:
          LOG_SIM("CL_DEVICE_TYPE_GPU");
        break;
      }

      LOG_SIM("  Device ID:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_VENDOR_ID>());
      LOG_SIM("  Max compute units:\t\t\t\t " << (*i).getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>());
      LOG_SIM("  Max work items dimensions:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>());

      
      cl::detail::param_traits<cl::detail::cl_device_info,CL_DEVICE_MAX_WORK_ITEM_SIZES>::
        param_type witems = (*i).getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
      for (cl_uint x = 0; x < (*i).getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>(); x++) 
      {
        LOG_SIM("    Max work items[" << x << "]:\t\t\t\t " << witems[x]);
      }

      LOG_SIM("  Max work group size:\t\t\t\t " << (*i).getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
      LOG_SIM("  Preferred vector width char:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR>());
      LOG_SIM("  Preferred vector width short:\t\t " 
        << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT>());
      LOG_SIM("  Preferred vector width int:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT>());
      LOG_SIM("  Preferred vector width long:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG>());
      LOG_SIM("  Preferred vector width float:\t\t " 
        << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT>());
      LOG_SIM("  Preferred vector width double:\t\t " 
        << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE>());
      LOG_SIM("  Max clock frequency:\t\t\t\t " << (*i).getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() 
        << "Mhz");
      LOG_SIM("  Address bits:\t\t\t\t " << (*i).getInfo<CL_DEVICE_ADDRESS_BITS>());
      LOG_SIM("  Max memeory allocation:\t\t\t " << (*i).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>());
      LOG_SIM("  Image support:\t\t\t\t " 
        << ((*i).getInfo<CL_DEVICE_IMAGE_SUPPORT>() ? "Yes" : "No"));
      
      if ((*i).getInfo<CL_DEVICE_IMAGE_SUPPORT>()) 
      {
        LOG_SIM("  Max number of images read arguments:\t\t " 
          << (*i).getInfo<CL_DEVICE_MAX_READ_IMAGE_ARGS>());
        LOG_SIM("  Max number of images write arguments:\t " 
          << (*i).getInfo<CL_DEVICE_MAX_WRITE_IMAGE_ARGS>());
        LOG_SIM("  Max image 2D width:\t\t\t\t " 
          << (*i).getInfo<CL_DEVICE_IMAGE2D_MAX_WIDTH>());
        LOG_SIM("  Max image 2D height:\t\t\t\t " 
          << (*i).getInfo<CL_DEVICE_IMAGE2D_MAX_HEIGHT>());
        LOG_SIM("  Max image 3D width:\t\t\t\t " 
          << (*i).getInfo<CL_DEVICE_IMAGE3D_MAX_WIDTH>());
        LOG_SIM("  Max image 3D height:\t\t\t\t " 
          << (*i).getInfo<CL_DEVICE_IMAGE3D_MAX_HEIGHT>());
        LOG_SIM("  Max image 3D depth:\t\t\t\t " 
          << (*i).getInfo<CL_DEVICE_IMAGE3D_MAX_DEPTH>());
        LOG_SIM("  Max samplers within kernel:\t\t\t " 
          << (*i).getInfo<CL_DEVICE_MAX_SAMPLERS>());      
      }

      LOG_SIM("  Max size of kernel argument:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_MAX_PARAMETER_SIZE>());
      LOG_SIM("  Alignment (bits) of base address:\t\t " 
        << (*i).getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>());
      LOG_SIM("  Minimum alignment (bytes) for any datatype:\t " 
        << (*i).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>());
      LOG_SIM("  Single precision floating point capability");
      LOG_SIM("    Denorms:\t\t\t\t\t " << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() & 
        CL_FP_DENORM ? "Yes" : "No"));
      LOG_SIM("    Quiet NaNs:\t\t\t\t " << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() & 
        CL_FP_INF_NAN ? "Yes" : "No"));
      LOG_SIM("    Round to nearest even:\t\t\t " << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
        CL_FP_ROUND_TO_NEAREST ? "Yes" : "No"));
      LOG_SIM("    Round to zero:\t\t\t\t " << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
        CL_FP_ROUND_TO_ZERO ? "Yes" : "No"));
      LOG_SIM("    Round to +ve and infinity:\t\t\t " << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
        CL_FP_ROUND_TO_INF ? "Yes" : "No"));
      LOG_SIM("    IEEE754-2008 fused multiply-add:\t\t " 
        << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() & CL_FP_FMA ? "Yes" : "No"));
      LOG_SIM("    Correctly rounded div and sqrt:\t\t " << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
        CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT ? "Yes" : "No"));
      LOG_SIM("  Cache type:\t\t\t\t\t ");

      switch ((*i).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>()) 
      {
        case CL_NONE:
          LOG_SIM("None");
        break;
        case CL_READ_ONLY_CACHE:
          LOG_SIM("Read only");
        break;
        case CL_READ_WRITE_CACHE:
          LOG_SIM("Read/Write");
        break;
      }
      LOG_SIM("  Cache line size:\t\t\t\t " << (*i).getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>());
      LOG_SIM("  Cache size:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>());
      LOG_SIM("  Global memory size:\t\t\t\t " << (*i).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>());
      LOG_SIM("  Constant buffer size:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>());
      LOG_SIM("  Max number of constant args:\t\t\t " << (*i).getInfo<CL_DEVICE_MAX_CONSTANT_ARGS>());
      LOG_SIM("  Local memory type:\t\t\t\t ");

      switch ((*i).getInfo<CL_DEVICE_LOCAL_MEM_TYPE>()) 
      {
        case CL_LOCAL:
          LOG_SIM("Scratchpad");
        break;
        case CL_GLOBAL:
          LOG_SIM("Scratchpad");
        break;
      }
      
      LOG_SIM("  Local memory size:\t\t\t\t " << (*i).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>());
      LOG_SIM("  Profiling timer resolution:\t\t\t " 
        << (*i).getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>());
      LOG_SIM("  Device endianess:\t\t\t\t " 
        << ((*i).getInfo<CL_DEVICE_ENDIAN_LITTLE>() ? "Little" : "Big"));
      LOG_SIM("  Available:\t\t\t\t\t " 
        << ((*i).getInfo<CL_DEVICE_AVAILABLE>() ? "Yes" : "No"));
      LOG_SIM("  Compiler available:\t\t\t\t " 
        << ((*i).getInfo<CL_DEVICE_COMPILER_AVAILABLE>() ? "Yes" : "No"));
      LOG_SIM("  Execution capabilities:\t\t\t\t ");
      LOG_SIM("    Execute OpenCL kernels:\t\t\t " 
        << ((*i).getInfo<CL_DEVICE_EXECUTION_CAPABILITIES>() & CL_EXEC_KERNEL ? "Yes" : "No"));
      LOG_SIM("    Execute native function:\t\t\t " << ((*i).getInfo<CL_DEVICE_EXECUTION_CAPABILITIES>() 
        & CL_EXEC_NATIVE_KERNEL ? "Yes" : "No"));
      LOG_SIM("  Queue properties:\t\t\t\t ");
      LOG_SIM("    Out-of-Order:\t\t\t\t " << ((*i).getInfo<CL_DEVICE_QUEUE_PROPERTIES>() & 
        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE ? "Yes" : "No"));
      LOG_SIM("    Profiling :\t\t\t\t\t " << ((*i).getInfo<CL_DEVICE_QUEUE_PROPERTIES>() & 
        CL_QUEUE_PROFILING_ENABLE ? "Yes" : "No"));
      LOG_SIM("  Platform ID:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_PLATFORM>());
      LOG_SIM("  Name:\t\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_NAME>().c_str());
      LOG_SIM("  Vendor:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_VENDOR>().c_str());
      LOG_SIM("  Driver version:\t\t\t\t " << (*i).getInfo<CL_DRIVER_VERSION>().c_str());
      LOG_SIM("  Profile:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_PROFILE>().c_str());
      LOG_SIM("  Version:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_VERSION>().c_str());
      LOG_SIM("  Extensions:\t\t\t\t\t " << (*i).getInfo<CL_DEVICE_EXTENSIONS>().c_str());
    }
  }
}



void 
Neurosim::setup()
{
#if STATISTICS_ENABLE
  double startSetupTime;
  GET_TIME_NS(startSetupTime);
#endif

  cl_int err = CL_SUCCESS;
  cl_device_type dType;
  cl_uint myDeviceId;
  bool found = false;

  err = cl::Platform::get(&platforms);
  ASSERT_CL_SUCCESS(err, "Neurosim::setupCL: Platform::get() failed.")

  /*Find target platform*/
  found = false;
  std::vector<cl::Platform>::iterator i;
  for (i = platforms.begin(); i != platforms.end(); ++i) 
  {
    /*if(!strcmp((*i).getInfo<CL_PLATFORM_VENDOR>().c_str(), TARGET_PLATFORM_VENDOR))
    {
      found = true; break;
    }*/
    found = true; break;
  }
  
  if(!found)
  {
    std::stringstream ss;
    ss << "Neurosim::setup: Unable to find target platform with vendor " << 
      (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
    throw SimException(ss.str());
  }
  
  (*i).getDevices(CL_DEVICE_TYPE_ALL, &devices);

  std::vector<cl::Device>::iterator d;

  FIND_TARGET_DEVICE(devices, TARGET_DEVICE_NAME, d, found);
  
  if(!found)
  {
    std::stringstream ss;
    ss << "Neurosim::setup: Unable to find target devices " << TARGET_DEVICE_NAME 
      << " on platform from vendor "  << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
    throw SimException(ss.str());
  }
  
  LOG_SIM("Found device " << (*d).getInfo<CL_DEVICE_NAME>());
      
  dType = (*d).getInfo<CL_DEVICE_TYPE>();
  myDeviceId = (*d).getInfo<CL_DEVICE_VENDOR_ID>();

  cl_context_properties cps[3] = 
  { 
      CL_CONTEXT_PLATFORM, 
      (cl_context_properties)(*i)(),
      0 
  };

  context = cl::Context(dType, cps, NULL, NULL, &err);
  ASSERT_CL_SUCCESS(err, "Neurosim::setupCL: Context::Context() failed.")

  devices = context.getInfo<CL_CONTEXT_DEVICES>();
  ASSERT_CL_SUCCESS(err, "Neurosim::setupCL: Context::getInfo() failed.")
  
  FIND_TARGET_DEVICE(devices, TARGET_DEVICE_NAME, d, found);
  
  if(!found)
  {
    std::stringstream ss;
    ss << "Neurosim::setup: Unable to find target devices " << TARGET_DEVICE_NAME 
      << " within created context on platform from vendor "  
      << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
    throw SimException(ss.str());
  }
  
  LOG_SIM("Creating command queue for device " << (*d).getInfo<CL_DEVICE_NAME>());
  commandQueue = cl::CommandQueue(context, *d, 0, &err);
  ASSERT_CL_SUCCESS(err, "Neurosim::setupCL: CommandQueue::CommandQueue() failed.")

  /* Allocate host memory objects*/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  this->allocateHostData(*d);
#endif

    /*Initialize data objects*/
    
#if EXPAND_EVENTS_ENABLE
  synapticEvents = new Data_SynapticEvents
  (
    context,
    *d,
    commandQueue,
    CL_FALSE,
    kernelStats,
    dataToSimulationLogFile,
    dataToReportLogFile,
    EXPAND_EVENTS_TIME_SLOTS,
    EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS,
    EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE,
    EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS,
    EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS
  );
#elif GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || \
  GROUP_EVENTS_ENABLE_V03
  synapticEvents = new Data_SynapticEvents
  (
    context,
    *d,
    commandQueue,
    CL_FALSE,
    kernelStats,
    dataToSimulationLogFile,
    dataToReportLogFile,
    GROUP_EVENTS_TIME_SLOTS,
    GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
    GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE
  );
#endif

#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes);
  }
#endif

#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataDebugHostGroupEventsBuffer, 
      dataDebugHostGroupEventsSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes);
#endif
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
  {
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes);
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  {
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes);
  }
#endif

#if EXPAND_EVENTS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  /* allocate memory for spike data */
  {
#if EXPAND_EVENTS_ENABLE
    spikeEvents = new Data_SpikeEvents
    (
      context,
      *d,
      commandQueue,
      CL_FALSE,
      kernelStats,
      dataToSimulationLogFile,
      dataToReportLogFile,
      SIMULATION_STEP_SIZE,
      EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS,
      EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS,
      EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE,
      EXPAND_EVENTS_SPIKE_PACKETS,
      EXPAND_EVENTS_TOTAL_NEURONS
    );
    
    connectome = new Data_Connectome
    (
      context,
      *d,
      commandQueue,
      CL_FALSE,
      kernelStats,
      dataToSimulationLogFile,
      dataToReportLogFile,
      EXPAND_EVENTS_TOTAL_NEURONS,
      MAX_SYNAPSES_PER_NEURON,
      SYNAPSE_DEVIATION_RATIO,
      SYNAPSE_GABA_PERCENT,
      EXPAND_EVENTS_MIN_DELAY,
      EXPAND_EVENTS_MAX_DELAY
    );
#elif UPDATE_NEURONS_ENABLE_V00
    spikeEvents = new Data_SpikeEvents
    (
      context,
      *d,
      commandQueue,
      CL_FALSE,
      kernelStats,
      dataToSimulationLogFile,
      dataToReportLogFile,
      SIMULATION_STEP_SIZE,
      UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS,
      UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS,
      UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE,
      UPDATE_NEURONS_SPIKE_PACKETS_V00,
      UPDATE_NEURONS_TOTAL_NEURONS
    );
#endif
  }
#endif

#if GROUP_EVENTS_ENABLE_V00
  {
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelGroupEventsV00,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V00=1",
      blockSizeX_kernelGroupEventsV00,
      blockSizeY_kernelGroupEventsV00
    );
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 
  {
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelGroupEventsV01,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V01=1",
      blockSizeX_kernelGroupEventsV01,
      blockSizeY_kernelGroupEventsV01
    );
  }
#endif

#if GROUP_EVENTS_ENABLE_V02
  {
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelGroupEventsV02,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V02=1",
      blockSizeX_kernelGroupEventsV02,
      blockSizeY_kernelGroupEventsV02
    );
  }
#endif

#if GROUP_EVENTS_ENABLE_V03
  {
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelGroupEventsV03,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D GROUP_EVENTS_DEVICE_V03=1",
      blockSizeX_kernelGroupEventsV03,
      blockSizeY_kernelGroupEventsV03
    );
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataMakeEventPtrsStructBuffer, 
      dataMakeEventPtrsStructSizeBytes);
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE
  {
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataMakeEventPtrsDebugHostBuffer, 
      dataMakeEventPtrsDebugHostSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataMakeEventPtrsDebugDeviceBuffer, 
      dataMakeEventPtrsDebugDeviceSizeBytes);
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataMakeEventPtrsErrorBuffer, dataMakeEventPtrsErrorSizeBytes);
#endif
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelMakeEventPtrs,
      MAKE_EVENT_PTRS_KERNEL_FILE_NAME,
      MAKE_EVENT_PTRS_KERNEL_NAME,
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF),
      blockSizeX_kernelMakeEventPtrs,
      blockSizeY_kernelMakeEventPtrs
    );
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelGlueEventPtrs,
      GLUE_EVENT_PTRS_KERNEL_FILE_NAME,
      GLUE_EVENT_PTRS_KERNEL_NAME,
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF),
      blockSizeX_kernelGlueEventPtrs,
      blockSizeY_kernelGlueEventPtrs
    );
#endif
  }
#endif

#if UPDATE_NEURONS_ENABLE_V00
  {
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataUpdateNeuronsDebugHostBuffer, 
      dataUpdateNeuronsDebugHostSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataUpdateNeuronsDebugDeviceBuffer, 
      dataUpdateNeuronsDebugDeviceSizeBytes);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataUpdateNeuronsErrorBuffer, dataUpdateNeuronsErrorSizeBytes);
#endif
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    CREATE_BUFFER(CL_MEM_READ_ONLY, psToleranceBuffer, psToleranceSizeBytes);
#endif
    CREATE_BUFFER(CL_MEM_READ_ONLY, constantCoefficientsBuffer, constantCoefficientsSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_WRITE, modelVariablesBuffer, modelVariablesSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_ONLY, modelParametersBuffer, modelParametersSizeBytes);
    
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelUpdateNeuronsV00,
      UPDATE_NEURONS_KERNEL_FILE_NAME,
      UPDATE_NEURONS_KERNEL_NAME,
#if COMPILER_FLAGS_OPTIMIZE_ENABLE
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -cl-unsafe-math-optimizations "
      "-cl-finite-math-only -cl-fast-relaxed-math",
#else
      /*"-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt -cl-opt-disable",*/
      /*"-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1",*/
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1",
#endif
      blockSizeX_kernelUpdateNeuronsV00,
      blockSizeY_kernelUpdateNeuronsV00
    );
    
    createKernel
    (
#if LOG_SIMULATION
      dataToSimulationLogFile,
#endif
      context,
      *d,
      kernelUpdateSpikedNeuronsV00,
      UPDATE_NEURONS_KERNEL_FILE_NAME,
      UPDATE_NEURONS_SPIKED_KERNEL_NAME,
#if COMPILER_FLAGS_OPTIMIZE_ENABLE
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -cl-unsafe-math-optimizations "
      "-cl-finite-math-only -cl-fast-relaxed-math",
#else
      /*"-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt -cl-opt-disable",*/
      /*"-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt",*/
      "-D COMPILER=" TOSTRING(COMPILER_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt",
#endif
      blockSizeX_kernelUpdateNeuronsV00,
      blockSizeY_kernelUpdateNeuronsV00
    );
  }
#endif

    /*Initialize operator objects*/
    
#if ENABLE_OPERATOR_SCAN
  operatorScan = new Operator_Scan
  (
#if SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || \
    ENABLE_UNIT_TEST_SCAN_V00 || ENABLE_UNIT_TEST_SCAN_V01
    CL_FALSE,
    commandQueue,
#endif
    context,
    *d,
    kernelStats,
    dataToSimulationLogFile,
    dataToReportLogFile
  );
#endif
    
#if ENABLE_OPERATOR_SORT
  operatorSort = new Operator_Sort
  (
    context,
    *d,
    commandQueue,
    CL_FALSE,
    kernelStats,
    dataToSimulationLogFile,
    dataToReportLogFile,
    GROUP_EVENTS_GRID_SIZE_WG,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT
  );
#endif
    
#if ENABLE_OPERATOR_EXPAND
  operatorExpand = new Operator_Expand
  (
#if (EXPAND_EVENTS_DEBUG_ENABLE || EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    CL_FALSE,
    commandQueue,
#endif
    context,
    *d,
    kernelStats,
    dataToSimulationLogFile,
    dataToReportLogFile
  );
#endif

  /* Register local memory for stats */
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  this->registerLocalMemory(*d);
#endif
  
  /* Log memory for stats */
  this->getMemoryUsageStats();
  
#if STATISTICS_ENABLE
  GET_TIME_NS(this->setupTime);
  this->setupTime -= startSetupTime;
#endif
}



void
Neurosim::getMemoryUsageStats
()
/**************************************************************************************************/
{
  set<std::string>::iterator kernel_name;

  LOG_SIM("Memory Allocations:");
  
  for
  (
    kernel_name = kernelStats.kernelNames.begin();
    kernel_name != kernelStats.kernelNames.end(); 
    ++kernel_name
  ){
    map<std::string, size_t>::iterator m;

    LOG_SIM("  Kernel " << *kernel_name << ":");

    /*GM allocations*/
    if(*kernel_name == KERNEL_ALL)
    {
      LOG_SIM("    Global Memory (global scope):");
    }
    map<std::string, size_t> gmSizes = kernelStats.gmSizes[*kernel_name];
    cl_ulong gmAllSizeBytes = 0, gmMaxSizeBytes = 0;
    for
    (
      m = gmSizes.begin(); 
      m != gmSizes.end(); 
      ++m
    ){
      std::string name = (m->first);
      size_t size = (m->second);
      
      if(*kernel_name == KERNEL_ALL)
      {
        if(name.compare("TOTAL") != 0)
        {
          if(gmMaxSizeBytes < size){gmMaxSizeBytes = size;}
          LOG_SIM("        " << name << ": " << ((double)size)/(1024.0*1024.0) << " MB");
        }
        else
        {
          gmAllSizeBytes = size;
        }
      }
      else
      {
        if(size != 0)
        {
          std::stringstream ss;
          ss << "Neurosim::getMemoryUsageStats: detected non-zero GM entry for " << *kernel_name 
            << ": " << name << ": " << ((double)size)/(1024.0*1024.0) << " MB" << std::endl;
          throw SimException(ss.str());
        }
      }
    }

    if(*kernel_name == KERNEL_ALL)
    {
      LOG_SIM("        TOTAL: " << ((double)gmAllSizeBytes)/(1024.0*1024.0) << " MB");
      LOG_REP("Total GM:" << (((double)gmAllSizeBytes)/(1024.0*1024.0)));
      LOG_REP("Max GM:" << (((double)gmMaxSizeBytes)/(1024.0*1024.0)));
    }

    /*CM allocations*/
    if(*kernel_name == KERNEL_ALL)
    {
      LOG_SIM("    Constant Memory (global scope):");
    }
    map<std::string, size_t> cmSizes = kernelStats.cmSizes[*kernel_name];
    cl_ulong cmAllSizeBytes = 0;
    for
    (
      m = cmSizes.begin(); 
      m != cmSizes.end(); 
      ++m
    ){
      std::string name = (m->first);
      size_t size = (m->second);

      if(*kernel_name == KERNEL_ALL)
      {
        if(name.compare("TOTAL") != 0)
        {
          LOG_SIM("        " << name << ": " << ((double)size)/(1024.0) << " KB");
        }
        else
        {
          cmAllSizeBytes = size;
        }
      }
      else
      {
        if(size != 0)
        {
          std::stringstream ss;
          ss << "Neurosim::getMemoryUsageStats: detected non-zero CM entry for " << *kernel_name 
            << ": " << name << ": " << ((double)size)/(1024.0) << " KB" << std::endl;
          throw SimException(ss.str());
        }
      }
    }

    if(*kernel_name == KERNEL_ALL)
    {
      LOG_SIM("        TOTAL: " << ((double)cmAllSizeBytes)/(1024.0) << " KB");
    }
    
    /*LM allocations*/
    if(*kernel_name != KERNEL_ALL)
    {
      LOG_SIM("    Local Memory (WG scope):");
    }
    map<std::string, size_t> lmSizes = kernelStats.lmSizes[*kernel_name];
    cl_ulong lmAllSizeBytes = 0;
    for
    (
      m = lmSizes.begin(); 
      m != lmSizes.end(); 
      ++m
    ){
      std::string name = (m->first);
      size_t size = (m->second);

      if(*kernel_name != KERNEL_ALL)
      {
        if(name.compare("TOTAL") != 0)
        {
          LOG_SIM("        " << name << ": " << ((double)size)/(1024.0) << " KB");
        }
        else
        {
          lmAllSizeBytes = size;
        }
      }
      else
      {
        if(size != 0)
        {
          std::stringstream ss;
          ss << "Neurosim::getMemoryUsageStats: detected non-zero LM entry for " << *kernel_name 
            << ": " << name << ": " << ((double)size)/(1024.0) << " KB" << std::endl;
          throw SimException(ss.str());
        }
      }
    }
    
    if(*kernel_name != KERNEL_ALL)
    {
      LOG_SIM("        TOTAL: " << ((double)lmAllSizeBytes)/(1024.0) << " KB");
    }
  }
}
/**************************************************************************************************/



void 
Neurosim::run()
{
#if STATISTICS_ENABLE
  double startRunTime;
  GET_TIME_NS(startRunTime);
#endif

  cl_int status;
  cl_int eventStatus = CL_QUEUED;
  cl::Event writeEvt = NULL;

  /*
    Initialize network
  */
#if EXPAND_EVENTS_ENABLE
    (*synapticEvents).clearUnsortedEvents
    (
      commandQueue, 
      CL_FALSE,
      0x3
    );
#if PREINITIALIZE_NETWORK_STATE
    (*spikeEvents).setEvents
    (
      commandQueue,
      CL_FALSE,
      PREINITIALIZE_NETWORK_MIN_SPIKE_PERCENT,
      PREINITIALIZE_NETWORK_MAX_SPIKE_PERCENT,
      NULL
    );
#else
    (*spikeEvents).setEvents
    (
      commandQueue,
      CL_FALSE,
      INITIALIZE_SPIKES_MIN_MAX_PERCENT
    );
#endif
#endif

#if PREINITIALIZE_NETWORK_STATE && GROUP_EVENTS_ENABLE_V00 && EXPAND_EVENTS_ENABLE
    {
      (*synapticEvents).setUnsortedEvents
      (
#if CLASS_VALIDATION_ENABLE
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
#endif
        commandQueue,
        CL_TRUE,
        false,
        PREINITIALIZE_NETWORK_TIME_SLOT_DELTA,
        GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00,
        GROUP_EVENTS_HISTOGRAM_BIN_MASK,
        GROUP_EVENTS_TOTAL_NEURONS,
        PREINITIALIZE_NETWORK_TIME_SLOT_DELTA_DEVIATION,
        PREINITIALIZE_NETWORK_PERCENT_INHIBITORY,
        GROUP_EVENTS_MIN_DELAY,
        1.0
      );
    }
#endif

#if MAKE_EVENT_PTRS_ENABLE
    if(initializeDataForKernelMakeEventPtrs(0, 0) != 1)
    {
      std::cout << "Failed initializeDataForKernelMakeEventPtrs" << std::endl; 
    }
#endif

#if UPDATE_NEURONS_ENABLE_V00
    (*spikeEvents).clearEvents(commandQueue, CL_FALSE);
#if PREINITIALIZE_NETWORK_STATE
    if(initializeDataForKernelUpdateNeurons(0, 1, 1, 0, 
      PREINITIALIZE_NETWORK_NEURON_VARIABLES_SAMPLE_FILE) != 1)
#else
    if(initializeDataForKernelUpdateNeurons(0, 1, 1, 0, NULL) != 1)
#endif
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelUpdateNeurons");
    }
#if (LOG_MODEL_VARIABLES)
    dataToTraceFile << startTimeStamp;
    dataToTraceFile << LOG_MODEL_VARIABLES_FILE_HEADER << std::endl;
#endif
#endif

#if PREINITIALIZE_NETWORK_STATE && GROUP_EVENTS_ENABLE_V00 && EXPAND_EVENTS_ENABLE && \
    UPDATE_NEURONS_ENABLE_V00
    {
      int res = injectUnsortedEvents
      (
        GROUP_EVENTS_TIME_SLOTS,
        GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
        GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
        REFERENCE_EVENT_QUEUE_SIZE,
        (*synapticEvents).dataUnsortedEventCounts,
        (*synapticEvents).dataUnsortedEventTargets,
        (*synapticEvents).dataUnsortedEventWeights,
        (*synapticEvents).dataUnsortedEventDelays,
        nrn_ps
      );
      
      if(res != 1)
      {
        throw SimException("Neurosim::run: Failed injectUnsortedEvents");
      }
    }
#endif

#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, dataGroupEventsTikBuffer, 
    dataGroupEventsTikSizeBytes, dataGroupEventsTik);
  }
#endif

#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
    dataDebugHostGroupEvents);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, commandQueue, dataDebugDeviceGroupEventsBuffer, 
    dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, commandQueue, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
    dataErrorGroupEvents);
#endif
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, dataHistogramGroupEventsTokBuffer, 
    dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, dataGroupEventsTokBuffer, 
    dataGroupEventsTokSizeBytes, dataGroupEventsTok);
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, dataMakeEventPtrsStructBuffer, 
    dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
  }
#endif

#if UPDATE_NEURONS_ENABLE_V00
  {
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, psToleranceBuffer, psToleranceSizeBytes, psTolerance);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsErrorBuffer, dataUpdateNeuronsErrorSizeBytes, 
      dataUpdateNeuronsError);
#endif
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, constantCoefficientsBuffer, constantCoefficientsSizeBytes, 
    constantCoefficients);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
    modelVariables);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, modelParametersBuffer, modelParametersSizeBytes, 
    modelParameters);
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsErrorBuffer, dataMakeEventPtrsErrorSizeBytes, 
      dataMakeEventPtrsError);
#endif
#endif

  /* Set arguments to the kernels */

#if GROUP_EVENTS_ENABLE_V00
  cl_uint argNumGrupEventsV00 = 0;
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV00, dataDebugHostGroupEventsBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, dataDebugDeviceGroupEventsBuffer, argNumGrupEventsV00++);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV00, dataErrorGroupEventsBuffer, argNumGrupEventsV00++);
#endif
  SET_KERNEL_ARG(kernelGroupEventsV00, (*operatorSort).dataHistogramGroupEventsTikBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, (*synapticEvents).dataUnsortedEventCountsBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, (*synapticEvents).dataUnsortedEventDelaysBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, dataGroupEventsTikBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, (*synapticEvents).dataHistogramBuffer, argNumGrupEventsV00++);
  }
#endif

#if GROUP_EVENTS_ENABLE_V01
  cl_uint argNumGrupEventsV01 = 0;
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV01, dataDebugHostGroupEventsBuffer, argNumGrupEventsV01++);
  SET_KERNEL_ARG(kernelGroupEventsV01, dataDebugDeviceGroupEventsBuffer, argNumGrupEventsV01++);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV01, dataErrorGroupEventsBuffer, argNumGrupEventsV01++);
#endif
  }
#endif

#if GROUP_EVENTS_ENABLE_V02
  cl_uint argNumGrupEventsV02 = 0;
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV02, dataDebugHostGroupEventsBuffer, argNumGrupEventsV02++);
  SET_KERNEL_ARG(kernelGroupEventsV02, dataDebugDeviceGroupEventsBuffer, argNumGrupEventsV02++);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV02, dataErrorGroupEventsBuffer, argNumGrupEventsV02++);
#endif
  SET_KERNEL_ARG(kernelGroupEventsV02, (*synapticEvents).dataUnsortedEventTargetsBuffer, argNumGrupEventsV02++);
  SET_KERNEL_ARG(kernelGroupEventsV02, (*synapticEvents).dataUnsortedEventDelaysBuffer, argNumGrupEventsV02++);
  SET_KERNEL_ARG(kernelGroupEventsV02, (*synapticEvents).dataUnsortedEventWeightsBuffer, argNumGrupEventsV02++);
  }
#endif

#if GROUP_EVENTS_ENABLE_V03
  cl_uint argNumGrupEventsV03 = 0;
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV03, dataDebugHostGroupEventsBuffer, argNumGrupEventsV03++);
  SET_KERNEL_ARG(kernelGroupEventsV03, dataDebugDeviceGroupEventsBuffer, argNumGrupEventsV03++);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelGroupEventsV03, dataErrorGroupEventsBuffer, argNumGrupEventsV03++);
#endif
  SET_KERNEL_ARG(kernelGroupEventsV03, (*synapticEvents).dataUnsortedEventTargetsBuffer, argNumGrupEventsV03++);
  SET_KERNEL_ARG(kernelGroupEventsV03, (*synapticEvents).dataUnsortedEventDelaysBuffer, argNumGrupEventsV03++);
  SET_KERNEL_ARG(kernelGroupEventsV03, (*synapticEvents).dataUnsortedEventWeightsBuffer, argNumGrupEventsV03++);
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE
  cl_uint argNumMakeEventPtrs = 0;
  {
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelMakeEventPtrs, dataMakeEventPtrsDebugHostBuffer, argNumMakeEventPtrs++);
  SET_KERNEL_ARG(kernelMakeEventPtrs, dataMakeEventPtrsDebugDeviceBuffer, argNumMakeEventPtrs++);
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelMakeEventPtrs, dataMakeEventPtrsErrorBuffer, argNumMakeEventPtrs++);
#endif
  }
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  cl_uint argNumGlueEventPtrs = 0;
  {
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelGlueEventPtrs, dataMakeEventPtrsDebugHostBuffer, argNumGlueEventPtrs++);
  SET_KERNEL_ARG(kernelGlueEventPtrs, dataMakeEventPtrsDebugDeviceBuffer, argNumGlueEventPtrs++);
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelGlueEventPtrs, dataMakeEventPtrsErrorBuffer, argNumGlueEventPtrs++);
#endif
  SET_KERNEL_ARG(kernelGlueEventPtrs, dataMakeEventPtrsStructBuffer, argNumGlueEventPtrs++);
  }
#endif
#endif

#if UPDATE_NEURONS_ENABLE_V00
  cl_uint argNumUpdateNeuronsV00 = 0;
  {
#if (UPDATE_NEURONS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataUpdateNeuronsDebugHostBuffer, 
    argNumUpdateNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataUpdateNeuronsDebugDeviceBuffer, 
    argNumUpdateNeuronsV00++);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataUpdateNeuronsErrorBuffer, argNumUpdateNeuronsV00++);
#endif
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, psToleranceBuffer, argNumUpdateNeuronsV00++);
#endif
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, constantCoefficientsBuffer, argNumUpdateNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, modelParametersBuffer, argNumUpdateNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, modelVariablesBuffer, argNumUpdateNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, (*spikeEvents).dataSpikePacketCountsBuffer, 
    argNumUpdateNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, (*spikeEvents).dataSpikePacketsBuffer, 
    argNumUpdateNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataMakeEventPtrsStructBuffer, argNumUpdateNeuronsV00++);
  }
  
  cl_uint argNumUpdateSpikedNeuronsV00 = 0;
  {
#if (UPDATE_NEURONS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataUpdateNeuronsDebugHostBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataUpdateNeuronsDebugDeviceBuffer, 
    argNumUpdateSpikedNeuronsV00++);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataUpdateNeuronsErrorBuffer, 
    argNumUpdateSpikedNeuronsV00++);
#endif
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, psToleranceBuffer, 
    argNumUpdateSpikedNeuronsV00++);
#endif
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, constantCoefficientsBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, modelParametersBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, modelVariablesBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, (*spikeEvents).dataSpikePacketCountsBuffer, 
    argNumUpdateSpikedNeuronsV00++);  
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, (*spikeEvents).dataSpikePacketsBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataMakeEventPtrsStructBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  }
#endif

  /* 
  * Enqueue a kernel run call.
  */
  
  cl::Event ndrEvt;

#if GROUP_EVENTS_ENABLE_V00
  cl::NDRange globalThreadsGroupEventsV00(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG);
  cl::NDRange localThreadsGroupEventsV00(blockSizeX_kernelGroupEventsV00, 
    blockSizeY_kernelGroupEventsV00);
#endif

#if GROUP_EVENTS_ENABLE_V01
  cl::NDRange globalThreadsGroupEventsV01(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG);
  cl::NDRange localThreadsGroupEventsV01(blockSizeX_kernelGroupEventsV01, 
    blockSizeY_kernelGroupEventsV01);
#endif

#if GROUP_EVENTS_ENABLE_V02
  cl::NDRange globalThreadsGroupEventsV02(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG);
  cl::NDRange localThreadsGroupEventsV02(blockSizeX_kernelGroupEventsV02, 
    blockSizeY_kernelGroupEventsV02);
#endif

#if GROUP_EVENTS_ENABLE_V03
  cl::NDRange globalThreadsGroupEventsV03(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG);
  cl::NDRange localThreadsGroupEventsV03(blockSizeX_kernelGroupEventsV03, 
    blockSizeY_kernelGroupEventsV03);
#endif

#if MAKE_EVENT_PTRS_ENABLE
  cl::NDRange globalThreadsMakeEventPtrs(MAKE_EVENT_PTRS_WG_SIZE_WI*MAKE_EVENT_PTRS_GRID_SIZE_WG);
  cl::NDRange localThreadsMakeEventPtrs(blockSizeX_kernelMakeEventPtrs, 
    blockSizeY_kernelMakeEventPtrs); 
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  cl::NDRange globalThreadsGlueEventPtrs(GLUE_EVENT_PTRS_WG_SIZE_WI*GLUE_EVENT_PTRS_GRID_SIZE_WG);
  cl::NDRange localThreadsGlueEventPtrs(blockSizeX_kernelGlueEventPtrs, 
    blockSizeY_kernelGlueEventPtrs);
#endif
#endif

#if UPDATE_NEURONS_ENABLE_V00
  cl::NDRange globalThreadsUpdateNeuronsV00(UPDATE_NEURONS_WG_SIZE_WI_V00*
    UPDATE_NEURONS_GRID_SIZE_WG_V00);
  cl::NDRange localThreadsUpdateNeuronsV00(blockSizeX_kernelUpdateNeuronsV00, 
    blockSizeY_kernelUpdateNeuronsV00);
#endif

#if PROFILING_MODE == 1 && START_PROFILING_AT_STEP > -1
  double startAppTime = 0, endAppTime = 0;
#endif

#if SIMULATION_SNAPSHOT
  bool takeSimSnapshot = false;
#endif

  /*Iterate through steps*/
  for(currentTimeStep = 0; currentTimeStep < SIMULATION_TIME_STEPS; currentTimeStep++)
  {
    currentTimeSlot = currentTimeStep%EVENT_TIME_SLOTS;
    
#if SIMULATION_SNAPSHOT
    if(currentTimeStep == SIMULATION_TIME_STEPS-1){takeSimSnapshot = true;}
    else{takeSimSnapshot = false;}
#endif

#if (SIMULATION_MODE == 0 || ERROR_TRACK_ENABLE)
#if ERROR_TRACK_ENABLE
#if ERROR_TRACK_ACCESS_EVERY_STEPS > KERNEL_ENDSTEP_VERIFY_EVERY_STEPS
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
#else
    if(!((currentTimeStep+1)%KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
#endif
#endif
    {
      std::cout << "\nExecuting step " << currentTimeStep << "(" << currentTimeSlot << ") out of " 
        << SIMULATION_TIME_STEPS << std::endl;
    }
#endif
#if PROFILING_MODE == 1 && START_PROFILING_AT_STEP > -1
    if(currentTimeStep == START_PROFILING_AT_STEP)
    {
      std::cout << "\nStarted device timing at step " << currentTimeStep << std::endl;
      GET_TIME_NS(startAppTime);
    }
#endif
/**************************************************************************************************/
#if ENABLE_OPERATOR_EXPAND
    {
#if ENABLE_UNIT_TEST_EXPAND_EVENTS

    (*operatorExpand).expandUnitTest
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      currentTimeStep,
      EXPAND_EVENTS_TIME_SLOTS,
      EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT,
      EXPAND_EVENTS_TEST_MODE,
      INITIALIZE_SPIKES_MIN_MAX_PERCENT,
      SYNAPSE_GABA_PERCENT,
      EXPAND_EVENTS_MIN_DELAY,
      EXPAND_EVENTS_MAX_DELAY,
      (*spikeEvents),
      (*synapticEvents),
      (*connectome)
    );
    
#elif OVERWRITE_SPIKES_UNTILL_STEP > 0

    (*operatorExpand).expand
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      currentTimeStep,
      EXPAND_EVENTS_TIME_SLOTS,
      OVERWRITE_SPIKES_UNTILL_STEP,
      INITIALIZE_SPIKES_MIN_MAX_PERCENT,
      (*spikeEvents),
      (*synapticEvents),
      (*connectome)
    );
    
#else

    (*operatorExpand).expand
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      currentTimeStep,
      EXPAND_EVENTS_TIME_SLOTS,
      (*spikeEvents),
      (*synapticEvents),
      (*connectome)
    );
#endif

#if EXPAND_EVENTS_VERIFY_ENABLE
    (*operatorExpand).verifyExpand
    (
      commandQueue,
      false,
      currentTimeStep, 
#if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
      (EXPAND_EVENTS_SPIKE_PACKETS_PER_WF),
#else
      (EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF),
#endif
      EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT,
      EXPAND_EVENTS_HISTOGRAM_BIN_MASK,
      EXPAND_EVENTS_MAX_DELAY,
      EXPAND_EVENTS_MIN_DELAY,
      (*spikeEvents),
      (*synapticEvents),
      (*connectome)
    );
#endif
    }
#endif
/**************************************************************************************************/
#if ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00
    {
#if ENABLE_UNIT_TEST_SCAN_V00

    (*operatorScan).scanUnitTest_v00
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      currentTimeStep
    );
    
#else

#if SCAN_VERIFY_ENABLE
    (*synapticEvents).refresh(commandQueue, CL_FALSE);
#endif

    (*operatorScan).scan_v00
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      currentTimeStep,
      SAFE_GET((*synapticEvents).timeSlots),
      (*synapticEvents).dataHistogramBuffer
    );

#if DEVICE_HOST_DATA_COHERENCE
    (*synapticEvents).invalidateHistogram();
#endif

#if SCAN_VERIFY_ENABLE
    cl_uint scanTimeStep = currentTimeStep % SAFE_GET((*synapticEvents).timeSlots);
    
    (*operatorScan).verifyScan_v00
    (
      commandQueue,
      (void*)synapticEvents, 
      Data_SynapticEvents::getPreviousHistogramItem, 
      Data_SynapticEvents::getCurrentHistogramItem,
      SAFE_GET((*synapticEvents).timeSlots),
      SAFE_GET((*synapticEvents).histogramBinCount),
      SAFE_GET((*synapticEvents).histogramBinSize),
      scanTimeStep
    );
#endif
#endif
    }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00
    {
      cl_uint groupEventsTimeStep = currentTimeStep%GROUP_EVENTS_TIME_SLOTS;
/*Unit test initialization*/
#if !(EXPAND_EVENTS_ENABLE && SCAN_ENABLE_V00)
    /*Reset the data with new values every some number of steps for better represenation*/
#if !EXPAND_EVENTS_ENABLE && !SCAN_ENABLE_V00
    if(!groupEventsTimeStep)
#endif
    {
#if (GROUP_EVENTS_DEBUG_ENABLE)
    memset(dataDebugHostGroupEvents, 0, dataDebugHostGroupEventsSizeBytes);
    memset(dataDebugDeviceGroupEvents, 0, dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    memset(dataErrorGroupEvents, 0, dataErrorGroupEventsSizeBytes);
#endif
    /*Set it with non-zero since kernel should zero it out*/
    memset((*operatorSort).dataHistogramGroupEventsTik, 0xF, (*operatorSort).dataHistogramGroupEventsTikSizeBytes);
    memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
    
#if GROUP_EVENTS_TEST_MODE == -1
    double perecentInh = 5.0;
    double deltaDev = 50.0; 
    cl_uint detla = cl_uint(((GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE*deltaDev)/100)/
      GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS);
#elif (GROUP_EVENTS_TEST_MODE >= 0) && (GROUP_EVENTS_TEST_MODE <= 100)
    double perecentInh = 5.0;
    double deltaDev = GROUP_EVENTS_TEST_MODE; 
    cl_uint detla = NULL;
#endif

    (*synapticEvents).setUnsortedEvents
    (
#if CLASS_VALIDATION_ENABLE
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
#endif
      commandQueue,
      CL_TRUE,
      true,
      detla,
      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK,
      GROUP_EVENTS_TOTAL_NEURONS,
      deltaDev,
      perecentInh,
      GROUP_EVENTS_MIN_DELAY,
      0.0
    );

    (*operatorSort).storeData(commandQueue, CL_TRUE, 0x1); /*TODO: temporary, get rid of it*/

    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    }
/*End unit test initialization*/

#elif GROUP_EVENTS_VERIFY_ENABLE || SORT_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventTargetsBuffer, (*synapticEvents).dataUnsortedEventTargetsSizeBytes, 
      (*synapticEvents).dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventDelaysBuffer, (*synapticEvents).dataUnsortedEventDelaysSizeBytes, 
      (*synapticEvents).dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
      (*synapticEvents).dataUnsortedEventWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataHistogramBuffer, (*synapticEvents).dataHistogramSizeBytes, (*synapticEvents).dataHistogram);
#endif

#if SORT_VERIFY_ENABLE && \
    (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 && \
     GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE)
    if(captureUnsortedEvents((*synapticEvents).dataUnsortedEventCounts, (*synapticEvents).dataUnsortedEventTargets, 
       (*synapticEvents).dataUnsortedEventDelays, (*synapticEvents).dataUnsortedEventWeights) != 1)
    {
      throw SimException("Neurosim::run: Failed captureUnsortedEvents.");
    }
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

    SET_KERNEL_ARG(kernelGroupEventsV00, groupEventsTimeStep, argNumGrupEventsV00);
    
    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelGroupEventsV00, 
      globalThreadsGroupEventsV00, localThreadsGroupEventsV00, NULL, ndrEvt, 
      "kernelGroupEventsV00");

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::stringstream ss;
        ss << "Neurosim::run: GROUP_EVENTS_ERROR_TRACK_ENABLE, kernelGroupEventsV00: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        throw SimException(ss.str());
      }
    }
#endif

#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    (*operatorSort).getData(commandQueue, CL_TRUE, 0x1); /*TODO: temporary, get rid of it*/
    
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
      
    if(verifyKernelGroupEventsV00() != 1)
    {
      throw SimException("Neurosim::run: Failed verifyKernelGroupEventsV00");
    }
#endif
    }
#endif
/**************************************************************************************************/
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
      sortStep < (32/GROUP_EVENTS_HISTOGRAM_BIN_BITS); 
      sortStep++
    )
#endif
    {
/**************************************************************************************************/
#if ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01
    {
#if !(GROUP_EVENTS_ENABLE_V00)

    (*operatorScan).scanUnitTest_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      SCAN_HISTOGRAM_MAX_COUNT_FOR_TEST_V01
    );
    
#else

#if SCAN_VERIFY_ENABLE
    (*operatorSort).refreshHistogram(commandQueue, CL_TRUE);
#endif

    (*operatorScan).scan_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      (*operatorSort).dataHistogramGroupEventsTikBuffer
    );
    
#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if SCAN_VERIFY_ENABLE
    (*operatorScan).verifyScan_v01
    (
      commandQueue,
      (void*)operatorSort, 
      Operator_Sort::getPreviousHistogramItem,
      Operator_Sort::getCurrentHistogramItem,
      SCAN_HISTOGRAM_BIN_BACKETS,
      SCAN_HISTOGRAM_TOTAL_BINS_V01,
      SCAN_HISTOGRAM_BIN_SIZE_V01
    );
#endif
#endif
    }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02
    if(sortStep < (32/GROUP_EVENTS_HISTOGRAM_BIN_BITS)-1)
#endif
    {
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01) 
    {
#if GROUP_EVENTS_TEST_MODE == -1
    if(initializeDataForKernelGroupEventsV01(((int)sortStep)-1, 1, 0) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV01");
    }
#elif (GROUP_EVENTS_TEST_MODE >= 0) && (GROUP_EVENTS_TEST_MODE <= 100)
    if(initializeDataForKernelGroupEventsV01(((int)sortStep)-1, 1, GROUP_EVENTS_TEST_MODE) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV01");
    }
#endif

    (*operatorSort).storeData(commandQueue, CL_TRUE, 0x1); /*TODO: temporary, get rid of it*/

    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    }
/*End unit test initialization*/

#elif GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*operatorSort).dataHistogramGroupEventsTikBuffer, 
      (*operatorSort).dataHistogramGroupEventsTikSizeBytes, (*operatorSort).dataHistogramGroupEventsTik);
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

    SET_KERNEL_ARG(kernelGroupEventsV01, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV01++);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTikBuffer, argNumGrupEventsV01);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTokBuffer, argNumGrupEventsV01+1);
    SET_KERNEL_ARG(kernelGroupEventsV01, (*operatorSort).dataHistogramGroupEventsTikBuffer, argNumGrupEventsV01+2);
    SET_KERNEL_ARG(kernelGroupEventsV01, sortStep, argNumGrupEventsV01+3);
    argNumGrupEventsV01--;

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelGroupEventsV01, 
      globalThreadsGroupEventsV01, localThreadsGroupEventsV01, NULL, ndrEvt, 
      "kernelGroupEventsV01");

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::stringstream ss;
        ss << "Neurosim::run: GROUP_EVENTS_ERROR_TRACK_ENABLE, kernelGroupEventsV01: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        throw SimException(ss.str());
      }
    }
#endif

#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV01(sortStep) != 1)
    {
      throw SimException("Neurosim::run: Failed verifyKernelGroupEventsV01");
    }
#endif
    }
#endif
/**************************************************************************************************/
    }
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02
    else
#endif
    {
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V02
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01) 
    {
#if GROUP_EVENTS_TEST_MODE == -1
    if(initializeDataForKernelGroupEventsV02_V03(((int)sortStep)-1, 1, 0) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV02_V03");
    }
#elif (GROUP_EVENTS_TEST_MODE >= 0) && (GROUP_EVENTS_TEST_MODE <= 100)
    if(initializeDataForKernelGroupEventsV02_V03(((int)sortStep)-1, 1, GROUP_EVENTS_TEST_MODE) 
      != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV02_V03");
    }
#endif

    (*operatorSort).storeData(commandQueue, CL_TRUE, 0x1); /*TODO: temporary, get rid of it*/
    
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventTargetsBuffer, (*synapticEvents).dataUnsortedEventTargetsSizeBytes, 
      (*synapticEvents).dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventDelaysBuffer, (*synapticEvents).dataUnsortedEventDelaysSizeBytes, 
      (*synapticEvents).dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
      (*synapticEvents).dataUnsortedEventWeights);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    }
/*End unit test initialization*/

#elif GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventTargetsBuffer, (*synapticEvents).dataUnsortedEventTargetsSizeBytes, 
      (*synapticEvents).dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventDelaysBuffer, (*synapticEvents).dataUnsortedEventDelaysSizeBytes, 
      (*synapticEvents).dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
      (*synapticEvents).dataUnsortedEventWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*operatorSort).dataHistogramGroupEventsTikBuffer, 
      (*operatorSort).dataHistogramGroupEventsTikSizeBytes, (*operatorSort).dataHistogramGroupEventsTik);
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

    SET_KERNEL_ARG(kernelGroupEventsV02, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV02++);
    SET_KERNEL_ARG(kernelGroupEventsV02, dataGroupEventsTikBuffer, argNumGrupEventsV02);
    SET_KERNEL_ARG(kernelGroupEventsV02, dataGroupEventsTokBuffer, argNumGrupEventsV02+1);
    SET_KERNEL_ARG(kernelGroupEventsV02, (*operatorSort).dataHistogramGroupEventsTikBuffer, argNumGrupEventsV02+2);
    SET_KERNEL_ARG(kernelGroupEventsV02, sortStep, argNumGrupEventsV02+3);
    argNumGrupEventsV02--;

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelGroupEventsV02, 
      globalThreadsGroupEventsV02, localThreadsGroupEventsV02, NULL, ndrEvt, 
      "kernelGroupEventsV02");
      
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::stringstream ss;
        ss << "Neurosim::run: GROUP_EVENTS_ERROR_TRACK_ENABLE, kernelGroupEventsV02: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        throw SimException(ss.str());
      }
    }
#endif

#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV02(sortStep, GROUP_EVENTS_REPLACEMENT_KEY_OFFSET_V02) != 1)
    {
      throw SimException("Neurosim::run: Failed verifyKernelGroupEventsV02");
    }
#endif
    }
#endif
/**************************************************************************************************/
    }/*if*/
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02
    /*Swap buffers for next iteration*/
    swap2((*operatorSort).dataHistogramGroupEventsTikBuffer, dataHistogramGroupEventsTokBuffer);
    swap2(dataGroupEventsTikBuffer, dataGroupEventsTokBuffer);
#endif
    }/*for*/
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V03
    /*
      Sort positive integers (neuron IDs).
    */
    for
    (
      cl_uint sortStep = 0; 
      sortStep < GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS; 
      sortStep++
    )
#endif
    {
/**************************************************************************************************/
#if ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V01
    {
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02)

    (*operatorScan).scanUnitTest_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      SCAN_HISTOGRAM_MAX_COUNT_FOR_TEST_V01
    );
    
#else

#if SCAN_VERIFY_ENABLE
    (*operatorSort).refreshHistogram(commandQueue, CL_TRUE);
#endif

    (*operatorScan).scan_v01
    (
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE
      currentTimeStep,
#endif
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
      commandQueue,
      ndrEvt,
      (*operatorSort).dataHistogramGroupEventsTikBuffer
    );
    
#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if SCAN_VERIFY_ENABLE
    (*operatorScan).verifyScan_v01
    (
      commandQueue,
      (void*)operatorSort, 
      Operator_Sort::getPreviousHistogramItem, 
      Operator_Sort::getCurrentHistogramItem,
      SCAN_HISTOGRAM_BIN_BACKETS,
      SCAN_HISTOGRAM_TOTAL_BINS_V01,
      SCAN_HISTOGRAM_BIN_SIZE_V01
    );
#endif
#endif
    }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V03
    if(sortStep < (GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS-1))
#endif
    {
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && SCAN_ENABLE_V01 &&\
      GROUP_EVENTS_ENABLE_V02)
    {
#if GROUP_EVENTS_TEST_MODE == -1
    if(initializeDataForKernelGroupEventsV01(((int)sortStep)-1, 1, 0) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV01");
    }
#elif (GROUP_EVENTS_TEST_MODE >= 0) && (GROUP_EVENTS_TEST_MODE <= 100)
    if(initializeDataForKernelGroupEventsV01(((int)sortStep)-1, 1, GROUP_EVENTS_TEST_MODE) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV01");
    }
#endif

    (*operatorSort).storeData(commandQueue, CL_TRUE, 0x1); /*TODO: temporary, get rid of it*/

    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    }
/*End unit test initialization*/

#elif GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*operatorSort).dataHistogramGroupEventsTikBuffer, 
      (*operatorSort).dataHistogramGroupEventsTikSizeBytes, (*operatorSort).dataHistogramGroupEventsTik);
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

    SET_KERNEL_ARG(kernelGroupEventsV01, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV01++);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTikBuffer, argNumGrupEventsV01);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTokBuffer, argNumGrupEventsV01+1);
    SET_KERNEL_ARG(kernelGroupEventsV01, (*operatorSort).dataHistogramGroupEventsTikBuffer, argNumGrupEventsV01+2);
    SET_KERNEL_ARG(kernelGroupEventsV01, sortStep, argNumGrupEventsV01+3);
    argNumGrupEventsV01--;

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelGroupEventsV01, 
      globalThreadsGroupEventsV01, localThreadsGroupEventsV01, NULL, ndrEvt, 
      "kernelGroupEventsV01");
      
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::stringstream ss;
        ss << "Neurosim::run: GROUP_EVENTS_ERROR_TRACK_ENABLE, kernelGroupEventsV01: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        throw SimException(ss.str());
      }
    }
#endif

#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV01(sortStep) != 1)
    {
      throw SimException("Neurosim::run: Failed verifyKernelGroupEventsV01");
    }
#endif
    }
#endif
/**************************************************************************************************/
    }
#if GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V03
    else
#endif
    {
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V03
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && SCAN_ENABLE_V01 &&\
      GROUP_EVENTS_ENABLE_V02)
    {
#if GROUP_EVENTS_TEST_MODE == -1
    if(initializeDataForKernelGroupEventsV02_V03(((int)sortStep)-1, 1, 0) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV02_V03");
    }
#elif (GROUP_EVENTS_TEST_MODE >= 0) && (GROUP_EVENTS_TEST_MODE <= 100)
    if(initializeDataForKernelGroupEventsV02_V03(((int)sortStep)-1, 1, GROUP_EVENTS_TEST_MODE) 
      != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelGroupEventsV02_V03");
    }
#endif

    (*operatorSort).storeData(commandQueue, CL_TRUE, 0x1); /*TODO: temporary, get rid of it*/
    
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventTargetsBuffer, (*synapticEvents).dataUnsortedEventTargetsSizeBytes, 
      (*synapticEvents).dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventDelaysBuffer, (*synapticEvents).dataUnsortedEventDelaysSizeBytes, 
      (*synapticEvents).dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
      (*synapticEvents).dataUnsortedEventWeights);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    }
/*End unit test initialization*/

#elif GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventTargetsBuffer, (*synapticEvents).dataUnsortedEventTargetsSizeBytes, 
      (*synapticEvents).dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventDelaysBuffer, (*synapticEvents).dataUnsortedEventDelaysSizeBytes, 
      (*synapticEvents).dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
      (*synapticEvents).dataUnsortedEventWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*operatorSort).dataHistogramGroupEventsTikBuffer, 
      (*operatorSort).dataHistogramGroupEventsTikSizeBytes, (*operatorSort).dataHistogramGroupEventsTik);
#elif SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
        (*synapticEvents).dataUnsortedEventWeights);
    }
#endif

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

    SET_KERNEL_ARG(kernelGroupEventsV03, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV03++);
    SET_KERNEL_ARG(kernelGroupEventsV03, dataGroupEventsTikBuffer, argNumGrupEventsV03);
    SET_KERNEL_ARG(kernelGroupEventsV03, dataGroupEventsTokBuffer, argNumGrupEventsV03+1);
    SET_KERNEL_ARG(kernelGroupEventsV03, (*operatorSort).dataHistogramGroupEventsTikBuffer, argNumGrupEventsV03+2);
    SET_KERNEL_ARG(kernelGroupEventsV03, sortStep, argNumGrupEventsV03+3);
    argNumGrupEventsV03--;

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelGroupEventsV03, 
      globalThreadsGroupEventsV03, localThreadsGroupEventsV03, NULL, ndrEvt, 
      "kernelGroupEventsV03");

#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif

#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::stringstream ss;
        ss << "Neurosim::run: GROUP_EVENTS_ERROR_TRACK_ENABLE, verifyKernelGroupEventsV03: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        throw SimException(ss.str());
      }
    }
#endif

#if DEVICE_HOST_DATA_COHERENCE
    (*operatorSort).invalidateHistogram();
#endif

#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV03(sortStep) != 1)
    {
      std::stringstream ss;
      ss << "Neurosim::run: Failed verifyKernelGroupEventsV03, sort step " << sortStep 
      << std::endl; 
      throw SimException(ss.str());
    }
#endif
    }
#endif
/**************************************************************************************************/
    }/*if*/
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V03
    /*Swap buffers for next iteration*/
    swap2((*operatorSort).dataHistogramGroupEventsTikBuffer, dataHistogramGroupEventsTokBuffer);
    swap2(dataGroupEventsTikBuffer, dataGroupEventsTokBuffer);
#endif
    }/*for*/
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02 &&\
      GROUP_EVENTS_ENABLE_V03 && SCAN_ENABLE_V01)
    /*Use milder mode if UPDATE_NEURONS_ENABLE_V00*/
    int mode = (currentTimeStep == 0) ? 0 : (UPDATE_NEURONS_ENABLE_V00 ? 1 : 2);
    
    if(initializeDataForKernelMakeEventPtrs(mode, currentTimeStep-1) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelMakeEventPtrs");
    }
    
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
/*End unit test initialization*/

#elif MAKE_EVENT_PTRS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
#endif

#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsDebugHostBuffer, 
      dataMakeEventPtrsDebugHostSizeBytes, dataMakeEventPtrsDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsDebugDeviceBuffer, 
      dataMakeEventPtrsDebugDeviceSizeBytes, dataMakeEventPtrsDebugDevice);
#endif

    SET_KERNEL_ARG(kernelMakeEventPtrs, dataHistogramGroupEventsTokBuffer, argNumMakeEventPtrs);
    SET_KERNEL_ARG(kernelMakeEventPtrs, dataGroupEventsTikBuffer, argNumMakeEventPtrs+1);
    SET_KERNEL_ARG(kernelMakeEventPtrs, dataMakeEventPtrsStructBuffer, argNumMakeEventPtrs+2); /*TODO: instead of a new buffer it could be a reuse*/

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelMakeEventPtrs, 
      globalThreadsMakeEventPtrs, localThreadsMakeEventPtrs, NULL, ndrEvt, 
      "kernelMakeEventPtrs");
      
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelGlueEventPtrs, 
      globalThreadsGlueEventPtrs, localThreadsGlueEventPtrs, NULL, ndrEvt, 
      "kernelGlueEventPtrs");
#endif

#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsDebugHostBuffer, 
      dataMakeEventPtrsDebugHostSizeBytes, dataMakeEventPtrsDebugHost);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsDebugDeviceBuffer, 
      dataMakeEventPtrsDebugDeviceSizeBytes, dataMakeEventPtrsDebugDevice);
#endif

#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsErrorBuffer, dataMakeEventPtrsErrorSizeBytes, 
        dataMakeEventPtrsError);
      if(dataMakeEventPtrsError[0] != 0)
      {
        std::stringstream ss;
        ss << "Neurosim::run: MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE: Received error code "
          "from device: " << dataMakeEventPtrsError[0] << std::endl;
        throw SimException(ss.str());
      }
    }
#endif

#if MAKE_EVENT_PTRS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, 
      dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
    if(verifyKernelMakeEventPtrs() != 1)
    {
      throw SimException("Neurosim::run: Failed verifyKernelMakeEventPtrs");
    }
#endif
    }
#endif
/**************************************************************************************************/
#if (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 && \
     GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE)
#if SORT_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, 
      dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
    if(verifySortedEvents(dataGroupEventsTik, dataMakeEventPtrsStruct, 2) != 1)
    {
      throw SimException("Neurosim::run: Failed verifySortedEvents");
    }
#elif SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, 
        dataGroupEventsTikSizeBytes, dataGroupEventsTik);
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, 
        dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
    }
#endif
#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02 &&\
      GROUP_EVENTS_ENABLE_V03 && SCAN_ENABLE_V00 && SCAN_ENABLE_V01 && EXPAND_EVENTS_ENABLE)
    bool resetVarsAndParams = (currentTimeStep==0); /*TODO: fix, doesnt work with !(currentTimeStep%17);*/

    (*spikeEvents).clearEvents(commandQueue, CL_TRUE);

    if(initializeDataForKernelUpdateNeurons(!(MAKE_EVENT_PTRS_ENABLE), 
      resetVarsAndParams, resetVarsAndParams, 5.0*(!(currentTimeStep%3)), NULL) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelUpdateNeurons");
    }
    if(resetVarsAndParams)
    {
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
      ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, psToleranceBuffer, psToleranceSizeBytes, 
        psTolerance);
#endif
      ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, constantCoefficientsBuffer, constantCoefficientsSizeBytes, 
        constantCoefficients);
      ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, modelParametersBuffer, modelParametersSizeBytes, 
        modelParameters);
    }
#if !(MAKE_EVENT_PTRS_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
      dataMakeEventPtrsStruct);
# else
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
      dataMakeEventPtrsStruct);
#endif
#endif
/*END unit test initialization*/

/*Debugging*/
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsDebugHostBuffer, 
      dataUpdateNeuronsDebugHostSizeBytes, dataUpdateNeuronsDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsDebugDeviceBuffer, 
      dataUpdateNeuronsDebugDeviceSizeBytes, dataUpdateNeuronsDebugDevice);
#endif
/*END debugging*/

/*Get dataMakeEventPtrsStruct from device for profiling: pre-profiling steps in mode 1 require 
dataMakeEventPtrsStruct to have the same contents for device and host. Host reads it before
device modifies it.*/
#if (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 &&\
    GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE &&\
    EXPAND_EVENTS_ENABLE)
#if (UPDATE_NEURONS_VERIFY_ENABLE || (KERNEL_ENDSTEP_VERIFY_EVERY_STEPS > 0)) ||\
    (PROFILING_MODE == 1 && START_PROFILING_AT_STEP > -1) || SIMULATION_SNAPSHOT
#if UPDATE_NEURONS_VERIFY_ENABLE || (KERNEL_ENDSTEP_VERIFY_EVERY_STEPS > 0)
    bool readPtrs = !((currentTimeStep+1)%KERNEL_ENDSTEP_VERIFY_EVERY_STEPS) || 
      (currentTimeStep == 0);
#else
    bool readPtrs = false;
#endif
#if (PROFILING_MODE == 1 && START_PROFILING_AT_STEP > -1)
    readPtrs |= (currentTimeStep < START_PROFILING_AT_STEP);
#endif
#if OVERWRITE_SPIKES_UNTILL_STEP > 0
    readPtrs |= (currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP);
#endif
#if INJECT_CURRENT_UNTILL_STEP > 0
    readPtrs |= (currentTimeStep < INJECT_CURRENT_UNTILL_STEP);
#endif
#if SIMULATION_SNAPSHOT
    readPtrs |= (takeSimSnapshot);
#endif
    if(readPtrs)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
        dataMakeEventPtrsStruct);
    }
#endif
#endif
/*END Get dataMakeEventPtrsStruct*/

/*Set loop-dependent arguments for update kernel and equeue it*/
    SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataGroupEventsTikBuffer, argNumUpdateNeuronsV00);
    SET_KERNEL_ARG(kernelUpdateNeuronsV00, currentTimeStep, argNumUpdateNeuronsV00+1);

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelUpdateNeuronsV00, 
      globalThreadsUpdateNeuronsV00, localThreadsUpdateNeuronsV00, NULL, ndrEvt, 
      "kernelUpdateNeuronsV00");
      
#if DEVICE_HOST_DATA_COHERENCE
    (*spikeEvents).invalidateEvents();
#endif
/*END enqueue update kernel*/

/*Set loop-dependent arguments for update kernel for spiked neurons and equeue it*/
    SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataGroupEventsTikBuffer, 
      argNumUpdateSpikedNeuronsV00);
    SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, currentTimeStep, 
      argNumUpdateSpikedNeuronsV00+1);

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelUpdateSpikedNeuronsV00, 
      globalThreadsUpdateNeuronsV00, localThreadsUpdateNeuronsV00, NULL, ndrEvt, 
      "kernelUpdateSpikedNeuronsV00");
      
#if DEVICE_HOST_DATA_COHERENCE
    (*spikeEvents).invalidateEvents();
#endif
/*END enqueue update kernel for spiked neurons*/

/*Debugging: buffer exchange*/
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsDebugHostBuffer, 
      dataUpdateNeuronsDebugHostSizeBytes, dataUpdateNeuronsDebugHost);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsDebugDeviceBuffer, 
      dataUpdateNeuronsDebugDeviceSizeBytes, dataUpdateNeuronsDebugDevice);
#endif
/*END debugging*/

/*Error tracking: read error mask*/
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      bool ignoreSolverFailuresDevice = IGNORE_SOLVER_EXCEPTIONS;
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsErrorBuffer, dataUpdateNeuronsErrorSizeBytes, 
        dataUpdateNeuronsError);
      if(dataUpdateNeuronsError[0] != 0)
      {
        std::cout << "UPDATE_NEURONS_ERROR_TRACK_ENABLE: received error code from the device: ";
        PRINT_HEX(4, dataUpdateNeuronsError[0]); std::cout << std::endl;
        
        if((!ignoreSolverFailuresDevice && dataUpdateNeuronsError[0] != 0) || 
           (ignoreSolverFailuresDevice && 
           ((UPDATE_NEURONS_ERROR_NON_SOLVER_FAILURE_MASK&dataUpdateNeuronsError[0]) != 0))
        ){
          throw SimException("Neurosim::run: Failed ignoreSolverFailuresDevice condition");
        }
        else
        {
          std::cout << "UPDATE_NEURONS_ERROR_TRACK_ENABLE: ignoring solver failures\n";
          memset(dataUpdateNeuronsError, 0, UPDATE_NEURONS_ERROR_BUFFER_SIZE_WORDS*sizeof(cl_uint));
          ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataUpdateNeuronsErrorBuffer, 
            dataUpdateNeuronsErrorSizeBytes, dataUpdateNeuronsError);
        }
      }
    }
#endif
/*END error tracking*/

/*Verification*/
#if (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 &&\
    GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE &&\
    EXPAND_EVENTS_ENABLE)
/*Verification without profiling*/
#if UPDATE_NEURONS_VERIFY_ENABLE || (KERNEL_ENDSTEP_VERIFY_EVERY_STEPS > 0)
    {
    bool verify = !((currentTimeStep+1)%KERNEL_ENDSTEP_VERIFY_EVERY_STEPS) || 
      (currentTimeStep == 0);
#if OVERWRITE_SPIKES_UNTILL_STEP > 0
    verify |= (currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP);
#endif
#if INJECT_CURRENT_UNTILL_STEP > 0
    verify |= (currentTimeStep < INJECT_CURRENT_UNTILL_STEP);
#endif
    if(verify)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);
    }
#if OVERWRITE_SPIKES_UNTILL_STEP > 0
    if(currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP)
    {
      if
      (
        verifyKernelUpdateNeurons
        (
          true,
          true,
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          modelVariables,
          (*spikeEvents),
          (*connectome),
          commandQueue
        ) != 1
      )
      {
        throw SimException("Neurosim::run: Failed verifyKernelUpdateNeurons");
      }
    }
    else
#endif
    {
      if
      (
        verifyKernelUpdateNeurons
        (
          verify,
          false,
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          modelVariables,
          (*spikeEvents),
          (*connectome),
          commandQueue
        ) != 1
      )
      {
        throw SimException("Neurosim::run: Failed verifyKernelUpdateNeurons");
      }
    }
    }
/*END: Verification without profiling*/

/*Profiling mode 1: verify host-device until steps reach OVERWRITE_SPIKES_UNTILL_STEP or 
  INJECT_CURRENT_UNTILL_STEP or INJECT_CURRENT_UNTILL_STEP*/
#elif (PROFILING_MODE == 1 && START_PROFILING_AT_STEP > -1)
#if OVERWRITE_SPIKES_UNTILL_STEP > 0
    if(currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);

      if
      (
        verifyKernelUpdateNeurons
        (
          true,
          true,
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          modelVariables,
          (*spikeEvents),
          (*connectome),
          commandQueue
        ) != 1
      )
      {
        throw SimException("Neurosim::run: Failed verifyKernelUpdateNeurons");
      }
    }
#elif INJECT_CURRENT_UNTILL_STEP > 0
    if(currentTimeStep < INJECT_CURRENT_UNTILL_STEP)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);

      if
      (
        verifyKernelUpdateNeurons
        (
          true,
          false,
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          modelVariables,
          (*spikeEvents),
          (*connectome),
          commandQueue
        ) != 1
      )
      {
        throw SimException("Neurosim::run: Failed verifyKernelUpdateNeurons");
      }
    }
#else
    if(currentTimeStep < START_PROFILING_AT_STEP)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);
        
      if
      (
        verifyKernelUpdateNeurons
        (
          true,
          false,
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          modelVariables,
          (*spikeEvents),
          (*connectome),
          commandQueue
        ) != 1
      )
      {
        throw SimException("Neurosim::run: Failed verifyKernelUpdateNeurons");
      }
    }
#endif
#endif
/*END verification during inject of spikes or currents*/

/*Unit test verification*/
#elif UPDATE_NEURONS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
      modelVariables);
      
    if
    (
      verifyKernelUpdateNeurons
      (
        true,
        false,
        false,
        currentTimeStep,
        dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
        dataMakeEventPtrsStruct,
        dataGroupEventsTik,
        modelVariables,
        (*spikeEvents),
        (*connectome),
        commandQueue
      ) != 1
    )
    {
      throw SimException("Neurosim::run: Failed verifyKernelUpdateNeurons");
    }
#endif
/*END: verification*/
    }
#endif
/**************************************************************************************************/
#if SIMULATION_SNAPSHOT && \
    (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 &&\
    GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE &&\
    EXPAND_EVENTS_ENABLE && UPDATE_NEURONS_ENABLE_V00)
  if(takeSimSnapshot)
  {
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, modelVariablesBuffer, modelVariablesSizeBytes, 
      modelVariables);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
        
    takeSimulationSnapshot
    (
      currentTimeStep,
      1000,
      dataMakeEventPtrsStruct,
      dataGroupEventsTik,
      modelVariables,
      (*spikeEvents),
      (*synapticEvents),
      commandQueue
    );
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V03
#if ((GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS%2) == 0)
    /*Restore original buffer pointer for the next step*/
    swap2((*operatorSort).dataHistogramGroupEventsTikBuffer, dataHistogramGroupEventsTokBuffer);
    swap2(dataGroupEventsTikBuffer, dataGroupEventsTokBuffer);
#endif
#endif
/**************************************************************************************************/
  }/*for SIMULATION_TIME_STEPS*/
  
  std::cout << "\nDispatched all kernels" << std::endl;
  status = commandQueue.flush();
  ASSERT_CL_SUCCESS(status, "Neurosim::run: cl::CommandQueue.flush failed")
  
  std::cout << "\nEnqueued all kernels" << std::endl;
  std::cout << "\nWaiting for the command queue to complete..." << std::endl;
  eventStatus = CL_QUEUED;
  while(eventStatus != CL_COMPLETE)
  {
    status = ndrEvt.getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS, &eventStatus);
    ASSERT_CL_SUCCESS(status, "Neurosim::run: cl:Event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS)"
      " failed")
  }
  
#if PROFILING_MODE == 1 && START_PROFILING_AT_STEP > -1
  status = commandQueue.finish();
  ASSERT_CL_SUCCESS(status, "Neurosim::run: cl::CommandQueue.finish failed.")
  
  status = ndrEvt.wait();
  GET_TIME_NS(endAppTime);
  ASSERT_CL_SUCCESS(status, "Neurosim::run: commandQueue, cl:Event.wait() failed.")
  
  cl_uint profiledSteps = SIMULATION_TIME_STEPS - START_PROFILING_AT_STEP;
  double appTimeSec = (endAppTime - startAppTime);
  std::cout << "\nDevice timing: " << profiledSteps << " steps, " << appTimeSec << " sec, " 
    << appTimeSec/profiledSteps << " sec/step" << std::endl;
  
  /*Host equivalent computation*/
#if (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 &&\
    GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE &&\
    EXPAND_EVENTS_ENABLE && UPDATE_NEURONS_ENABLE_V00)
  {
    double startAppTime = 0, endAppTime = 0;
    
    std::cout << "\nStarted host timing at step " << START_PROFILING_AT_STEP << std::endl;
    
    GET_TIME_NS(startAppTime);
    for(cl_uint step = START_PROFILING_AT_STEP; step < SIMULATION_TIME_STEPS; step++)
    {
      int result = propagateSpikes
      (
        UPDATE_NEURONS_TOTAL_NEURONS,
        REFERENCE_EVENT_QUEUE_SIZE,
        nrn_ps,
        ne,
        te_ps,
        (*connectome)
      );
      if(result != 1)
      {
        throw SimException("Neurosim::run: propagateSpikes failed.");
      }

      result = updateStep
      (
        IGNORE_SOLVER_EXCEPTIONS,
        UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP,
        step,
        UPDATE_NEURONS_TOTAL_NEURONS,
        UPDATE_NEURONS_PS_ORDER_LIMIT,
        UPDATE_NEURONS_NR_ORDER_LIMIT,
#if (UPDATE_NEURONS_TOLERANCE_MODE == 0)
        UPDATE_NEURONS_NR_ZERO_TOLERANCE
#elif (UPDATE_NEURONS_TOLERANCE_MODE == 1)
        UPDATE_NEURONS_NR_TOLERANCE
#elif (UPDATE_NEURONS_TOLERANCE_MODE > 1)
        UPDATE_NEURONS_NR_ZERO_TOLERANCE
#endif
      );
#if !IGNORE_SOLVER_EXCEPTIONS
      if(result != 1)
      {
        throw SimException("Neurosim::run: updateStep failed.");
      }
#endif
    }
    GET_TIME_NS(endAppTime);
    
    cl_uint profiledSteps = SIMULATION_TIME_STEPS - START_PROFILING_AT_STEP;
    double appTimeSec = (endAppTime - startAppTime);
    std::cout << "\nHost timing: " << profiledSteps << " steps, " << appTimeSec << " sec, " 
      << appTimeSec/profiledSteps << " sec/step" << std::endl;
  }
#endif
#elif KERNEL_LEVEL_PROFILING
  set<std::string>::iterator kernel_name;
  std::cout << "Kernel execution time profile" << std::endl;
  std::cout << "Kernel Name\t\t\tTotal Time (s)\tCount\tAverage Time (ms)" << std::endl 
    << std::endl;
  std::cout << std::fixed;
  double totalExecTime = 0;
  for
  (
    kernel_name = kernelStats.kernelNamesExecTime.begin();
    kernel_name != kernelStats.kernelNamesExecTime.end(); 
    ++kernel_name
  ){
    if(kernelStats.execTime.find(*kernel_name) == kernelStats.execTime.end()){continue;}
    map<std::string, double> execTime = kernelStats.execTime[*kernel_name];
    double averageTime = (1000*execTime["Time"]/execTime["Count"]);
    std::cout << *kernel_name << "\t\t" << std::setprecision(3) << execTime["Time"] << "\t\t" 
      << (cl_uint)execTime["Count"] << "\t" << averageTime << std::endl;
    totalExecTime += execTime["Time"];
    LOG_REP("Kernel " << *kernel_name << " Total Time:" << execTime["Time"]);
    LOG_REP("Kernel " << *kernel_name << " Average Time:" << averageTime);
    LOG_REP("Kernel " << *kernel_name << " Execution Count:" << (cl_uint)execTime["Count"]);
  }
  std::cout << "Total time: " << totalExecTime << std::endl;
  LOG_REP("Total Time:" << totalExecTime);
  std::cout.setf(std::ios::scientific,std::ios::floatfield);
#endif

#if STATISTICS_ENABLE
  GET_TIME_NS(this->runTime);
  this->runTime -= startRunTime;
#endif
}
/**************************************************************************************************/



int 
Neurosim::verifyKernelGroupEventsV00
()
/**************************************************************************************************/
{
  bool result = 1;
#if GROUP_EVENTS_ENABLE_V00
  {
  cl_uint print = 0;
  std::stringstream ss;
  
  cl_uint offset = currentTimeSlot*
    (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + 1);
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_TIME_SLOTS*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE+1);
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, (*synapticEvents).dataHistogram, (*synapticEvents).dataHistogramSizeBytes);
  
  /*Init data for verification*/
  for(cl_uint i = 0; i < dataHistogramGroupEventsVerifySize; i++)
  {
    dataHistogramGroupEventsVerify[i] = 0;
  }
  
  cl_uint total_synaptic_events = 
    (*synapticEvents).dataHistogram[offset+GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
  
  print ? ss << "Total synaptic events: " << total_synaptic_events << std::endl,true : false;
  
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = total_synaptic_events/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI) < 
     total_synaptic_events)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = total_synaptic_event_chunks/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  if(wg_chunk_size*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);

  for(cl_uint i = 0; i < dataGroupEventsTikVerifySize; i++)
  {
    dataGroupEventsTikVerify[i] = 0;
  }
  
  /*Group events for verification*/
  for(cl_uint b = 0; b < GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS; b++)
  {
    cl_uint event_total = (*synapticEvents).dataUnsortedEventCounts[GROUP_EVENTS_TIME_SLOTS*b + currentTimeSlot];
    
    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Compute pointer to event data*/
      cl_uint ptr = 
        /*Event data buffers*/
        b * GROUP_EVENTS_TIME_SLOTS * GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE +
        /*Current event data buffer*/
        currentTimeSlot * GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE +
        /*Current event*/
        e;

      /*Access event*/
      cl_uint key = (*synapticEvents).dataUnsortedEventDelays[ptr];
      
      /*Compute offset key for target neuron*/
      cl_uint bin = (key >> GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00)&
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
        
      /*Offset is based on time slot, bin, WG*/
      cl_uint bin_offset = 
        /*Time slot*/
        currentTimeSlot*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + 1) + 
        /*WG offset*/
        b +
        /*time slot + bin with GROUP_EVENTS_HISTOGRAM_BIN_SIZE as a pitch*/
        bin*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

      /*Check for offset overlap*/
      if(dataOffsetGroupEventsCopy[bin_offset] >= (*synapticEvents).dataHistogram[bin_offset+1])
      {
        std::cout << "verifyKernelGroupEventsV00: Destination event bin pointer overlaps "
          << "with next pointer for bin " << 
          bin << ", time slot " << currentTimeSlot << ", buffer " << b << std::endl;
        result = 0; 
        break;
      }
      
      /*Calculate offset in the grouped data space*/
      cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
        
      /*Compute histogram key for target neuron*/
      cl_uint hist_out_ptr = 
        /*WG offset*/
        (b/GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG) * 
        (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
        /*WG offset for histogram out*/
        (dest_offset/wg_chunk_size) + 
        /*bin*/
        GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
        ((key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V00)&
        GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT);
      /*Verify*/
      if(hist_out_ptr > dataHistogramGroupEventsVerifySize)
      {
        std::cout << "verifyKernelGroupEventsV00: Pointer to an element in output histogram "
          << "is outside of its range" << std::endl;
        result = 0;
        break;
      }
      /*Increment counter for this bin.*/
      dataHistogramGroupEventsVerify[hist_out_ptr]++;
      
      /*Store event at its group location (grouped by bins)*/
      dataGroupEventsTikVerify[dest_offset] = key;
#if GROUP_EVENTS_VALUES_MODE_V00
      dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE+dest_offset] = ptr;
#endif
      /*Increment ptr for next data item*/
      dataOffsetGroupEventsCopy[bin_offset]++;
    }
    if(result != 1){break;}
  }
  
  free(dataOffsetGroupEventsCopy);
  
  print ? ss << "Time slot " << currentTimeSlot << ": " << std::endl, true : false;

  /*Verify event data*/
  for(cl_uint j = 0; j < GROUP_EVENTS_HISTOGRAM_TOTAL_BINS; j++)
  {
    /*Error counters*/
    unsigned long long verify_checksum_target_neuron = 0;
    unsigned long long actual_checksum_target_neuron = 0;
    unsigned long long verify_checksum_value_01 = 0;
    unsigned long long actual_checksum_value_01 = 0;
    unsigned long long verify_error_count_target_neuron = 0;
    unsigned long long actual_error_count_target_neuron = 0;
    
    /*Start and end bin ptr*/
    cl_uint start = (*synapticEvents).dataHistogram[offset + j*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    cl_uint end = (*synapticEvents).dataHistogram[offset + (j+1)*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];

    print ? ss << "Start-End: " << start << "-" << end << std::endl, true : false;

    /*Verify correct bin and checksum in that bin for current time slot*/
    for(cl_uint p = start; p < end; p++)
    {
      cl_uint v_key = (dataGroupEventsTikVerify[p] >> GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00) & 
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
      cl_uint a_key = (dataGroupEventsTik[p] >> GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00) & 
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;

      verify_error_count_target_neuron += (v_key != j);
      (print && (v_key != j)) ? ss << "V," << p << ": (" << v_key << " != " << j << ")" 
        << "; ", true : false;
      actual_error_count_target_neuron += (a_key != j);
      (print && (a_key != j)) ? ss << "A," << p << ": (" << a_key << " != " << j << ")" 
        << "; ", true : false;

      CHECKSUM01(verify_checksum_target_neuron, dataGroupEventsTikVerify[p]);
      /* verify_checksum_target_neuron += dataGroupEventsTikVerify[p]; */
      CHECKSUM01(actual_checksum_target_neuron, dataGroupEventsTik[p]);
      /*actual_checksum_target_neuron += dataGroupEventsTik[p]; */
#if GROUP_EVENTS_VALUES_MODE_V00
      CHECKSUM01(verify_checksum_value_01, dataGroupEventsTikVerify[p + 
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]);
      /*verify_checksum_value_01 += dataGroupEventsTikVerify[p + 
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]; */
      CHECKSUM01(actual_checksum_value_01, dataGroupEventsTik[p + 
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]);
      /*actual_checksum_value_01 += dataGroupEventsTik[p + 
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]; */
#endif
    }
    
    print ? ss << std::endl, true : false;
    
    if(verify_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to validate correct bin for " 
        << verify_error_count_target_neuron 
        << " keys out of " << (end-start) << " in verification data in bin " << j 
        << ", time slot " << currentTimeSlot << std::endl;
      result = 0;
    }
    
    if(actual_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to validate correct bin for " 
        << actual_error_count_target_neuron 
        << " keys out of " << (end-start) << " in actual data in bin " << j 
        << ", time slot " << currentTimeSlot << std::endl;
      result = 0;
    }
    
    if(verify_checksum_target_neuron != actual_checksum_target_neuron)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to match neuron checksum in bin " 
        << j << ", time slot " 
        << currentTimeSlot << std::endl;
      result = 0;
    }
#if GROUP_EVENTS_VALUES_MODE_V00
    if(verify_checksum_value_01 != actual_checksum_value_01)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to match value 01 checksum in bin " 
        << j << ", time slot " 
        << currentTimeSlot << std::endl;
      result = 0;
    }
#endif
    if(print && (result != 1))
    {
      ss << "WG Pointers for bin " << j 
        << "\nWG,Ptr Start,Ptr End,Index,Key Verified,Key Actual,F/P" << std::endl;
        
      for(cl_uint w = 0; w < GROUP_EVENTS_HISTOGRAM_BIN_SIZE; w++)
      {
        cl_uint start = (*synapticEvents).dataHistogram[offset + j*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + w];
        cl_uint end = (*synapticEvents).dataHistogram[offset + j*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + w + 1] - 1;
          
        for(cl_uint p = start; p < end; p++)
        {
          cl_uint v_key = dataGroupEventsTikVerify[p];
          cl_uint a_key = dataGroupEventsTik[p];
          bool pass = (v_key == a_key);
          
          ss << w << "," << start << "," << end << "," << p << "," << v_key << "," << a_key << ","; 
          if(!pass){ss << "F";}
          ss << std::endl;
        }
      }
    }
    
    if(result != 1){break;}
  }
  
  print ? ss << std::endl, true : false;
  
  /*Verify event data*/
  for(cl_uint w = 0; w < (GROUP_EVENTS_GRID_SIZE_WG); w++)
  {
    for(cl_uint b = 0; b < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT); b++)
    {
      /*Checksums have to be the same across each bin*/
      unsigned long long verify_checksum_histogram_out = 0;
      unsigned long long actual_checksum_histogram_out = 0;
      for(cl_uint j = 0; j < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE); j++)
      {
        cl_uint p = 
          /*WG offset*/
          w*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
        CHECKSUM01(verify_checksum_histogram_out, dataHistogramGroupEventsVerify[p]);
        /*verify_checksum_histogram_out += dataHistogramGroupEventsVerify[p]; */
        CHECKSUM01(actual_checksum_histogram_out, (*operatorSort).dataHistogramGroupEventsTik[p]);
        /*actual_checksum_histogram_out += (*operatorSort).dataHistogramGroupEventsTik[p]; */
      }
      
      /*Verify*/
      if(verify_checksum_histogram_out != actual_checksum_histogram_out)
      {
        std::cout << "verifyKernelGroupEventsV00: Failed to match output histogram checksum in bin " 
          << b << ", WG " << w << ", time slot " << currentTimeSlot << std::endl;
        result = 0;
        print ? ss << "Output histogram:" << std::endl, true : false;
        for(cl_uint j = 0; j < (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE); j++)
        {
          cl_uint p = 
            /*WG offset*/
            w*(GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
            /*WG offset for histogram out*/
            j +
            /*bin*/
            b*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
          print ? ss << p << "," << dataHistogramGroupEventsVerify[p] << ","
            << (*operatorSort).dataHistogramGroupEventsTik[p] << std::endl, true : false;
        }
        break;
      }
    }
    if(result != 1){break;}
  }

  if(print){LOG_SIM(ss.str());}
  }
#endif
  return result;
}
/**************************************************************************************************/



int 
Neurosim::verifyKernelGroupEvents
(
  cl_uint verifyKeyBins,
  cl_uint timeSlot,
  cl_uint step,
  cl_uint value1CarryEnable,
  cl_uint value2CarryEnable,
  cl_uint stepShiftEnable,
  cl_uint destinationBufferSize,
  cl_uint histogramBinSize,
  cl_uint histogramTotalBins,
  cl_uint histogramBitMask,
  cl_uint histogramBitShift,
  cl_uint histogramOutTotalBins,
  cl_uint histogramOutTotalGroups,
  cl_uint groupSize,
  cl_uint *dataHistogram,
  cl_uint *dataHistogramOut,
  cl_uint *dataHistogramOutVerify,
  cl_uint *dataGroupedEvents,
  cl_uint *dataGroupedEventsVerify
)
/**************************************************************************************************/
{
  int result = 1;
  bool print = false;
  
  print ? std::cout << "Time slot " << timeSlot << ": " << std::endl, true : false;

  /*Verify event data*/
  for(cl_uint j = 0; j < histogramTotalBins; j++)
  {
    /*Error counters*/
    unsigned long long verify_checksum_target_neuron = 0;
    unsigned long long actual_checksum_target_neuron = 0;
    unsigned long long verify_checksum_value_01 = 0;
    unsigned long long actual_checksum_value_01 = 0;
    unsigned long long verify_checksum_value_02 = 0;
    unsigned long long actual_checksum_value_02 = 0;
    unsigned long long verify_error_count_target_neuron = 0;
    unsigned long long actual_error_count_target_neuron = 0;
    
    /*Start and end bin ptr*/
    cl_uint start = dataHistogram[j*histogramBinSize];
    cl_uint end = dataHistogram[(j+1)*histogramBinSize];
    
    print ? std::cout << "Start-End: " << start << "-" << end << std::endl, true : false;
    
    /*Verify correct bin and checksum in that bin for current time slot*/
    for(cl_uint p = start; p < end; p++)
    {
      cl_uint v_key = 0;
      cl_uint a_key = 0;
      if(stepShiftEnable)
      {
        v_key = (dataGroupedEventsVerify[p]>>(histogramBitShift*step)) & 
        histogramBitMask;
        a_key = (dataGroupedEvents[p]>>(histogramBitShift*step)) & 
        histogramBitMask;
      }
      else
      {
        v_key = (dataGroupedEventsVerify[p]>>histogramBitShift) & 
        histogramBitMask;
        a_key = (dataGroupedEvents[p]>>histogramBitShift) & 
        histogramBitMask;
      }
      
      if(verifyKeyBins)
      {
        verify_error_count_target_neuron += (v_key != j);
        actual_error_count_target_neuron += (a_key != j);
      }
      CHECKSUM01(verify_checksum_target_neuron, dataGroupedEventsVerify[p]);
      /*verify_checksum_target_neuron += dataGroupedEventsVerify[p]; */
      CHECKSUM01(actual_checksum_target_neuron, dataGroupedEvents[p]);
      /*actual_checksum_target_neuron += dataGroupedEvents[p]; */
      if(value1CarryEnable)
      {
        CHECKSUM01(verify_checksum_value_01, dataGroupedEventsVerify[destinationBufferSize + p]);
        /*verify_checksum_value_01 += dataGroupedEventsVerify[destinationBufferSize + p]; */
        CHECKSUM01(actual_checksum_value_01, dataGroupedEvents[destinationBufferSize + p]);
        /*actual_checksum_value_01 += dataGroupedEvents[destinationBufferSize + p]; */
      }
      if(value2CarryEnable)
      {
        CHECKSUM01(verify_checksum_value_02, dataGroupedEventsVerify[2*destinationBufferSize + p]);
        /*verify_checksum_value_02 += dataGroupedEventsVerify[2*destinationBufferSize + p]; */
        CHECKSUM01(actual_checksum_value_02, dataGroupedEvents[2*destinationBufferSize + p]);
        /*actual_checksum_value_02 += dataGroupedEvents[2*destinationBufferSize + p]; */
      }
    }
    
    if(verifyKeyBins)
    {
    if(verify_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEvents: Failed to validate correct bin for " 
        << verify_error_count_target_neuron 
        << " keys out of " << (end-start) << " in verification data in bin " << j 
        << ", time slot " << timeSlot << std::endl;
      result = 0;
      break;
    }
    
    if(actual_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEvents: Failed to validate correct bin for " 
        << actual_error_count_target_neuron 
        << " keys out of " << (end-start) << " in actual data in bin " << j 
        << ", time slot " << timeSlot << std::endl;
      result = 0;
      break;
    }
    }
    
    if(verify_checksum_target_neuron != actual_checksum_target_neuron)
    {
      std::cout << "verifyKernelGroupEvents: Failed to match neuron checksum in bin " 
        << j << ", time slot " << timeSlot << std::endl;
      result = 0;
      break;
    }
    if(value1CarryEnable && (verify_checksum_value_01 != actual_checksum_value_01))
    {
      std::cout << "verifyKernelGroupEvents: Failed to match value 01 checksum in bin " 
        << j << ", time slot " 
        << timeSlot << std::endl;
      result = 0;
      break;
    }
    if(value2CarryEnable && (verify_checksum_value_02 != actual_checksum_value_02))
    {
      std::cout << "verifyKernelGroupEvents: Failed to match value 02 checksum in bin " 
        << j << ", time slot " 
        << timeSlot << std::endl;
      result = 0;
      break;
    }
    /*
    print ? std::cout << verify_checksum_target_neuron << ", ", true : false;
    */
  }
  
  if(result != 1){return result;}
  
  print ? std::cout << std::endl, true : false;
  
  /*Verify event data*/
  if(dataHistogramOutVerify != NULL)
  {
  for(cl_uint w = 0; w < (groupSize); w++)
  {
    for(cl_uint b = 0; b < (histogramOutTotalBins); b++)
    {
      unsigned long long verify_checksum_histogram_out = 0;
      unsigned long long actual_checksum_histogram_out = 0;
      for(cl_uint j = 0; j < (histogramOutTotalGroups); j++)
      {
        cl_uint p = 
          /*WG offset*/
          w*(histogramOutTotalGroups*histogramOutTotalBins) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*histogramOutTotalGroups;
        CHECKSUM01(verify_checksum_histogram_out, dataHistogramOutVerify[p]);
        /*verify_checksum_histogram_out += dataHistogramOutVerify[p]; */
        CHECKSUM01(actual_checksum_histogram_out, dataHistogramOut[p]);
        /*actual_checksum_histogram_out += dataHistogramOut[p]; */
      }
      /*Verify*/
      if(verify_checksum_histogram_out != actual_checksum_histogram_out)
      {
        std::cout << "verifyKernelGroupEvents: Failed to match output histogram checksum in bin "
          << b << ", WG " << w << ", time slot " << currentTimeSlot << std::endl;
        result = 0;
        break;
      }
    }
    if(result != 1){break;}
  }
  }

  return result;
}
/**************************************************************************************************/



#if GROUP_EVENTS_ENABLE_V01
int 
Neurosim::verifyKernelGroupEventsV01(cl_uint step)
/**************************************************************************************************/
{
  int result = 1;
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, (*operatorSort).dataHistogramGroupEventsTik, 
    (*operatorSort).dataHistogramGroupEventsTikSizeBytes);
  
  cl_uint event_total = (*operatorSort).dataHistogramGroupEventsTik[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    
  /*Init data for verification*/
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = event_total/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI) 
    < event_total)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks for WG*/
  cl_uint total_wg_chunks = total_synaptic_event_chunks/GROUP_EVENTS_GRID_SIZE_WG;
  if(total_wg_chunks*GROUP_EVENTS_GRID_SIZE_WG  < total_synaptic_event_chunks)
  {
    total_wg_chunks++;
  }
  
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = total_synaptic_event_chunks/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  if(wg_chunk_size*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);

  for(cl_uint i = 0; i < dataHistogramGroupEventsVerifySize; i++)
  {
    dataHistogramGroupEventsVerify[i] = 0;
  }

  for(cl_uint i = 0; i < dataGroupEventsTikVerifySize; i++)
  {
    dataGroupEventsTikVerify[i] = 0;
  }
    
  /*Group events for verification*/
  for(cl_uint e = 0; e < event_total; e++)
  {
    /*Access event*/
    cl_uint key = dataGroupEventsTik[e];
    cl_uint value = dataGroupEventsTik[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + e];
    
    /*Compute WG which is working on this element*/
    cl_uint wg_id = e/(total_wg_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI));

    /*Compute bin for neuron*/
#if GROUP_EVENTS_ENABLE_STEP_SHIFT_V01
    cl_uint bin = (key>>(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01*step))&
      GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#else
    cl_uint bin = (key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01)&
      GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#endif

    /*Offset is based on time slot, bin, WG*/
    cl_uint bin_offset = 
    /*WG offset*/
    wg_id +
    /*time slot + bin with GROUP_EVENTS_HISTOGRAM_BIN_SIZE as a pitch*/
    bin*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

    /*Check for offset overlap*/
    if(dataOffsetGroupEventsCopy[bin_offset] >= (*operatorSort).dataHistogramGroupEventsTik[bin_offset+1])
    {
      std::cout << "verifyKernelGroupEventsV01: Destination event bin pointer overlaps "
      << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot << ", WG " 
      << wg_id << ", pointer " << e << ", sort step " << step << std::endl;
      result = 0; 
      break;
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];

    /*Compute histogram key for target neuron*/
    cl_uint hist_out_ptr = 
      /*WG offset*/
      (wg_id/GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG) * 
      (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
      /*WG offset for histogram out*/
      dest_offset/wg_chunk_size + 
      /*bin*/
      GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
#if GROUP_EVENTS_ENABLE_STEP_SHIFT_V01
      ((key>>(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01*(step+1)))&
#else
      ((key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01)&
#endif
      GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT);
    /*Verify overflow*/
    if(hist_out_ptr > dataHistogramGroupEventsVerifySize)
    {
      std::cout << "verifyKernelGroupEventsV01: Pointer to an element in output histogram is "
        << "outside of its range" << std::endl;
      result = 0;
      break;
    }
    /*Increment counter for this bin.*/
    dataHistogramGroupEventsVerify[hist_out_ptr]++;
    
    /*Store event at its group location (grouped by bins)*/
    dataGroupEventsTikVerify[dest_offset] = key;
#if GROUP_EVENTS_VALUES_MODE_V01 == 1
    dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = e;
#elif GROUP_EVENTS_VALUES_MODE_V01 == 2
    dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = value;
#endif
    /*Increment ptr for next data item*/
    dataOffsetGroupEventsCopy[bin_offset]++;
  }
  
  free(dataOffsetGroupEventsCopy);
  if(result != 1){return result;}

  result = verifyKernelGroupEvents
  (
    1,
    currentTimeSlot,
    step,
    GROUP_EVENTS_VALUES_MODE_V01,
    0,
    GROUP_EVENTS_ENABLE_STEP_SHIFT_V01,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_GRID_SIZE_WG,
    (*operatorSort).dataHistogramGroupEventsTik,
    dataHistogramGroupEventsTok,
    dataHistogramGroupEventsVerify,
    dataGroupEventsTok,
    dataGroupEventsTikVerify
  );
  
  return result;
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V02
int 
Neurosim::verifyKernelGroupEventsV02(cl_uint step, cl_uint keyOffset)
/**************************************************************************************************/
{
  int result = 1;
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, (*operatorSort).dataHistogramGroupEventsTik, 
    (*operatorSort).dataHistogramGroupEventsTikSizeBytes);
    
  cl_uint event_total = (*operatorSort).dataHistogramGroupEventsTik[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    
  /*Init data for verification*/
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = event_total/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI) 
    < event_total)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks for WG*/
  cl_uint total_wg_chunks = total_synaptic_event_chunks/GROUP_EVENTS_GRID_SIZE_WG;
  if(total_wg_chunks*GROUP_EVENTS_GRID_SIZE_WG  < total_synaptic_event_chunks)
  {
    total_wg_chunks++;
  }
  
  /*Compute total chunks per WG for output histogram*/
  cl_uint wg_chunk_size = total_synaptic_event_chunks/GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  if(wg_chunk_size*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE  < total_synaptic_event_chunks)
  {
    wg_chunk_size++;
  }
  wg_chunk_size *= (GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI);

  for(cl_uint i = 0; i < dataHistogramGroupEventsVerifySize; i++)
  {
    dataHistogramGroupEventsVerify[i] = 0;
  }

  for(cl_uint i = 0; i < dataGroupEventsTikVerifySize; i++)
  {
    dataGroupEventsTikVerify[i] = 0;
  }
    
  /*Group events for verification*/
  for(cl_uint e = 0; e < event_total; e++)
  {
    /*Access event*/
    cl_uint key = dataGroupEventsTik[e];
    cl_uint value = dataGroupEventsTik[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + e];
    cl_uint new_key = (*synapticEvents).dataUnsortedEventTargets[value+keyOffset];
    
    /*Compute WG which is working on this element*/
    cl_uint wg_id = e/(total_wg_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI));

    /*Compute bin for neuron*/
#if GROUP_EVENTS_ENABLE_STEP_SHIFT_V02
    cl_uint bin = (key>>(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02*step))&
      GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#else
    cl_uint bin = (key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02)&
      GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#endif

    /*Offset is based on time slot, bin, WG*/
    cl_uint bin_offset = 
    /*WG offset*/
    wg_id +
    /*time slot + bin with GROUP_EVENTS_HISTOGRAM_BIN_SIZE as a pitch*/
    bin*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

    /*Check for offset overlap*/
    if(dataOffsetGroupEventsCopy[bin_offset] >= (*operatorSort).dataHistogramGroupEventsTik[bin_offset+1])
    {
      std::cout << "verifyKernelGroupEventsV02: Destination event bin pointer overlaps "
      << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot << ", WG " 
      << wg_id << ", pointer " << e << std::endl;
      result = 0; 
      break;
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
      
    /*Compute histogram key for target neuron*/
    cl_uint hist_out_ptr = 
      /*WG offset*/
      (wg_id/GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG) * 
      (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT) + 
      /*WG offset for histogram out*/
      dest_offset/wg_chunk_size + 
      /*bin*/
      GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
#if GROUP_EVENTS_ENABLE_STEP_SHIFT_V02
      ((new_key>>(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V02*(step+1)))&
#else
      ((new_key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V02)&
#endif
      GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT);
    /*Verify overflow*/
    if(hist_out_ptr > dataHistogramGroupEventsVerifySize)
    {
      std::cout << "verifyKernelGroupEventsV02: Pointer to an element in output histogram is "
        << "outside of its range" << std::endl;
      result = 0;
      break;
    }
    /*Increment counter for this bin.*/
    dataHistogramGroupEventsVerify[hist_out_ptr]++;
    
    /*Store event at its group location (grouped by bins)*/
    dataGroupEventsTikVerify[dest_offset] = new_key;
    dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = value;

    /*Increment ptr for next data item*/
    dataOffsetGroupEventsCopy[bin_offset]++;
  }
  
  free(dataOffsetGroupEventsCopy);
  if(result != 1){return result;}

  result = verifyKernelGroupEvents
  (
    0,
    currentTimeSlot,
    step,
    GROUP_EVENTS_VALUES_MODE_V02,
    0,
    GROUP_EVENTS_ENABLE_STEP_SHIFT_V02,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_GRID_SIZE_WG,
    (*operatorSort).dataHistogramGroupEventsTik,
    dataHistogramGroupEventsTok,
    dataHistogramGroupEventsVerify,
    dataGroupEventsTok,
    dataGroupEventsTikVerify
  );

  return result;
}
/**************************************************************************************************/
#endif



#if GROUP_EVENTS_ENABLE_V03
int 
Neurosim::verifyKernelGroupEventsV03(cl_uint step)
/**************************************************************************************************/
{
  int result = 1;
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, (*operatorSort).dataHistogramGroupEventsTik, 
    (*operatorSort).dataHistogramGroupEventsTikSizeBytes);
    
  cl_uint event_total = (*operatorSort).dataHistogramGroupEventsTik[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    
  /*Init data for verification*/
  for(cl_uint i = 0; i < dataGroupEventsTikVerifySize; i++)
  {
    dataGroupEventsTikVerify[i] = 0;
  }
  
  /*Compute total chunks for grid*/
  cl_uint total_synaptic_event_chunks = event_total/(GROUP_EVENTS_WG_SIZE_WI*
    GROUP_EVENTS_ELEMENTS_PER_WI);
  if(total_synaptic_event_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI) 
    < event_total)
  {
    total_synaptic_event_chunks++;
  }
  
  /*Compute total chunks for WG*/
  cl_uint total_wg_chunks = total_synaptic_event_chunks/GROUP_EVENTS_GRID_SIZE_WG;
  if(total_wg_chunks*GROUP_EVENTS_GRID_SIZE_WG  < total_synaptic_event_chunks)
  {
    total_wg_chunks++;
  }
    
  /*Group events for verification*/
  for(cl_uint e = 0; e < event_total; e++)
  {
    /*Access event*/
    cl_uint key = dataGroupEventsTik[e];
    cl_uint valuePtr = dataGroupEventsTik[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + e];

    /*Compute WG which is working on this element*/
    cl_uint wg_id = e/(total_wg_chunks*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI));

    /*Compute bin for neuron*/
#if GROUP_EVENTS_ENABLE_STEP_SHIFT_V03
    cl_uint bin = (key>>(GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03*step))&
      GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#else
    cl_uint bin = (key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03)&
      GROUP_EVENTS_HISTOGRAM_BIN_MASK;
#endif

    /*Offset is based on time slot, bin, WG*/
    cl_uint bin_offset = 
    /*WG offset*/
    wg_id +
    /*time slot + bin with GROUP_EVENTS_HISTOGRAM_BIN_SIZE as a pitch*/
    bin*GROUP_EVENTS_HISTOGRAM_BIN_SIZE;

    /*Check for offset overlap*/
    if(dataOffsetGroupEventsCopy[bin_offset] >= (*operatorSort).dataHistogramGroupEventsTik[bin_offset+1])
    {
      std::cout << "verifyKernelGroupEventsV03: Destination event bin pointer overlaps "
      << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot << ", WG " 
      << wg_id << ", pointer " << e << std::endl;
      result = 0; 
      break;
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
    
    /*Store event at its group location (grouped by bins)*/
    dataGroupEventsTikVerify[dest_offset] = key;
    dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = 
      (*synapticEvents).dataUnsortedEventDelays[valuePtr];
    dataGroupEventsTikVerify[2*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = 
      (*synapticEvents).dataUnsortedEventWeights[valuePtr];
    
    /*Increment ptr for next data item*/
    dataOffsetGroupEventsCopy[bin_offset]++;
  }
  
  free(dataOffsetGroupEventsCopy);
  if(result != 1){return result;}

  result = verifyKernelGroupEvents
  (
    1,
    currentTimeSlot,
    step,
    GROUP_EVENTS_VALUES_MODE_V03,
    1,
    GROUP_EVENTS_ENABLE_STEP_SHIFT_V03,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_GRID_SIZE_WG,
    (*operatorSort).dataHistogramGroupEventsTik,
    dataHistogramGroupEventsTok,
    NULL,
    dataGroupEventsTok,
    dataGroupEventsTikVerify
  );

  return result;
}
/**************************************************************************************************/
#endif



#if MAKE_EVENT_PTRS_ENABLE
int 
Neurosim::verifyKernelMakeEventPtrs()
/**************************************************************************************************/
{
  int result = 1;
  
  /*Load total event for the test*/
  cl_uint totalEvents = dataHistogramGroupEventsTok[MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET];

  /*Compute total chunks in the grid*/
  cl_uint totalEventChunks = totalEvents/
    (MAKE_EVENT_PTRS_WF_SIZE_WI*MAKE_EVENT_PTRS_ELEMENTS_PER_WI);
  if(totalEventChunks*(MAKE_EVENT_PTRS_WF_SIZE_WI*MAKE_EVENT_PTRS_ELEMENTS_PER_WI) < totalEvents)
  {
    totalEventChunks++;
  }
  
  /*Compute total chunks per a WF*/
  cl_uint chunksPerWf = totalEventChunks/(MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF);
  if(chunksPerWf*(MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF)  < totalEventChunks)
  {
    chunksPerWf++;
  }

#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0

  if(result == 1)
  {
    if(!totalEvents)
    {
      for(cl_uint j = 0; j < MAKE_EVENT_PTRS_TOTAL_NEURONS; j++)
      {
        cl_uint ptr = dataMakeEventPtrsStruct[j*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE];
        cl_uint count = dataMakeEventPtrsStruct[j*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1];
        if(count)
        {
          std::cerr << "ERROR, verifyKernelMakeEventPtrs, found non-zero count for neuron ID " 
            << j << ": " << ptr << "->" << count << std::endl;
          result = 0;
          break;
        }
      }

      for(cl_uint j = 0; j < (2*MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF); j++)
      {
        cl_uint nrnId = dataMakeEventPtrsStruct[(MAKE_EVENT_PTRS_TOTAL_NEURONS + j)*
          MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE];
        cl_uint ptr = dataMakeEventPtrsStruct[(MAKE_EVENT_PTRS_TOTAL_NEURONS + j)*
          MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1];
        if(nrnId != MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE || ptr != MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE)
        {
          std::cerr << "ERROR, verifyKernelMakeEventPtrs, found non-zero neuron ID or pointer " 
            << "in inter-WF boundary handler structure for element " << j << ": " << nrnId 
            << "->" << ptr << std::endl;
          result = 0;
          break;
        }
      }
    }
    else
    {
      bool error = false;
      cl_uint totalCountsVerifyOriginal = 0;
      cl_uint previousElementOffset = 0;
      cl_uint previousElement = dataGroupEventsTik[previousElementOffset];
      
      for(cl_uint j = 0; j < totalEvents; j++)
      {
        /*Detect boundary*/
        if(dataGroupEventsTik[j] != previousElement)
        {
          cl_uint offsetVerify = dataMakeEventPtrsStruct[dataGroupEventsTik[j]*
            MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE];
            
          cl_uint countVerify = dataMakeEventPtrsStruct[previousElement*
            MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1];
          
          /*Verify pointer*/
          error = (j != offsetVerify);
          /*Verify count*/
          error |= (j-previousElementOffset != countVerify);

          if(error)
          {
            std::cerr << "ERROR, verifyKernelMakeEventPtrs, mismatched data for neuron ID " 
              << dataGroupEventsTik[j] << ": pointer (" << offsetVerify << " vs " << j 
              << "), count (" << countVerify << " vs " << j-previousElementOffset << ")"
              << std::endl;
            result = 0;
          }
          
          previousElement = dataGroupEventsTik[j];
          previousElementOffset = j;
          totalCountsVerifyOriginal++;
          if(error){break;}
        }
      }
    }
  }

#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

  cl_uint *wfDataIter = (cl_uint *)calloc((MAKE_EVENT_PTRS_GRID_SIZE_WG*MAKE_EVENT_PTRS_WG_SIZE_WF), 
    sizeof(cl_uint));
  
  /*Verify that empty structs (including dummy struct) have 0 counts and limit address stored*/
  cl_uint totalCountsVerifyStruct = 0;
  cl_uint limitAddress = totalEvents;
  
  for(cl_uint i = 0; i < MAKE_EVENT_PTRS_STRUCTS+1; i++)
  {
    cl_uint count = dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*i];
    totalCountsVerifyStruct += count;
    
    cl_uint firstAddress = dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*i + 
      MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE];
    
    if(count == 0 && firstAddress != limitAddress)
    {
      std::cerr << "ERROR, verifyKernelMakeEventPtrs, stored limit address doesn't match actual: "
        << firstAddress << " vs " << limitAddress << " in struct # " << i << " with pointer count " 
        << count << std::endl;
      result = 0;
      break;
    }
  }

  /*Verify that each key in the data has a pointer in pointer structure*/
  cl_uint totalCountsVerifyOriginal = 0;
  cl_uint j = 0;
  cl_uint wfWorkSize = chunksPerWf*(MAKE_EVENT_PTRS_WF_SIZE_WI*MAKE_EVENT_PTRS_ELEMENTS_PER_WI);
  cl_uint wf = j/wfWorkSize;
  cl_uint count = dataMakeEventPtrsStruct[0];
  cl_uint offsetPtr = (wfDataIter[wf]++)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
  cl_uint error = (dataGroupEventsTik[j] != dataMakeEventPtrsStruct[1 + offsetPtr]) && count;
  error += (j != dataMakeEventPtrsStruct[1 + offsetPtr + 1]) && count;
  
  if(!error && (result != 0))
  {
    for(j = 1; j < totalEvents; j++)
    {
      /*Detect boundary*/
      if(dataGroupEventsTik[j] != dataGroupEventsTik[j-1])
      {
        totalCountsVerifyOriginal++;
        
        /*The next pointer struct (WF) is detected*/
        if(wf != j/wfWorkSize)
        {
          /*Verify count for current pointer struct*/
          if(count != wfDataIter[wf])
          {
            std::cerr << "ERROR, verifyKernelMakeEventPtrs, stored count doesn't match actual "
              << "count for WF " << wf << ": " << count << " != " << wfDataIter[wf]
              << std::endl;
            result = 0;
            break;
          }
          
          /*Init count for the next pointer struct*/
          wf = j/wfWorkSize;
          count = dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*wf];
        }
        
        /*Verify pointer*/
        wf = j/wfWorkSize;
        offsetPtr = (wfDataIter[wf]++)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
        error = (dataGroupEventsTik[j] != 
          dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*wf + 1 + offsetPtr]);
        error += (j != 
          dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*wf + 1 + offsetPtr + 1]);
        if(error){break;}
      }
    }
    /*Count first pointer*/
    if(totalCountsVerifyOriginal){totalCountsVerifyOriginal++;}
  }
  
  if(error)
  {
    offsetPtr = (wfDataIter[wf]-1)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
    std::cerr << "ERROR, verifyKernelMakeEventPtrs, unmatched pointer for WF " << wf
      << ", global element ID " << j << ": key " << dataGroupEventsTik[j] << " vs " 
      << dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*wf + 1 + offsetPtr] 
      << ", value " << j << " vs " 
      << dataMakeEventPtrsStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*wf + 1 + offsetPtr + 1] 
      << std::endl;
    result = 0;
  }
  
  /*Verify total count match*/
  if(result != 0 && totalCountsVerifyStruct != totalCountsVerifyOriginal)
  {
    std::cerr << "ERROR, verifyKernelMakeEventPtrs, sum of pointer counters in all structs "
      << "doesn't match total in original data: " << totalCountsVerifyStruct << " vs " 
      << totalCountsVerifyOriginal << std::endl;
    result = 0;
  }

  free(wfDataIter);
#endif

  return result;
}
/**************************************************************************************************/
#endif



/*
  Inject sorted synaptic events into event queue
*/
int 
Neurosim::injectSortedEvents
(
  bool          verify,
  bool          resetEventBuffer,
  cl_uint       totalNeurons,
  cl_uint       eventQueueSize,
  cl_uint       pointersPitch,
  cl_uint       sortedEventsSize,
  cl_uint       *sortedEvents,
  cl_uint       *pointersToEvents,
  neuron_iz_ps  *nrn
)
/**************************************************************************************************/
{
  int result = 1;

  /*Clear input buffers*/
  if(resetEventBuffer)
  {
    for(cl_uint i = 0; i < totalNeurons; i++){nrn[i].n_in = 0;}
  }
  
  for(cl_uint nId = 0; nId < totalNeurons; nId++)
  {
    cl_uint ptr = pointersToEvents[nId*pointersPitch];
    cl_uint count = pointersToEvents[nId*pointersPitch+1];
    if(!count){continue;}
    
    if(verify && count >= eventQueueSize)
    {
      std::cerr << "ERROR, injectSortedEvents, detected event queue overflow for nID " 
        << nId << ": " << count << " >= " << eventQueueSize << std::endl;
      result = 0;
      break;
    }
    
    if(verify && nrn[nId].n_in != 0)
    {
      std::cerr << "ERROR, injectSortedEvents, detected another entry for  " 
        << nId << " with count " << nrn[nId].n_in << std::endl;
      result = 0;
      break;
    }
    
    nrn[nId].n_in = count;

    for(cl_uint j = 0; j < count; j++)
    {
      cl_uint  p = ptr + j;

      if(verify && (p > sortedEventsSize))
      {
        std::cerr 
          << "ERROR, injectSortedEvents, detected address violation in access to sortedEvents: " 
          << p << " > " << sortedEventsSize << std::endl;
        result = 0;
        break;
      }

      cl_float t = (*((cl_float *)(&sortedEvents[p + sortedEventsSize])));
      cl_float w = (*((cl_float *)(&sortedEvents[p + 2*sortedEventsSize])));
      
      if(verify)
      {
        if((p != ptr) && (t < *((cl_float *)(&sortedEvents[p-1 + sortedEventsSize]))))
        {
          std::cerr << "ERROR, injectSortedEvents, detected violation of time sort order for neuron ID " 
            << nId << ": " << t << " < " << *((cl_float *)(&sortedEvents[p-1 + sortedEventsSize])) 
            << std::endl;
          result = 0;
          break;
        }
        
        if(t > 1.0 || t < 0.0)
        {
          std::cerr << "ERROR, injectSortedEvents, detected event time outside of bounds for neuron "
            << nId << ": " << t << " is not within (0.0, 1.0)" << std::endl;
          result = 0;
          break;
        }
      }
      
      nrn[nId].in_t[j] = t;
      nrn[nId].in_w[j] = w;
    }
    if(result != 1){break;}
  }

  return result;
}
/**************************************************************************************************/



/*
  Inject unsorted synaptic events into event queue
*/
int 
Neurosim::injectUnsortedEvents
(
  cl_uint       timeSlots,
  cl_uint       buffers,
  cl_uint       bufferSize,
  cl_uint       eventQueueSize,
  cl_uint       *dataUnsortedEventCounts,
  cl_uint       *dataUnsortedEventTargets,
  cl_uint       *dataUnsortedEventWeights,
  cl_uint       *dataUnsortedEventDelays,
  neuron_iz_ps  *nrn
)
/**************************************************************************************************/
{
  int result = 1;

  for(cl_uint s = 0; s < timeSlots; s++)
  {
    for(cl_uint b = 0; b < buffers; b++)
    {
      cl_uint total = dataUnsortedEventCounts[timeSlots*b + s];
      
      cl_uint ptr = 
        /*Event data buffers*/
        b * timeSlots * 
        bufferSize +
        /*Current event data buffer*/
        s * bufferSize;

      for(cl_uint e = 0; e < total; e++)
      {
        cl_float w = *((cl_float *)(&dataUnsortedEventWeights[ptr + e]));
        cl_float t = *((cl_float *)(&dataUnsortedEventDelays[ptr + e])) + s + 1;
        
        /*Get target neuron and its event count*/
        cl_uint k = dataUnsortedEventTargets[ptr + e];
        unsigned int n_in = nrn[k].n_in;
        
        if(n_in<eventQueueSize)
        {
	        unsigned int j=n_in; 

          /*Use insertion sort to maintain ordered synaptic events*/
          while ((j > 0) && (nrn[k].in_t[j-1] > t))
          {
	          nrn[k].in_t[j] = nrn[k].in_t[j-1]; /*shift*/
	          nrn[k].in_w[j] = nrn[k].in_w[j-1]; 
            j--;
	        }
          
          nrn[k].in_t[j] = t;
          nrn[k].in_w[j] = w;
          nrn[k].n_in++;
        }
        else
        {
          std::cerr << "ERROR, injectUnsortedEvents, detected event queue overflow for nID " 
            << k << ": " << n_in << " >= " << eventQueueSize << std::endl;
          result = 0;
          break;
        }
      }
      if(result != 1){break;}
    }
    if(result != 1){break;}
  }
  
  return result;
}
/**************************************************************************************************/



/*
  Propagation of spike events to synaptic events for PS method
*/
int 
Neurosim::propagateSpikes
(
  unsigned int      totalNeurons,
  unsigned int      eventQueueSize,
  neuron_iz_ps      *nrn,
  int               *ne,
  DATA_TYPE         *te,
  Data_Connectome   &connectome
)
/**************************************************************************************************/
{
  int result = 1;
  
#if !FLIXIBLE_DELAYS_ENABLE
  /*Clear input buffers*/
  for(unsigned int i = 0; i < totalNeurons; i++){nrn[i].n_in = 0;}
#endif
  
  /*Iterate over source neurons*/
  for(unsigned int i = 0; i < totalNeurons; i++)
  {
    /*Detect if source neuron spiked*/
    if(--ne[i] == 0)
    {
      unsigned int ptrEnd = connectome.getSynapseCount(i);

      /*Iterate over target neurons of this source neuron*/
      for(unsigned int s = 0; s < ptrEnd; s++)
      {
        /*Get synapse data*/
        unsigned int n_event = 0; DATA_TYPE t_event = 0; DATA_TYPE w_event = 0;
        connectome.getSynapse(i, s, n_event, t_event, w_event);
        
        unsigned int n_in = nrn[n_event].n_in;
        
        if(n_in<eventQueueSize)
        {
	        unsigned int j=n_in;
          
#if FLIXIBLE_DELAYS_ENABLE
          t_event += te[i];
#else
          t_event = te[i];
#endif
          /*Use insertion sort to maintain ordered synaptic events*/
          while ((j > 0) && (nrn[n_event].in_t[j-1] > t_event))
          {
	          nrn[n_event].in_t[j] = nrn[n_event].in_t[j-1]; /*shift*/
	          nrn[n_event].in_w[j] = nrn[n_event].in_w[j-1]; 
            j--;
	        }
          
          nrn[n_event].in_t[j] = t_event;
          nrn[n_event].in_w[j] = w_event;
          nrn[n_event].n_in++;
        }
        else
        {
          std::cerr << "ERROR, propagateSpikes, detected event queue overflow for neuron ID " 
            << n_event << ": " << n_in << " >= " << eventQueueSize << std::endl;
          result = 0;
          break;
        }
      }
      if(result != 1){break;}
    }
  }
  
#if FLIXIBLE_DELAYS_ENABLE
  if(result == 1)
  {
    /*Decrement synaptic events: */
    for(unsigned int i = 0; i < totalNeurons; i++)
    {
      for(unsigned int j = 0; j < nrn[i].n_in; j++)
      {
        nrn[i].in_t[j] = nrn[i].in_t[j] - 1.0f;
      }
    }
  }
#endif

  return result;
}
/**************************************************************************************************/



/*
  Verification of synaptic event between host and device.
*/
int 
Neurosim::verifyEvents
(
  bool          ignoreWarnings,
  bool          correctWeightPositionMismatch,
  unsigned int  totalNeurons,
  unsigned int  structElementSize,
  unsigned int  sortedEventsSize,
  unsigned int  *pointerStruct,
  unsigned int  *sortedEvents,
  neuron_iz_ps  *nrn
)
/**************************************************************************************************/
{
  int result = 1;

  /*Iterate over neurons*/
  for(unsigned int nH = 0; nH < totalNeurons; nH++)
  {
    /*Get target neuron and its event count*/
    unsigned int eventCount = nrn[nH].n_in;
    
    unsigned int entryAddress = *(pointerStruct + nH*structElementSize);
    unsigned int entryCount = *(pointerStruct + nH*structElementSize + 1);
    
    /*Iterate over events*/
    unsigned int e;
    for(e = 0; e < eventCount; e++)
    {
      DATA_TYPE tH = nrn[nH].in_t[e];
      DATA_TYPE wH = nrn[nH].in_w[e];
      
      if(tH > 1.0){break;}
      
      if(e >= entryCount)
      {
        std::cerr << "ERROR, verifyEvents, event count mismatch 1 for neuron " << nH << ": " 
          << e+1 << " != " << entryCount << std::endl;
        result = 0;
        break;
      }
      
      unsigned int nD = sortedEvents[entryAddress + e];
      DATA_TYPE tD = *((DATA_TYPE *)(&sortedEvents[entryAddress + e + sortedEventsSize]));
      DATA_TYPE wD = *((DATA_TYPE *)(&sortedEvents[entryAddress + e + 2*sortedEventsSize]));
      
      /*Verify neuron IDs*/
      if(nH != nD)
      {
        std::cerr << "ERROR, verifyEvents, neuron ID mismatch: " << nH << " != " << nD << std::endl;
        result = 0;
        break;
      }
      /*Verify event time*/
      if(tH != tD)
      {
        std::cerr << "ERROR, verifyEvents, event time mismatch for neuron " << nH << ": " 
          << tH << " != " << tD << " ("; PRINT_HEX(4, tH); std::cerr << " != "; PRINT_HEX(4, tD);
          std::cerr << ")" << std::endl;
        result = 0;
        break;
      }
      if(result == 0){break;}
      /*Verify event weight*/
      if(wH != wD)
      {
        /*Count events with equal time*/
        unsigned int timeCount = 0;
        for(unsigned int i = e; i < entryCount; i++)
        {
          unsigned int nD1 = sortedEvents[entryAddress + i];
          DATA_TYPE tH1 = nrn[nH].in_t[i];
          DATA_TYPE tD1 = *((DATA_TYPE *)(&sortedEvents[entryAddress + i + sortedEventsSize]));
            
          if((nD1 == nH) && (tD1 == tD) && (tH1 == tH)){timeCount++;}
          else{break;}
        }
        
        /*Verify all weights within the range of event count for the events with equal time*/
        if(timeCount > 1)
        {
          char *weightCheck = (char *)calloc(timeCount, sizeof(char));

          for(unsigned int i = 0; i < timeCount; i++)
          {
            unsigned int p1 = e+i;
            DATA_TYPE wH1 = nrn[nH].in_w[p1];
            DATA_TYPE tH1 = nrn[nH].in_t[p1];
            
            for(unsigned int j = 0; j < timeCount; j++)
            {
              unsigned int p2 = e+j;
              DATA_TYPE wD1 = *((DATA_TYPE *)(&sortedEvents[entryAddress + p2 + 
                2*sortedEventsSize]));
                
              if((wH1 == wD1) && !weightCheck[j])
              {
                weightCheck[j] = 1;
                if(!ignoreWarnings)
                {
                  std::cerr << "WARNING, verifyEvents, event weight position mismatch for neuron " 
                    << nH << ", event time " << tH1 << ": D(" << p2 << " -> " << wD1 << "), H(" 
                    << p1 << " -> " << wH1 << ")" << std::endl;
                }
              }
            }
          }
          
          for(unsigned int i = 0; i < timeCount; i++)
          {
            if(!weightCheck[i])
            {
              DATA_TYPE w = *((DATA_TYPE *)(&sortedEvents[entryAddress + (e+i) + 
                2*sortedEventsSize]));
                
              std::cerr << "ERROR, verifyEvents, unable to find matched weight for neuron " << nH 
                << ": " << w << std::endl;
              result = 0;
              break;
            }
          }

          free(weightCheck);
          
          if(result != 0)
          {
            if(correctWeightPositionMismatch)
            {
              for(unsigned int i = 0; i < timeCount; i++)
              {
                unsigned int p = e+i;
                nrn[nH].in_w[p] = 
                  *((DATA_TYPE *)(&sortedEvents[entryAddress + p + 2*sortedEventsSize]));
              }
            }
            
            e += (timeCount-1);
          }
          else{break;}
        }
        else
        {
          std::cerr << "ERROR, verifyEvents, event weight mismatch for neuron " << nH << ": " 
            << wH << " != " << wD << std::endl;
          result = 0;
          break;
        }
      }
    }
    if(result == 0){break;}
    
    if(result == 1 && e != entryCount)
    {
      std::cerr << "ERROR, verifyEvents, event count mismatch 2 for neuron " << nH << ": " 
        << e << " != " << entryCount << std::endl;
      result = 0;
      break;
    }
  }
  
  return result;
}
/**************************************************************************************************/



/*
  Main stepper routine for PS method on Izhikevich neuron - runs full step
*/
int 
Neurosim::stepIzPs
(
  DATA_TYPE **yp, 
  DATA_TYPE **co, 
  DATA_TYPE *yold, 
  DATA_TYPE *ynew, 
  neuron_iz_ps *nrnp, 
  int *ne, 
  DATA_TYPE *te, 
  int *ip, 
  DATA_TYPE *fp, 
#if STATISTICS_ENABLE
  int *totalSpikes,
  unsigned long long int *psStepCount, 
  DATA_TYPE *mu, 
  int *max_order, 
#endif
  int ps_order_limit,
  int nr_order_limit,
  DATA_TYPE nrTolerance,
  bool variableDelaysEnalbe
)
/**************************************************************************************************/
{
#define PRINT_stepIzPs                                                        (SIMULATION_MODE == 0)
	int result = 1;
  
  DATA_TYPE v,u,g_ampa,g_gaba,chi,E_ampa,E_gaba,t,start;
	DATA_TYPE k,a,b,E,I,vnew,dv,dx,dx_old,dt_part,dt_full,dt; 
	int ps_order,i,j,nrn_ind,steps;
	int err_nv=4, nv=5;
	steps = ip[1];
  nrn_ind = ip[4];

  /*decay_ampa = fp[10];decay_gaba = fp[11];*/
	dt = fp[12]; dt_full = fp[13]; t = fp[15]; start = fp[16];
		
	/*extract variables from neuron structure*/
	v = nrnp->v; u = nrnp->u; g_ampa = nrnp->g_ampa; g_gaba = nrnp->g_gaba;
	E_ampa = nrnp->E_ampa; E_gaba = nrnp->E_gaba; k = nrnp->k;
	I = nrnp->I; E = nrnp->E; a = nrnp->a; b = nrnp->b;
	chi = k*v - g_ampa - g_gaba + nrnp->l;

	/*Set error tolerance */
#if(UPDATE_NEURONS_TOLERANCE_MODE != 0)
  DATA_TYPE tol = fp[17], eta[4];
  eta[0] = tol; eta[1] = tol; eta[2] = tol; eta[3] = tol;
#endif
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
  if(tol != 0){nrTolerance = tol;}
#endif

	yp[0][0] = v;yp[1][0] = u;yp[2][0] = g_ampa;yp[3][0] = g_gaba;yp[4][0] = chi;
	fp[0] = I; fp[1] = k; fp[3] = E_ampa; fp[4] = E_gaba;
	fp[5] = E; fp[6] = a;fp[7] = b; fp[99] = dt_full;

#if (UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD)
  fp[2] = nrnp->v_peak;
#endif
  
  int diverged = 0;
  
  /*integrated step function*/  
  ps_order = ps_step(
#if(UPDATE_NEURONS_TOLERANCE_MODE != 0)
    eta,
#endif
    yp,
    co,
    yold,
    ynew,
    fp,
    iz_first,
    iz_iter,
    ps_order_limit,
    nv,
    err_nv,
    diverged
  );
  
	if(diverged)
  {
#if PRINT_stepIzPs
    std::cerr << "WARNING, stepIzPs, PS step diverged for neuron " << nrn_ind << std::endl;
#endif
    result = 0;
  }
  
	if(ps_order >= ps_order_limit)
  {
#if PRINT_stepIzPs
    std::cerr << "WARNING, stepIzPs, PS order limit overflow: " << ps_order << std::endl;
#endif
    ps_order = ps_order_limit-1;
    result = 0;
  }

#if STATISTICS_ENABLE
  (*psStepCount)++;
	*mu += ((DATA_TYPE)ps_order - *mu)/(DATA_TYPE)*psStepCount;
	if(ps_order > *max_order)*max_order = ps_order;
#endif

  vnew=ynew[0]; /*New membrane voltage value*/

	if (vnew >= nrnp->v_peak) /*rare*/
  {
		yp[0][0] = v - nrnp->v_peak; /*shifted for root finding*/
		dt_part = -yp[0][0]/yp[0][1];	/*First step*/
    dx_old = 100.0f;

    /*Up to nr_order_limit NR iterations */
    for (i = 0; i<nr_order_limit; i++)
    {
			vnew=yp[0][ps_order]*dt_part+yp[0][ps_order-1];
			dv=yp[0][ps_order];

      for(j=ps_order-2;j>=0;j--)
      {
				dv = vnew + dv*dt_part;
				vnew=yp[0][j]+vnew*dt_part;
			}
      
			dx = vnew/dv; 
      dt_part -= dx; 

      if(fabs(dx)<nrTolerance)break;
      if(fabs(dx+dx_old)<nrTolerance)break; 
      dx_old=dx;/*For oscillations*/
		}

    if(i==nr_order_limit)
    {
#if PRINT_stepIzPs
      std::cerr << "WARNING, stepIzPs, failed NR tolerance for neuron ID: " << nrn_ind 
        << std::endl;
#endif
      result = 0;
    }

		if(dt_part>dt_full || dt_part<0)
    {
      dt_part=dt_full/2;
#if PRINT_stepIzPs
      std::cerr << "WARNING, stepIzPs, detected NR divergence for neuron ID: " << nrn_ind 
        << std::endl;
#endif
      result = 0;
    }
    
#if STATISTICS_ENABLE
		/*Record spike and schedule events*/
    (*totalSpikes)++;
#endif

    if(variableDelaysEnalbe)
    {
      /*Events with variable delays have to be scheduled immediatly: */
      ne[nrn_ind]=1; 
    }
    else
    {
      /*Events with fixed delay arrive after specified number of steps: */
      ne[nrn_ind]=steps; 
    }

    te[nrn_ind]=start+dt_part;

		/*Evaluate u, g_ampa, g_gaba at corrected spike time*/
		ps_update(yp,1,ps_order,dt_part,&u);
		ps_update(yp,2,ps_order,dt_part,&g_ampa);
		ps_update(yp,3,ps_order,dt_part,&g_gaba);

		/*post spike updates*/
		v = nrnp->v_reset; u += nrnp->u_step; chi = k*v - g_ampa - g_gaba + nrnp->l;
		yp[0][0] = v;yp[1][0] = u;yp[2][0] = g_ampa;yp[3][0] = g_gaba;
    yp[4][0] = chi;

		dt_part = dt_full-dt_part; fp[99] = dt_part;

    int diverged = 0;
    
    /*new integrated step function*/  
    ps_order = ps_step(
#if(UPDATE_NEURONS_TOLERANCE_MODE != 0)
    eta,
#endif
      yp,
      co,
      yold,
      ynew,
      fp,
      iz_first,
      iz_iter,
      ps_order_limit,
      nv,
      err_nv,
      diverged
    );
    
    if(diverged)
    {
#if PRINT_stepIzPs
      std::cerr << "WARNING, stepIzPs, PS step diverged for neuron " << nrn_ind << std::endl;
#endif
      result = 0;
    }
    
    if(ps_order >= ps_order_limit)
    {
#if PRINT_stepIzPs
      std::cerr << "WARNING, stepIzPs, PS order limit overflow: " << ps_order << std::endl;
#endif
      ps_order = ps_order_limit-1;
      result = 0;
    }

#if STATISTICS_ENABLE
    (*psStepCount)++; 
		*mu += ((DATA_TYPE)ps_order- *mu)/(DATA_TYPE)*psStepCount;
		if(ps_order > *max_order)*max_order = ps_order;
#endif

		nrnp->v=ynew[0]; nrnp->u=ynew[1]; nrnp->g_ampa=ynew[2]; 
    nrnp->g_gaba=ynew[3];
	}
	else
  {
    nrnp->v=ynew[0]; nrnp->u=ynew[1]; nrnp->g_ampa=ynew[2]; 
    nrnp->g_gaba=ynew[3];
	}	
  return result;
  
#undef PRINT_stepIzPs
}
/**************************************************************************************************/



#if UPDATE_NEURONS_ENABLE_V00
int 
Neurosim::updateStep
(
  bool    ignoreFailures,
  int     injectCurrentUntilStep,
  cl_uint currentTimeStep,
  cl_uint totalNeurons,
  int     psOrderLimit,
  int     nrOrderLimit,
  double  nrTolerance
)
/**************************************************************************************************/
{
#define PRINT_updateStep                                                           STATISTICS_ENABLE
  int result = 1;
  int nrn_ind;

#if STATISTICS_ENABLE
  int totalSpikes = 0, max_order_ps = 0;
  unsigned int max_events = 0;
  unsigned long long psStepCount = 0, eventCount = 0;
  DATA_TYPE   mu_order_ps = 0.0;
#endif

  DATA_TYPE t_ps = 0.0, start_ps, w_ps, fp_ps[100];
  neuron_iz_ps   *nrnp_ps, *nrnx_ps;

  /*Runtime storage for variables: */
  DATA_TYPE **yp, *yold, *ynew;
  yold = (DATA_TYPE *)malloc(5*sizeof(DATA_TYPE));  
  ynew = (DATA_TYPE *)malloc(5*sizeof(DATA_TYPE));
  yp = (DATA_TYPE **)malloc(5*sizeof(DATA_TYPE *));  
	for(int i=0;i<5;i++)
    {yp[i] = (DATA_TYPE *)malloc((psOrderLimit+1)*sizeof(DATA_TYPE));}

	int ip[100];
  ip[1]     = steps_ps;
  fp_ps[8]  = co_g_ampa_ps;
  fp_ps[9]  = co_g_gaba_ps;
  fp_ps[12] = dt_ps;
#if(UPDATE_NEURONS_TOLERANCE_MODE == 1)
  fp_ps[17] = tol_ps;
#endif

  nrnx_ps = nrn_ps + totalNeurons;
  
  if((injectCurrentUntilStep > 0) && (currentTimeStep == (cl_uint)injectCurrentUntilStep))
  {
    for(nrnp_ps=nrn_ps; nrnp_ps<nrnx_ps; nrnp_ps++)nrnp_ps->I=0;
#if PRINT_updateStep
  std::cout << "Stopping current injection\n"; 
#endif
  }
  
  t_ps += dt_ps;

  for(nrnp_ps=nrn_ps, nrn_ind=0; nrnp_ps<nrnx_ps; nrnp_ps++, nrn_ind++)
  {
    ip[4] = nrn_ind; 
    /*real time at start of step*/
    fp_ps[15] = t_ps; 
    
    /*Work through substeps separated by synaptic events*/
    /*start time of substep (in [0 1] of whole step)*/
    start_ps = 0; 
    fp_ps[16] = start_ps;
    
#if(UPDATE_NEURONS_TOLERANCE_MODE > 1)
    fp_ps[17] = psTolerance[nrn_ind/(UPDATE_NEURONS_TOTAL_NEURONS/UPDATE_NEURONS_TOLERANCE_CHUNKS)];
#endif
    if(nrnp_ps->n_in)
    {
#if STATISTICS_ENABLE
      if(max_events < nrnp_ps->n_in){max_events = nrnp_ps->n_in;}
#endif
      /*one substep per event*/
#if(FLIXIBLE_DELAYS_ENABLE)
      unsigned int substep = 0;
      DATA_TYPE event = nrnp_ps->in_t[substep];

      while((event <= 1.0f) && (substep < (nrnp_ps->n_in)))
#else
      for(substep = 0; substep < (nrnp_ps->n_in); substep++)
#endif
      {
#if STATISTICS_ENABLE
        eventCount++;
#endif
#if(FLIXIBLE_DELAYS_ENABLE)
        fp_ps[13] = event - start_ps;
#else
        fp_ps[13] = nrnp_ps->in_t[substep] - start_ps;
#endif
        if(fp_ps[13]>0)
        {
          result = stepIzPs(
            yp,
            co,
            yold,
            ynew,
            nrnp_ps,
            ne,
            te_ps,
            ip,
            fp_ps,
#if STATISTICS_ENABLE
            &totalSpikes,
            &psStepCount,
            &mu_order_ps,
            &max_order_ps,
#endif
            psOrderLimit,
            nrOrderLimit,
            (DATA_TYPE)nrTolerance,
            FLIXIBLE_DELAYS_ENABLE
          );
          if(result != 1 && !ignoreFailures){break;}

          start_ps = nrnp_ps->in_t[substep]; 
          fp_ps[16] = start_ps; 
          fp_ps[15] += dt_ps*fp_ps[13];
        }
        
        w_ps = nrnp_ps->in_w[substep];
        
        if(w_ps > 0)
        {
          nrnp_ps->g_ampa+=w_ps; /*AMPA*/
        }
        else 
        {
          nrnp_ps->g_gaba-=w_ps; /*GABA*/
        }
#if(FLIXIBLE_DELAYS_ENABLE)
        substep++;
        event = nrnp_ps->in_t[substep];
#endif
      }
      if(result != 1 && !ignoreFailures){break;}
      
#if(FLIXIBLE_DELAYS_ENABLE)
      /*Shift future synaptic events: */
      nrnp_ps->n_in -= substep;
      
      if(nrnp_ps->n_in > 0)
      {
        unsigned int ptr1 = 0, ptr2 = substep; 
        substep = nrnp_ps->n_in;
        
        while (substep > 0)
        {
          nrnp_ps->in_t[ptr1] = nrnp_ps->in_t[ptr2];
          nrnp_ps->in_w[ptr1] = nrnp_ps->in_w[ptr2];
          ptr1++; ptr2++; substep--;
        }
      }
#endif
      fp_ps[13] = 1-start_ps; /*remainder of time step*/
      
      result = stepIzPs(
        yp,
        co,
        yold,
        ynew,
        nrnp_ps,
        ne,
        te_ps,
        ip,
        fp_ps,
#if STATISTICS_ENABLE
        &totalSpikes,
        &psStepCount,
        &mu_order_ps,
        &max_order_ps,
#endif
        psOrderLimit,
        nrOrderLimit,
        (DATA_TYPE)nrTolerance,
        FLIXIBLE_DELAYS_ENABLE
      );
      if(result != 1 && !ignoreFailures){break;}
    }
    else
    {
      fp_ps[13] = 1;
      
      result = stepIzPs(
        yp,
        co,
        yold,
        ynew,
        nrnp_ps,
        ne,
        te_ps,
        ip,
        fp_ps,
#if STATISTICS_ENABLE
        &totalSpikes,
        &psStepCount,
        &mu_order_ps,
        &max_order_ps,
#endif
        psOrderLimit,
        nrOrderLimit,
        (DATA_TYPE)nrTolerance,
        FLIXIBLE_DELAYS_ENABLE
      );
      if(result != 1 && !ignoreFailures){break;}
    }
  }/*For neurons*/
  
#if PRINT_updateStep && STATISTICS_ENABLE
  std::cout << "\n Spikes: " 
            << "\n  Total: " << totalSpikes 
            << "\n  Average: " << ((double)totalSpikes/(double)totalNeurons)
            << "\n Synaptic events: " 
            << "\n  Max per neuron: " << max_events 
            << "\n ps_step calls: " << psStepCount 
            << "\n Average PS order: " << mu_order_ps 
            << "\n Max PS order: " << max_order_ps 
            << std::endl;
#endif

#if STATISTICS_ENABLE
    averageSpikesInNetworkCounter++;
    averageSpikesInNetwork = averageSpikesInNetwork + (totalSpikes - averageSpikesInNetwork)/
      (averageSpikesInNetworkCounter);

    averageEventsInNetworkCounter++;
    averageEventsInNetwork = averageEventsInNetwork + (eventCount - averageEventsInNetwork)/
      (averageEventsInNetworkCounter);
#endif

  free(yold); free(ynew);
	for(int i=0;i<5;i++){free(yp[i]);} free(yp);

  return result;
#undef PRINT_updateStep
}
/**************************************************************************************************/
#endif



#if UPDATE_NEURONS_ENABLE_V00
int 
Neurosim::verifyKernelUpdateNeurons
(
  bool              verify,
  bool              injectSpikes,
  bool              propagateEvents,
  cl_uint           step,
  cl_uint           sortedEventsSize,
  unsigned int      *pointerStruct,
  unsigned int      *sortedEvents,
  DATA_TYPE         *modelVariables,
  Data_SpikeEvents  &spikeEvents,
  Data_Connectome   &connectome,
  cl::CommandQueue  &queue
)
/**************************************************************************************************/
{
  int result = 1;
  bool breakOnFailure = 1;
  bool ignoreSolverFailuresHost = IGNORE_SOLVER_EXCEPTIONS;

  if(propagateEvents)
  {
    if(injectSpikes)
    {
      memset(ne, 0, UPDATE_NEURONS_TOTAL_NEURONS*sizeof(int));
      
      /*Iterate through spike packets*/
      for(cl_uint packet = 0; packet < UPDATE_NEURONS_SPIKE_PACKETS_V00; packet++)
      {
        cl_uint total_spikes = spikeEvents.getPastSpikeCount(queue, packet);

        /*Iterate through spikes in a current packet and inject spikes*/
        for(cl_uint i = 0; i < total_spikes; i++)
        {
          cl_uint spiked_neuron = 0;
          cl_float spike_time = 0;
          spikeEvents.getPastSpike(queue, packet, i, spiked_neuron, spike_time);
            
          if(ne[spiked_neuron] == 1)
          {
            std::cerr << "ERROR, verifyKernelUpdateNeurons, duplicate entry detected for neuron ID " 
              << spiked_neuron << " while injecting its spike: 1)" << spike_time << ", 2)" 
              << te_ps[spiked_neuron] << std::endl;
            return 0;
          }
          ne[spiked_neuron] = 1;
          te_ps[spiked_neuron] = spike_time;
        }
      }
    }
    
    result = propagateSpikes
    (
      UPDATE_NEURONS_TOTAL_NEURONS,
      REFERENCE_EVENT_QUEUE_SIZE,
      nrn_ps,
      ne,
      te_ps,
      connectome
    );
    if(result != 1){return result;}
    
#if FLIXIBLE_DELAYS_ENABLE
    if((pointerStruct != NULL) && (sortedEvents != NULL) && verify)
    {
      result = verifyEvents
      (
        true,
        true,
        UPDATE_NEURONS_TOTAL_NEURONS,
        UPDATE_NEURONS_STRUCT_ELEMENT_SIZE,
        sortedEventsSize,
        pointerStruct,
        sortedEvents,
        nrn_ps
      );
      
      if(result != 1){return result;}
    }
#endif
  }
  else if((pointerStruct != NULL) && (sortedEvents != NULL))
  {
    result = injectSortedEvents
    (
      true,
      true,
      UPDATE_NEURONS_TOTAL_NEURONS,
      REFERENCE_EVENT_QUEUE_SIZE,
      UPDATE_NEURONS_STRUCT_ELEMENT_SIZE,
      sortedEventsSize,
      sortedEvents,
      pointerStruct,
      nrn_ps
    );
    if(result != 1){return result;}
  }
  
  memset(te_ps, 0, UPDATE_NEURONS_TOTAL_NEURONS*sizeof(DATA_TYPE));
  
  result = updateStep
  (
    ignoreSolverFailuresHost,
    UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP,
    currentTimeStep,
    UPDATE_NEURONS_TOTAL_NEURONS,
    UPDATE_NEURONS_PS_ORDER_LIMIT,
    UPDATE_NEURONS_NR_ORDER_LIMIT,
#if (UPDATE_NEURONS_TOLERANCE_MODE == 0)
    UPDATE_NEURONS_NR_ZERO_TOLERANCE
#elif (UPDATE_NEURONS_TOLERANCE_MODE == 1)
    UPDATE_NEURONS_NR_TOLERANCE
#elif (UPDATE_NEURONS_TOLERANCE_MODE > 1)
    UPDATE_NEURONS_NR_ZERO_TOLERANCE
#endif
  );

  if(result != 1 && !ignoreSolverFailuresHost){return result;}

  if(verify)
  {
  for(cl_uint i=0; i<UPDATE_NEURONS_TOTAL_NEURONS; i++)
  {
    DATA_TYPE underTestType1, reference;
    underTestType1 = modelVariables[i];
    reference = nrn_ps[i].v;
    
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable v for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
    underTestType1 = modelVariables[UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].u;
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable u for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
    underTestType1 = modelVariables[2*UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].g_ampa;
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable g_ampa for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
    underTestType1 = modelVariables[3*UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].g_gaba;
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable g_gaba for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
#if (LOG_MODEL_VARIABLES)
    if(i == LOG_MODEL_VARIABLES_NEURON_ID)
    {
      dataToTraceFile LOG_MODEL_VARIABLES_FILE_BODY(LOG_MODEL_VARIABLES_NEURON_ID);
    }
#endif
    /*
    cl_uint underTestType2 = pointerStruct[i*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1];
    if(underTestType2 != 0)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, event count was not reset for neuron " << i
        << ": " << underTestType2 << std::endl;
      result = 0;
    }
    */
    if(result != 1 && breakOnFailure){break;}
  }

  /*Verify spikes*/
  char *spikeCheck = (char *)calloc(UPDATE_NEURONS_TOTAL_NEURONS, sizeof(char));
  
  /*All verified spikes have equivalents in the reference data*/
  for(cl_uint packet = 0; packet < UPDATE_NEURONS_SPIKE_PACKETS_V00; packet++)
  {
    cl_uint total_spikes = spikeEvents.getSpikeCount(queue, packet);

    /*Iterate through spikes in a current packet*/
    for(cl_uint i = 0; i < total_spikes; i++)
    {
      cl_uint spiked_neuron = 0;
      cl_float spike_time = 0;
      spikeEvents.getSpike(queue, packet, i, spiked_neuron, spike_time);  
        
      if(te_ps[spiked_neuron] != spike_time)
      {
        std::cerr << "ERROR, verifyKernelUpdateNeurons, spike time mismatch for neuron " 
          << spiked_neuron << ", packet " << packet << ": " << te_ps[spiked_neuron] << "!=" 
          << spike_time << std::endl;
        result = 0;
        if(breakOnFailure){break;}
      }
      else
      {
        spikeCheck[spiked_neuron] = 1;
      }
    }
    if(breakOnFailure && (result == 0)){break;}
  }
  
  /*There are no reference spikes, which are not present in the verified data*/
  if(result != 0)
  {
    for(cl_uint n = 0; n < UPDATE_NEURONS_TOTAL_NEURONS; n++)
    {
      if(te_ps[n] != 0.0 && spikeCheck[n] == 0)
      {
        std::cerr << "ERROR, verifyKernelUpdateNeurons, a spike from neuron " << n 
          << " and spike time " << te_ps[n] << " is absent"<< std::endl;
        result = 0;
        if(breakOnFailure){break;}
      }
    }
  }
  
  free(spikeCheck);
  }

  return result;
}
/**************************************************************************************************/
#endif



#if SORT_VERIFY_ENABLE
#define PRINT_VERIFY_SORTED_EVENTS  0
int 
Neurosim::verifySortedEvents
(
  cl_uint *sortedEvents, 
  cl_uint *pointerStruct, 
  cl_uint level
)
/**************************************************************************************************/
{
  int result = 1;
  cl_uint keyOffset = 0;
  cl_uint val1Offset = (keyOffset+1)%GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  cl_uint val2Offset = (keyOffset+2)%GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  char *sortedEventsCheck = NULL;
  
  if(level > 0)
  {
    sortedEventsCheck = (char *)calloc(dataUnsortedEventsSnapShotSize/
      GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS, sizeof(char));
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
    std::cout << "verifySortedEvents " << currentTimeStep << ": started" << std::endl;
    std::cout << "verifySortedEvents " << currentTimeStep << ": found " 
      << dataUnsortedEventsSnapShotSize/GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS 
      << " key-value(s) elements" << std::endl;
#endif
  
  /*Check sort order*/
  for
  (
    cl_uint i = 1; 
    i < dataUnsortedEventsSnapShotSize/GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS; 
    i++
  ){
    /*Verify sorted order for keys*/
    cl_uint v1 = sortedEvents[keyOffset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i];
    cl_uint v2 = sortedEvents[keyOffset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i-1];
    
    if(v1 < v2)
    {
      std::cout << "verifySortedEvents " << currentTimeStep 
        << ": failed to verify sort order for key element " << i << ": " << v1 << "<" << v2 
        << std::endl;
      result = 0;
      break;
    }
    
    /*Verify sorted order for values*/
    if(v1 == v2)
    {
      cl_float time1 = 
        *((cl_float *)(&sortedEvents[val1Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i]));
      cl_float time2 = 
        *((cl_float *)(&sortedEvents[val1Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i-1]));
      
      if(time1 < 0.0f || time1 > (float)(GROUP_EVENTS_MAX_DELAY-GROUP_EVENTS_MIN_DELAY))
      {
        std::cout << "verifySortedEvents: found an event in actual data with delay outside "
          << "of defined range: value " << time1 << ", range " << 0.0f << " - " 
          << (float)(GROUP_EVENTS_MAX_DELAY-GROUP_EVENTS_MIN_DELAY) << std::endl;
        result = 0;
        break;
      }
      
      if(time1 < time2)
      {
        std::cout << "verifySortedEvents " << currentTimeStep 
          << ": failed to verify sort order for value element " << i << ": " << time1 << "<" 
          << time2 << std::endl;
        result = 0;
        break;
      }
    }
    /*
    std::cout  << v1 << "->" 
      << *((cl_float *)(&sortedEvents[val1Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i])) 
      << ", ";
    */
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
  if(result != 0)
  {
    std::cout << "verifySortedEvents " << currentTimeStep << ": verified sort order" << std::endl;
  }
#endif

  if(result != 0 && level > 0)
  {
#if PRINT_VERIFY_SORTED_EVENTS
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1
    /*Determine how many avtive structs*/
    cl_uint activeStructsCount = 0xFFFFFFFF;
    for(cl_uint i = 0; i < MAKE_EVENT_PTRS_STRUCTS; i++)
    {
      if(pointerStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*i] == 0)
      {
        activeStructsCount = i;
        break;
      }
    }

    cl_uint printFreq = 100000;
    std::cout << "verifySortedEvents " << currentTimeStep << ": currently verifying key-value(s)/" 
      << printFreq << ": ";
#endif
#endif
    /*Each key-value(s) elements of dataUnsortedEventsSnapShot must be found in sortedEvents*/
    for
    (
      cl_uint p = 0; 
      p < dataUnsortedEventsSnapShotSize/GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS; 
      p++
    ){
      cl_uint testedElementKey = 
        dataUnsortedEventsSnapShot[keyOffset + p*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS];
      cl_uint entryCount = 0;
      cl_uint entryAddress = 0xFFFFFFFF;
      cl_uint nextEntryAddress = 0xFFFFFFFF;

#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0

      entryAddress = *(pointerStruct + testedElementKey*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE);
      entryCount = *(pointerStruct + testedElementKey*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1);
      nextEntryAddress = entryAddress + entryCount;

#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

      cl_uint structPtr = 0;
      
      /*Determine to which struct testedElementKey belongs to*/
      /*TODO: binary search on sorted array*/
      for(structPtr = 0; structPtr < MAKE_EVENT_PTRS_STRUCTS; structPtr++)
      {
        entryCount = pointerStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*structPtr];
        if(entryCount > 0)
        {
          cl_uint ptr = 
            /*base*/
            MAKE_EVENT_PTRS_STRUCT_SIZE*structPtr + 
            /*place for count*/
            1 + 
            /*last element*/
            (entryCount-1)*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
          if(pointerStruct[ptr] >= testedElementKey){break;}
        }
      }

      /*Search for entry for this element*/
      /*TODO: binary search on sorted array*/
      for(cl_uint s = 0; s < entryCount; s++)
      {
        cl_uint ptr = 
          /*base*/
          MAKE_EVENT_PTRS_STRUCT_SIZE*structPtr + 
          /*place for count*/
          1 + 
          /*key of the element*/
          s*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE;
        if(pointerStruct[ptr] == testedElementKey)
        {
          /*address*/
          entryAddress = pointerStruct[ptr+1];
          /*limit address for this element*/
          if(s < entryCount-1)
          {
            nextEntryAddress = pointerStruct[ptr+MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE + 1];
          }
          else
          {
            nextEntryAddress = pointerStruct[MAKE_EVENT_PTRS_STRUCT_SIZE*(structPtr + 1) + 
              MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE]; 
          }
          break;
        }
      }

      /*Verify that the entry was found*/
      if(entryAddress == 0xFFFFFFFF || nextEntryAddress == 0xFFFFFFFF)
      {
        std::cout << "verifySortedEvents " << currentTimeStep << ": not able to find a struct for "
          << "element " << testedElementKey << " at pointer " << p << std::endl;
        result = 0;
        break;
      }
#endif

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
           sortedEvents[keyOffset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + a])
        {
          if(dataUnsortedEventsSnapShot[val1Offset + p*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS] 
             == sortedEvents[val1Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + a])
          {
            if(dataUnsortedEventsSnapShot[val2Offset + p*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS] 
               == sortedEvents[val2Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + a])
            {
              result = 1;
              sortedEventsCheck[a] = 1;
              break;
            }
          }
        }
      }

      if(result == 0)
      {
        std::cout << "verifySortedEvents " << currentTimeStep << ": not able to find in "
          << "sorted data a combination of " << "key-value(s) for key " << p << "(" 
          << dataUnsortedEventsSnapShot[keyOffset + p*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS] 
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
    if(result != 0)
    {
      std::cout 
        << "\nverifySortedEvents" << currentTimeStep << ": verified mapping of unsorted "
        << "data in sorted data" << std::endl;
    }
#endif
  }
  
  /*Check for unrecognized keys*/
  if(result != 0 && level > 1)
  {
    for(cl_uint i = 0; i < dataUnsortedEventsSnapShotSize/GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS; 
      i++)
    {
      if(!sortedEventsCheck[i])
      {
        std::cerr << "verifySortedEvents " << currentTimeStep << ": found unrecognized element ("
        << i << "): (" << sortedEvents[keyOffset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i] 
        << ")->(" << sortedEvents[val1Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i] 
        << "," << sortedEvents[val2Offset*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + i] << ")" 
        << std::endl;
        /*result = 0;
        break; */
      }
    }
#if PRINT_VERIFY_SORTED_EVENTS
    if(result != 0)
    {
      std::cout << "verifySortedEvents " << currentTimeStep << ": checked for unrecognized keys" 
        << std::endl;
    }
#endif
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
  std::cout << "verifySortedEvents " << currentTimeStep << ": finished" << std::endl;
#endif

  if(dataUnsortedEventsSnapShot){free(dataUnsortedEventsSnapShot);} 
  dataUnsortedEventsSnapShot = NULL;
  
  if(level > 0)
  {
    free(sortedEventsCheck);
  }
  return result;
}
/**************************************************************************************************/
#endif



#if SORT_VERIFY_ENABLE
int 
Neurosim::captureUnsortedEvents
(
  cl_uint *unsortedEventCounts,
  cl_uint *unsortedEventTargets,
  cl_uint *unsortedEventDelays,
  cl_uint *unsortedEventWeights
)
/**************************************************************************************************/
{

  cl_uint ptr_store = 0;
  cl_uint event_total = 0;
  int result = 1;
  
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS; b++)
  {
    event_total += unsortedEventCounts[GROUP_EVENTS_TIME_SLOTS*b + currentTimeSlot];
  }
  
  if(dataUnsortedEventsSnapShot){free(dataUnsortedEventsSnapShot);}
  CALLOC(dataUnsortedEventsSnapShot, cl_uint, event_total*GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS);
  
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS; b++)
  {
    event_total = unsortedEventCounts[GROUP_EVENTS_TIME_SLOTS*b + currentTimeSlot];

    for(cl_uint e = 0; e < event_total; e++)
    {
      /*Compute pointer to event data*/
      cl_uint ptr = 
        /*Event data buffers*/
        b * GROUP_EVENTS_TIME_SLOTS * 
        (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE) +
        /*Current event data buffer*/
        currentTimeSlot * 
        (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE) +
        /*Current event*/
        e;

      /*Check if valid and copy event*/
      cl_uint neuron = unsortedEventTargets[ptr];
      cl_uint time = unsortedEventDelays[ptr];
      cl_uint weight = unsortedEventWeights[ptr];
      
      if(*((cl_float *)(&time)) < 0.0f || 
         *((cl_float *)(&time)) > (float)(GROUP_EVENTS_MAX_DELAY-GROUP_EVENTS_MIN_DELAY))
      {
        std::cout << "captureUnsortedEvents: found an event with delay outside of defined "
          << "range: value " << *((cl_float *)(&time)) << ", range " << 0.0f
          << " - " << (float)(GROUP_EVENTS_MAX_DELAY-GROUP_EVENTS_MIN_DELAY) << std::endl;
        result = 0;
        break;
      }
      
      dataUnsortedEventsSnapShot[ptr_store] = neuron;
      dataUnsortedEventsSnapShot[ptr_store+1] = time;
      dataUnsortedEventsSnapShot[ptr_store+2] = weight;
      ptr_store += GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
    }
  }
  
  if(result != 1)
  {
    free(dataUnsortedEventsSnapShot);
  }
  
  return result;
}
/**************************************************************************************************/
#endif



#if SIMULATION_SNAPSHOT
int 
Neurosim::takeSimulationSnapshot
(
  cl_uint           step,
  cl_uint           sampleSizeNeurons,
  cl_uint           *dataMakeEventPtrsStruct,
  cl_uint           *dataGroupEventsTik,
  cl_float          *modelVariables,
  Data_SpikeEvents       &spikeEvents,
  Data_SynapticEvents    &synapticEvents,
  cl::CommandQueue  &queue
)
/**************************************************************************************************/
{
  SET_RANDOM_SEED(srandSeed, srandCounter);
  LOG_SIM("takeSimulationSnapshot: set srand seed to " << srandSeed);
  
  (dataToSnapshotLogFile).str("");
  dataToSnapshotLogFile << "\n\nSNAPSHOT AT STEP: " << step << "\n\n";
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE && GROUP_EVENTS_ENABLE_V03
  {
    dataToSnapshotLogFile << "\n\nEVENT BUFFERS\n\n";
    dataToSnapshotLogFile << "Time Slot,Parameter Name,Parameter Value" << std::endl;
      
    for(cl_uint s = 0; s < GROUP_EVENTS_TIME_SLOTS; s++)
    {
      double totalEventsMean = 0, totalEventsVariance = 0, percentInh = 0, percentExc = 0;
      int totalEvents = 0, totalEventsMax = 0, 
        totalEventsMin = GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE;
      
      for(cl_uint b = 0; b < GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS; b++)
      {
        int total = synapticEvents.getEventCount(queue, b, s, Data_SynapticEvents::RECENT);
        totalEvents += total;
        if(totalEventsMax < total){totalEventsMax = total;}
        if(totalEventsMin > total){totalEventsMin = total;}
        double totalEventsMeanOld = totalEventsMean;
        totalEventsMean += (total - totalEventsMean)/(b+1);
        if(b > 0){totalEventsVariance += (total - totalEventsMeanOld)*(total - totalEventsMean);}
        
        cl_uint ptr = 
          /*Event data buffers*/
          b * GROUP_EVENTS_TIME_SLOTS * 
          GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE +
          /*Current event data buffer*/
          s * GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE;
        
        for(int e = 0; e < total; e++)
        {
          cl_float w = *((cl_float *)(&(synapticEvents.dataUnsortedEventWeights[ptr + e])));
          if(w < 0){percentInh++;}
          else{percentExc++;}
        }
      }
      totalEventsVariance /= GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS;
      totalEventsVariance = sqrt(totalEventsVariance);
      percentInh = 100.0*percentInh/totalEvents;
      percentExc = 100.0*percentExc/totalEvents;
      
      dataToSnapshotLogFile
        << s << ",Mean," << totalEventsMean << std::endl
        << s << ",Sigma," << totalEventsVariance << std::endl
        << s << ",Max," << totalEventsMax << std::endl
        << s << ",Min," << totalEventsMin << std::endl
        << s << ",Percent Inh," << percentInh << std::endl
        << s << ",Percent Ixc," << percentExc << std::endl;
    }
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02 &&\
    GROUP_EVENTS_ENABLE_V03 && MAKE_EVENT_PTRS_ENABLE && UPDATE_NEURONS_ENABLE_V00

  dataToSnapshotLogFile << "\n\nEVENTS\n\n";
  dataToSnapshotLogFile << "Parameter Name,Parameter Value" << std::endl;

  cl_uint totalEventsExc = 0, totalEventsInh = 0, eventsPerNeuronMax = 0, 
    eventsPerNeuronMin = GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE*
      GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
      
  /*Iterate over neurons*/
  for(unsigned int n = 0; n < UPDATE_NEURONS_TOTAL_NEURONS; n++)
  {
    unsigned int entryAddress = *(dataMakeEventPtrsStruct + n*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE);
    unsigned int entryCount = *(dataMakeEventPtrsStruct + n*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1);
    if(eventsPerNeuronMax < entryCount){eventsPerNeuronMax = entryCount;}
    if(eventsPerNeuronMin > entryCount){eventsPerNeuronMin = entryCount;}
    
    /*Iterate over events*/
    for(unsigned int e = 0; e < entryCount; e++)
    {
      unsigned int nD = dataGroupEventsTik[entryAddress + e];
      DATA_TYPE tD = *((DATA_TYPE *)(&dataGroupEventsTik[entryAddress + e +
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]));
      DATA_TYPE wD = *((DATA_TYPE *)(&dataGroupEventsTik[entryAddress + e +
        2*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]));
      
      if(wD > 0){totalEventsExc++;}
      else{totalEventsInh++;}
    }
  }
  
  dataToSnapshotLogFile << "Total Neurons," << UPDATE_NEURONS_TOTAL_NEURONS << std::endl
    << "Total Inhibitory Events," << totalEventsInh << std::endl
    << "Total Excitatory Events," << totalEventsExc << std::endl
    << "Events Per Neuron Max," << eventsPerNeuronMax << std::endl
    << "Events Per Neuron Min," << eventsPerNeuronMin << std::endl;
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE && UPDATE_NEURONS_ENABLE_V00
  dataToSnapshotLogFile << "\n\nSPIKES\n\n";
  dataToSnapshotLogFile << "Parameter Name,Parameter Value" << std::endl;
  
  cl_uint totalSpikes = 0, spikesPerPacketMax = 0, 
    spikesPerPacketMin = UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE;
  cl_float spikeTimeMax = 0.0f, spikeTimeMin = 100.0f;
  
  /*Iterate through spike packets*/
  for(cl_uint packet = 0; packet < UPDATE_NEURONS_SPIKE_PACKETS_V00; packet++)
  {
    cl_uint spikes = spikeEvents.getSpikeCount(queue, packet);
    
    totalSpikes += spikes;
    if(spikesPerPacketMax < spikes){spikesPerPacketMax = spikes;}
    if(spikesPerPacketMin > spikes){spikesPerPacketMin = spikes;}

    /*Iterate through spikes in a current packet*/
    for(cl_uint i = 0; i < spikes; i++)
    {
      cl_uint spiked_neuron = 0;
      cl_float spike_time = 0;
      spikeEvents.getSpike(queue, packet, i, spiked_neuron, spike_time);

      if(spikeTimeMax < spike_time){spikeTimeMax = spike_time;}
      if(spikeTimeMin > spike_time){spikeTimeMin = spike_time;}
    }
  }
  
  dataToSnapshotLogFile << "Total Spikes," << totalSpikes << std::endl
    << "Packet Spikes Max," << spikesPerPacketMax << std::endl
    << "Packet Spikes Min," << spikesPerPacketMin << std::endl
    << "Spike Time Max," << spikeTimeMax << std::endl
    << "Spike Time Min," << spikeTimeMin << std::endl;

#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
  dataToSnapshotLogFile << "\n\nNEURON MODEL VARIABLES\n\n";
  dataToSnapshotLogFile << "Neuron ID,v,u,g_ampa,g_gaba" << std::endl;
  cl_uint window = UPDATE_NEURONS_TOTAL_NEURONS/sampleSizeNeurons;
  
  for(cl_uint i = 0; i < sampleSizeNeurons; i++)
  {
    cl_uint nId = i*window + cl_uint(abs((window-1)*((double)rand()/((double)RAND_MAX))));
          
    dataToSnapshotLogFile << nId << "," << modelVariables[nId] << "," 
      << modelVariables[UPDATE_NEURONS_TOTAL_NEURONS+nId] << "," 
      << modelVariables[2*UPDATE_NEURONS_TOTAL_NEURONS+nId] << "," 
      << modelVariables[3*UPDATE_NEURONS_TOTAL_NEURONS+nId] << std::endl;
  }
#endif
/**************************************************************************************************/

  snapshotLogFile << (dataToSnapshotLogFile).str();
  (dataToSnapshotLogFile).str("");
  
  return 1;
}
#endif
/**************************************************************************************************/



#if UPDATE_NEURONS_ENABLE_V00
void 
Neurosim::psClean()
/**************************************************************************************************/
{
  for(int i=0;i<4;i++){free(co[i]);} free(co); free(nrn_ps); free(te_ps); free(ne);
}
/**************************************************************************************************/
#endif



#if STATISTICS_ENABLE
void 
Neurosim::printStats()
{
  LOG_REP("Average Events/Iteration:" << this->averageEventsInNetwork);
  LOG_REP("Average Spikes/Iteration:" << this->averageSpikesInNetwork);
  LOG_REP("Total Setup Time:" << this->setupTime);
  LOG_REP("Total Run Time:" << this->runTime);
}
#endif



int 
Neurosim::execute()
{
  std::stringstream excss;
  
  try
  {
    std::cout << "---Getting platform stats------------------------------------------" << std::endl;
    this->getPlatformStats();
      
    std::cout << "---Performing setup------------------------------------------------" << std::endl;
    this->setup();
      
    std::cout << "---Starting execution----------------------------------------------" << std::endl;
    this->run();
    std::cout << "---Finished execution----------------------------------------------" << std::endl;
    
#if (STATISTICS_ENABLE)
    this->printStats();
#endif
    
    LOG_REP("Result:PASS");
  }
  CATCH(excss, execute, LOG_REP("Result:FAIL:" << excss.str()); std::cerr << excss.str(); return 0;)

  return 1;
}



int 
/*main(int argc, char *argv[])*/
main()
{
  std::stringstream excss;
  
  try
  {
    _putenv("GPU_DUMP_DEVICE_KERNEL=3");
  
    Neurosim neuro;
  
    std::cout << "\n";
    
    if(!neuro.execute())
    {
      std::cout << "---RESULT: FAIL----------------------------------------------------\n\n";

      return 0;
    }
  }
  CATCH(excss, main, std::cerr << excss.str(); std::cout 
    << "---RESULT: FAIL----------------------------------------------------\n\n"; return 0;)
  
  std::cout << "---RESULT: PASS----------------------------------------------------\n\n";
  return 1;
}



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
