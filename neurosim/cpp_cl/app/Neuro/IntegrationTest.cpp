
/* ===============================================================================================

  TODO:
  - reduce/restructure: synaptic connectome data structure, time slot data structure, tik/tok
    structures.
  - review buffer flags fore read/write operations
  
  =============================================================================================== */



#include "IntegrationTest.hpp"



#if LOG_SIMULATION && LOG_REPORT
  #define LOG(message, type)\
    if(type == 0)\
    {\
      time_t t; time (&t);\
      *dataToSimulationLogFile << ctime(&t) << " " << message << std::endl;\
    }\
    if(type == 1)\
    {\
      *dataToReportLogFile << message << std::endl;\
    }
#elif LOG_SIMULATION
  #define LOG(message, type)\
    if(type == 0)\
    {\
      time_t t; time (&t);\
      *dataToSimulationLogFile << ctime(&t) << " " << message << std::endl;\
    }
#elif LOG_REPORT
  #define LOG(message, type)\
    if(type == 1)\
    {\
      *dataToReportLogFile << message << std::endl;\
    }
#else
  #define LOG(message, type)
#endif



int
IntegrationTest::allocateHostData()
{
  cl_uint size = 0;
  
/**************************************************************************************************/
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
  {
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
    (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
    (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)||\
    (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
  size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
#endif
#if SCAN_ENABLE_V01
  size = SCAN_HISTOGRAM_BIN_BACKETS*SCAN_HISTOGRAM_BIN_SIZE_V01*
    SCAN_HISTOGRAM_TOTAL_BINS_V01;
#endif

  /* allocate memory for histogram and its verification*/
  CALLOC(dataHistogramGroupEventsTik, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataHistogramGroupEventsTik);
  
  CALLOC(dataHistogramGroupEventsVerify, cl_uint, size);
  
#if (((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) && SCAN_ENABLE_V01) ||\
     ((GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) && SCAN_ENABLE_V01) ||\
     ((GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) && SCAN_ENABLE_V01) ||\
     ((GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) && SCAN_ENABLE_V01))
#if (GROUP_EVENTS_GRID_SIZE_WG != SCAN_HISTOGRAM_BIN_BACKETS)
  #error (GROUP_EVENTS_GRID_SIZE_WG != SCAN_HISTOGRAM_BIN_BACKETS)
#endif
#if (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V01)
  #error (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V01)
#endif
#if (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT != SCAN_HISTOGRAM_TOTAL_BINS_V01)
  #error (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT != SCAN_HISTOGRAM_TOTAL_BINS_V01)
#endif
#endif
  }
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
  {
  /* allocate memory for synaptic events */
#if EXPAND_EVENTS_ENABLE
  size = EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS*EXPAND_EVENTS_TIME_SLOTS;
#endif
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  size = GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS*GROUP_EVENTS_TIME_SLOTS;
#endif

  CALLOC(dataUnsortedEventCounts, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventCounts);
  
  /* allocate memory for synaptic events */
#if EXPAND_EVENTS_ENABLE
  size = EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS * EXPAND_EVENTS_TIME_SLOTS * 
    EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE ;
#endif
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  size = GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS * GROUP_EVENTS_TIME_SLOTS * 
    GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE;
#endif

  CALLOC(dataUnsortedEventTargets, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventTargets);
  CALLOC(dataUnsortedEventDelays, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventDelays);
  CALLOC(dataUnsortedEventWeights, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataUnsortedEventWeights);
  }
#if (EXPAND_EVENTS_ENABLE && GROUP_EVENTS_ENABLE_V00) ||\
    (EXPAND_EVENTS_ENABLE && GROUP_EVENTS_ENABLE_V02) ||\
    (EXPAND_EVENTS_ENABLE && GROUP_EVENTS_ENABLE_V03)
#if (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
  #error (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
#endif
#if (GROUP_EVENTS_TIME_SLOTS != EXPAND_EVENTS_TIME_SLOTS)
  #error (GROUP_EVENTS_TIME_SLOTS != EXPAND_EVENTS_TIME_SLOTS)
#endif
#if (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE != EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
  #error (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE != EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
#endif
#if (GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS != EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
  #error (GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS != EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
#endif
#endif
/**************************************************************************************************/
#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
  {
#if (EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  size = 
    EXPAND_EVENTS_TIME_SLOTS*(EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*EXPAND_EVENTS_GRID_SIZE_WG + 1);
#endif
#if SCAN_ENABLE_V00
  size = 
    SCAN_HISTOGRAM_TIME_SLOTS*(SCAN_HISTOGRAM_TOTAL_BINS_V00*SCAN_HISTOGRAM_BIN_SIZE_V00 + 1);
#endif
#if GROUP_EVENTS_ENABLE_V00
  size = 
    GROUP_EVENTS_TIME_SLOTS*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + 1);
#endif
  /* allocate memory for neuron count/offset histogram and its verification*/
  CALLOC(dataHistogram, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataHistogram);
  CALLOC(dataHistogramVerify, cl_uint, size);
  }
#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) && SCAN_ENABLE_V00)
#if (EXPAND_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
  #error (EXPAND_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
#endif
#if (EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
  #error (EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
#endif
#if (EXPAND_EVENTS_GRID_SIZE_WG != SCAN_HISTOGRAM_BIN_SIZE_V00)
  #error (EXPAND_EVENTS_GRID_SIZE_WG != SCAN_HISTOGRAM_BIN_SIZE_V00)
#endif
#endif
#if (SCAN_ENABLE_V00 && GROUP_EVENTS_ENABLE_V00)
#if (GROUP_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
  #error (GROUP_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
#endif
#if (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
  #error (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
#endif
#if (GROUP_EVENTS_HISTOGRAM_BIN_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V00)
  #error (GROUP_EVENTS_HISTOGRAM_BIN_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V00)
#endif
#endif
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
#if EXPAND_EVENTS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
#if EXPAND_EVENTS_ENABLE
  size = (EXPAND_EVENTS_SPIKE_PACKETS*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS);
#endif
#if UPDATE_NEURONS_ENABLE_V00
  size = (UPDATE_NEURONS_SPIKE_PACKETS_V00*UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS);
#endif
  /* allocate memory for spike data */
  CALLOC(dataSpikePackets, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataSpikePackets);
  }
#if EXPAND_EVENTS_ENABLE && (UPDATE_NEURONS_ENABLE_V00)
#if EXPAND_EVENTS_SPIKE_PACKETS != UPDATE_NEURONS_SPIKE_PACKETS_V00
  #error (EXPAND_EVENTS_SPIKE_PACKETS != UPDATE_NEURONS_SPIKE_PACKETS_V00)
#endif
#if EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS != UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS
  #error (EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS != UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS)
#endif
#endif
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
  {
  /* allocate memory for synaptic counters */
  CALLOC(dataSynapsePointer, cl_uint, EXPAND_EVENTS_SYNAPTIC_POINTER_SIZE);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataSynapsePointer);
  
  /*init synaptic pointer*/
  dataSynapsePointer[0] = 0;
  for(cl_uint i = 1; i < EXPAND_EVENTS_SYNAPTIC_POINTER_SIZE; i++)
  {
    dataSynapsePointer[i] = dataSynapsePointer[i-1] + 
      cl_uint(((double)EXPAND_EVENTS_MAX_SYNAPTIC_DATA_SIZE)*((1.0-SYNAPSE_DEVIATION_RATIO)+
      abs((SYNAPSE_DEVIATION_RATIO*((double)rand()/((double)RAND_MAX))))));
  }

  /* allocate memory for synaptic data */
  size = dataSynapsePointer[EXPAND_EVENTS_SYNAPTIC_POINTER_SIZE-1];
    
  CALLOC(dataSynapseTargets, cl_uint, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataSynapseTargets);
  
  CALLOC(dataSynapseDelays, cl_float, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataSynapseDelays);
  
  CALLOC(dataSynapseWeights, cl_float, size);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataSynapseWeights);

  /* allocate memory for synaptic events for verification*/
  size = EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS * EXPAND_EVENTS_TIME_SLOTS;
  CALLOC(dataUnsortedEventCountsVerify, cl_uint, size);
  
  /* allocate memory for synaptic events for verification*/
  size = EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS * EXPAND_EVENTS_TIME_SLOTS * 
    EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE;
  CALLOC(dataUnsortedEventsTargetsVerify, cl_uint, size);
  CALLOC(dataUnsortedEventsDelaysVerify, cl_uint, size);
  CALLOC(dataUnsortedEventsWeightsVerify, cl_uint, size);

#if (EXPAND_EVENTS_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataExpandEventsDebugHost, cl_uint, EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataExpandEventsDebugHost);
  /* allocate memory for debug device buffer */
  CALLOC(dataExpandEventsDebugDevice, cl_uint, EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataExpandEventsDebugDevice);
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataExpandEventsError, cl_uint, EXPAND_EVENTS_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataExpandEventsError);
#endif
  }
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
#if (SCAN_DEBUG_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataScanDebugHost, cl_uint, SCAN_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataScanDebugHost);
  /* allocate memory for debug device buffer */
  CALLOC(dataScanDebugDevice, cl_uint, SCAN_DEBUG_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataScanDebugDevice);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  /* allocate memory for debug host buffer */
  CALLOC(dataScanError, cl_uint, SCAN_ERROR_BUFFER_SIZE_WORDS);
  REGISTER_MEMORY(KERNEL_ALL, MEM_GLOBAL, dataScanError);
#endif
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
    std::cerr << "ERROR, allocateHostData: (MAKE_EVENT_PTRS_ENABLE and "
      << "UPDATE_NEURONS_ENABLE_V00 pointer struct size mismatch)" << std::endl;
    return SDK_FAILURE;
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
/**************************************************************************************************/

  return SDK_SUCCESS;
}



int 
IntegrationTest::initializeSpikeData
(
  double spikeBufferMinPercentFill,
  double spikeBufferMaxPercentFill
)
/**************************************************************************************************/
{
#if EXPAND_EVENTS_ENABLE
  memset(dataSpikePackets, 0, dataSpikePacketsSizeBytes);
  
  SET_RANDOM_SEED(srandSeed);
  LOG("initializeSpikeData: set srand seed to " << srandSeed, 0);
  
  cl_uint neuronsPerPacket = EXPAND_EVENTS_TOTAL_NEURONS/EXPAND_EVENTS_SPIKE_PACKETS;
  
  /* init spike data */
  for(cl_uint packet = 0; packet < EXPAND_EVENTS_SPIKE_PACKETS; packet++)
  {
    cl_uint packet_index = packet * EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS;
    
    int totalSpikes = EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE;
    GET_RANDOM_INT(totalSpikes, EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE, 
      spikeBufferMinPercentFill, spikeBufferMaxPercentFill);
    if(totalSpikes == -1){return SDK_FAILURE;}

    dataSpikePackets[packet_index] = totalSpikes;
    if(totalSpikes <= 0){continue;}
    
    cl_uint packetNeuronsPerSpikeCount = neuronsPerPacket/totalSpikes;
    
    if(packetNeuronsPerSpikeCount > 2)
    {
      cl_uint currentNeuronIdStart = (neuronsPerPacket*packet);
      for(int i = 0; i < totalSpikes; i++)
      {
        cl_uint spiked_neuron = currentNeuronIdStart +
          cl_uint(abs((packetNeuronsPerSpikeCount-1)*((double)rand()/((double)RAND_MAX))));
          
        cl_float spike_time = 
          cl_float(abs(cl_float(SIMULATION_STEP_SIZE)*((double)rand()/((double)RAND_MAX))));
          
        dataSpikePackets[packet_index + 
          EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * i] 
          = spiked_neuron;
          
        *((cl_float *)(&dataSpikePackets[packet_index + 
          EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * i + 1])) 
          = spike_time;
          
        currentNeuronIdStart += packetNeuronsPerSpikeCount;
      }
    }
    else
    {
      for(int i = 0; i < totalSpikes; i++)
      {
        cl_uint spiked_neuron = (neuronsPerPacket*packet) + i;
          
        cl_float spike_time = 
          cl_float(abs(cl_float(SIMULATION_STEP_SIZE)*((double)rand()/((double)RAND_MAX))));
        
        
        dataSpikePackets[packet_index + 
          EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * i] 
          = spiked_neuron;
          
        *((cl_float *)(&dataSpikePackets[packet_index + 
          EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * i + 1])) 
          = spike_time;
      }
    }
  }
  return SDK_SUCCESS;
#else
  std::cerr << "ERROR, initializeSpikeData: (EXPAND_EVENTS_ENABLE is not true)" 
    << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::initializeExpandEventsData()
/**************************************************************************************************/
{
#if EXPAND_EVENTS_ENABLE
  int result = SDK_SUCCESS;
  
  memset(dataUnsortedEventCounts, 0, dataUnsortedEventCountsSizeBytes);
  memset(dataUnsortedEventTargets, 0, dataUnsortedEventTargetsSizeBytes);
  memset(dataUnsortedEventDelays, 0, dataUnsortedEventDelaysSizeBytes);
  memset(dataUnsortedEventWeights, 0, dataUnsortedEventWeightsSizeBytes);
  memset(dataSynapseTargets, 0, dataSynapseTargetsSizeBytes);
  memset(dataSynapseDelays, 0, dataSynapseDelaysSizeBytes);
  memset(dataSynapseWeights, 0, dataSynapseWeightsSizeBytes);
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  memset(dataHistogram, 0, dataHistogramSizeBytes);
#endif
  
  SET_RANDOM_SEED(srandSeed);
  LOG("initializeExpandEventsData: set srand seed to " << srandSeed, 0);
  
  /*init synaptic data*/
  double gabaRatio = 1.0;
  for(cl_uint i = 0; i < EXPAND_EVENTS_SYNAPTIC_POINTER_SIZE-1; i++)
  {
    cl_uint ptrStart = dataSynapsePointer[i];
    cl_uint ptrEnd = dataSynapsePointer[i+1];
    cl_uint synapseCount = ptrEnd - ptrStart;

    /*Initialize synaptic structure*/
    for(cl_uint j = 0; j < synapseCount; j++)
    {
      cl_uint offset = (ptrStart + j);
      if(offset > dataSynapseTargetsSize)
      {
        std::cerr << "ERROR, initializeExpandEventsData: synapse pointer is outside of "
          << "synapse data range for neuronID " << i << ": " << offset << " > " 
          << dataSynapseTargetsSize << std::endl;
        return SDK_FAILURE;
      }
      
      /*target neuron (avoid direct feedback connections)*/
      cl_uint reinitCount = 100;
      dataSynapseTargets[offset] = i;
      while(dataSynapseTargets[offset] == i && reinitCount)
      {
        reinitCount--;
        dataSynapseTargets[offset] = 
          cl_uint(abs((EXPAND_EVENTS_TOTAL_NEURONS-1)*((double)rand()/((double)RAND_MAX))));
      }
      if(!reinitCount)
      {
        std::cerr << "ERROR, initializeExpandEventsData: failed to generate a connection "
          << "without direct feedback for neuron ID " << i << std::endl;
        return SDK_FAILURE;
      }
      
      /*weight*/
      double weightType = abs(100.0*((double)rand()/((double)RAND_MAX)));
      cl_float weight = 6.0f/1.4f;
      if(weightType < gabaRatio)
      {
        weight = -67.0f/1.4f;
      }
      dataSynapseWeights[offset] = cl_float(weight*((double)rand()/((double)RAND_MAX)));
      
      /*delay*/
      dataSynapseDelays[offset] = cl_float(EXPAND_EVENTS_MIN_DELAY + 
        abs((EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY)*((double)rand()/((double)RAND_MAX))));
    }
  }

  return result;
#else
  std::cerr << "ERROR, initializeExpandEventsData: (EXPAND_EVENTS_ENABLE is not true)" 
    << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int
IntegrationTest::initializeDataForKernelScanHistogramV00(cl_uint timeSlot, cl_uint mode)
/**************************************************************************************************/
{
#if SCAN_ENABLE_V00
  cl_uint runningSum = 0, element = 0;
  
  for(cl_uint j = 0; j < SCAN_HISTOGRAM_TOTAL_BINS_V00; j++)
  {
    cl_uint offset = timeSlot*(SCAN_HISTOGRAM_TOTAL_BINS_V00*SCAN_HISTOGRAM_BIN_SIZE_V00 + 1) + 
      j*SCAN_HISTOGRAM_BIN_SIZE_V00;
    
    for(cl_uint k = 0; k < SCAN_HISTOGRAM_BIN_SIZE_V00; k++)
    {
      if(mode)
      {
        element = dataHistogram[offset + k];
      }
      else
      {
        element = cl_uint(abs(SCAN_HISTOGRAM_MAX_COUNT*((double)rand()/((double)RAND_MAX))));
        dataHistogram[offset + k] = element;
      }

      if(j == 0 && k == 0)
      {
        runningSum = element;
        dataHistogramVerify[offset + k] = 0;
      }
      else
      {
        dataHistogramVerify[offset + k] = runningSum;
        runningSum += element;
      }
    }
  }
  dataHistogramVerify[timeSlot*(SCAN_HISTOGRAM_TOTAL_BINS_V00*SCAN_HISTOGRAM_BIN_SIZE_V00 + 1) + 
    SCAN_HISTOGRAM_TOTAL_BINS_V00*SCAN_HISTOGRAM_BIN_SIZE_V00] = runningSum;
#else
  std::cerr << "ERROR, initializeDataForKernelScanHistogramV00: (SCAN_ENABLE_V00 is not true)" 
    << std::endl;
  return SDK_FAILURE;
#endif

  return SDK_SUCCESS;
}
/**************************************************************************************************/



int 
IntegrationTest::initializeDataForKernelGroupEventsV00
(
  cl_uint   unsortedEventTimeSlotDelta,
  double    percentTimeSlotDeltaDeviation,
  double    percentInh
)
/**************************************************************************************************/
{
#if GROUP_EVENTS_ENABLE_V00
  cl_uint max_offset = 0;

#if (GROUP_EVENTS_DEBUG_ENABLE)
  memset(dataDebugHostGroupEvents, 0, dataDebugHostGroupEventsSizeBytes);
  memset(dataDebugDeviceGroupEvents, 0, dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  memset(dataErrorGroupEvents, 0, dataErrorGroupEventsSizeBytes);
#endif
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
  memset(dataHistogramGroupEventsTik, 0, dataHistogramGroupEventsTikSizeBytes);
#endif
  memset(dataUnsortedEventCounts, 0, dataUnsortedEventCountsSizeBytes);
  memset(dataUnsortedEventTargets, 0, dataUnsortedEventTargetsSizeBytes);
  memset(dataUnsortedEventDelays, 0, dataUnsortedEventDelaysSizeBytes);
  memset(dataUnsortedEventWeights, 0, dataUnsortedEventWeightsSizeBytes);
  memset(dataHistogram, 0, dataHistogramSizeBytes);
  memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);

  /*Initialize syn events and neuron counters*/
  int res = initializeEventBuffers
  (
    GROUP_EVENTS_TIME_SLOTS,
    GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
    GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
    GROUP_EVENTS_TOTAL_NEURONS,
    unsortedEventTimeSlotDelta,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_MIN_DELAY,
    0.0,
    percentTimeSlotDeltaDeviation,
    percentInh,
    dataUnsortedEventCounts,
    dataUnsortedEventTargets,
    dataUnsortedEventDelays,
    dataUnsortedEventWeights,
    dataHistogram
  );
  
  if(res != SDK_SUCCESS ){return SDK_FAILURE;}
  
  return initializeHistogram
  (
    GROUP_EVENTS_TIME_SLOTS,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    dataHistogram
  );

#else
  std::cerr << "ERROR, initializeDataForKernelGroupEventsV00: " 
    << "(GROUP_EVENTS_ENABLE_V00 is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::initializeEventBuffers
(
  cl_uint   timeSlots,
  cl_uint   buffers,
  cl_uint   bufferSize,
  cl_uint   neurons,
  cl_uint   unsortedEventTimeSlotDelta,
  cl_uint   histogramBitShift,
  cl_uint   histogramBinMask,
  cl_uint   histogramTotalBins,
  cl_uint   histogramBinSize,
  double    minDelay,
  double    percentDelayBoarder,
  double    percentTimeSlotDeltaDeviation,
  double    percentInh,
  cl_uint   *dataUnsortedEventCounts,
  cl_uint   *dataUnsortedEventTargets,
  cl_uint   *dataUnsortedEventDelays,
  cl_uint   *dataUnsortedEventWeights,
  cl_uint   *dataHistogram
)
/**************************************************************************************************/
{
  /*Initialize syn events and neuron counters*/
  for(cl_uint s = 0; s < timeSlots; s++)
  {
    for(cl_uint b = 0; b < buffers; b++)
    {
      cl_uint event_total = cl_uint(unsortedEventTimeSlotDelta*
        (timeSlots - s - 1)*(1 - percentTimeSlotDeltaDeviation/100.0) +
        abs(unsortedEventTimeSlotDelta*(timeSlots - s - 1)*
        (2*percentTimeSlotDeltaDeviation/100.0)*(double)rand()/((double)RAND_MAX)));
      
      if(event_total >= bufferSize){event_total = bufferSize;}

      dataUnsortedEventCounts[timeSlots*b + s] = event_total;
      
      cl_uint count_inhibitory = cl_uint(event_total*(percentInh/100.0));
      
      for(cl_uint e = 0; e < event_total; e++)
      {
        /*Compute pointer to event data*/
        cl_uint ptr = 
          /*Event data buffers*/
          b * timeSlots * 
          (bufferSize) +
          /*Current event data buffer*/
          s * 
          (bufferSize) +
          /*Current event*/
          e;

        /*Compute event data*/
        /*target neuron*/
        cl_uint target_neuron = 
          cl_uint(abs((neurons-1)*((double)rand()/((double)RAND_MAX))));
        /*weight*/
        cl_float weight = 6.0f/1.4f;
        if (count_inhibitory != 0)
        {
          count_inhibitory--;
          weight = -67.0f/1.4f;
        }
        cl_float ev_weight = cl_float(weight*((double)rand()/((double)RAND_MAX)));
        /*delay: avoid close proximity to the boarders because of decrement on the host*/
        cl_float ev_delay = 
          cl_float(minDelay*(percentDelayBoarder/100.0) + 
            abs(minDelay*(1.0-2*percentDelayBoarder/100.0)*((double)rand()/((double)RAND_MAX))));

        /*Store event*/
        dataUnsortedEventTargets[ptr] = target_neuron;
        *((cl_float *)(&dataUnsortedEventWeights[ptr])) = ev_weight;
        /*This reduction in accuracy is needed to match decrement operation on the host*/
        ev_delay += (s+1);
        *((cl_float *)(&dataUnsortedEventDelays[ptr])) = ev_delay;
        *((cl_float *)(&dataUnsortedEventDelays[ptr])) -= (s+1);
        
        /*Compute histogram key for target neuron based on MSBs*/
        cl_uint bin = (dataUnsortedEventDelays[ptr]>>histogramBitShift) & 
          histogramBinMask;
        /*Offset is based on time slot, bin, WG*/
        cl_uint offset = 
        /*Time slot*/
        s*(histogramTotalBins*histogramBinSize + 1) + 
        /*WG offset*/
        b +
        /*time slot + bin with histogramBinSize as a pitch*/
        bin*histogramBinSize;

        /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
        dataHistogram[offset]++;
      }
    }
  }

  return SDK_SUCCESS;
}
/**************************************************************************************************/



int 
IntegrationTest::initializeHistogram
(
  cl_uint   timeSlots,
  cl_uint   histogramTotalBins,
  cl_uint   histogramBinSize,
  cl_uint   eventDestinationBufferSize,
  cl_uint   *dataHistogram
)
/**************************************************************************************************/
{
  cl_uint print_bins = 0;
  
  /*Compute offsets based on counters.*/
  cl_uint max_offset = 0;
  
  print_bins ? std::cout << "initializeHistogram: " << 
    "Number of synaptic events in bins: " << std::endl, true : false;
  
  for(cl_uint i = 0; i < timeSlots; i++)
  {
    cl_uint offset = i*(histogramTotalBins*histogramBinSize + 1);
    cl_uint runningSum = 0;
    cl_uint runningSize = 0;
    
    print_bins ? std::cout << "initializeHistogram: Time slot " << i 
      << ": [", true : false;
    
    for(cl_uint j = 0; j < histogramTotalBins*histogramBinSize; j++)
    {
      cl_uint temp = dataHistogram[offset + j];
      dataHistogram[offset + j] = runningSum;
      runningSum += temp;

      if(j%histogramBinSize == 0 && j != 0)
      {
        print_bins ? std::cout << runningSum-runningSize << ", ", true : false;
        runningSize = runningSum;
      }
    }
    
    dataHistogram[offset + histogramTotalBins*histogramBinSize] = 
      runningSum;
    
    print_bins ? std::cout << runningSum-runningSize << "] = " << runningSum << 
      std::endl, true : false;
    
    if(max_offset < runningSum){max_offset = runningSum;}
  }

  if(max_offset > eventDestinationBufferSize)
  {
    std::cout << "initializeHistogram: Destination event buffer overflow. " << 
      "Need to increase eventDestinationBufferSize, which is currently " << 
      eventDestinationBufferSize << " above " << max_offset << std::endl;
    return SDK_FAILURE;
  }
  
  return SDK_SUCCESS;
}
/**************************************************************************************************/



int 
IntegrationTest::initializeDataForKernelScanHistogramV01()
/**************************************************************************************************/
{
#if SCAN_ENABLE_V01
#if (SCAN_DEBUG_ENABLE)
  memset(dataScanDebugHost, 0, dataScanDebugHostSizeBytes);
  memset(dataScanDebugDevice, 0, dataScanDebugDeviceSizeBytes);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  memset(dataScanError, 0, dataScanErrorSizeBytes);
#endif
  memset(dataHistogramGroupEventsTik, 0, dataHistogramGroupEventsTikSizeBytes);
  
  for(cl_uint j = 0; j < (SCAN_HISTOGRAM_BIN_SIZE_V01); j++)
  {
    for(cl_uint b = 0; b < (SCAN_HISTOGRAM_TOTAL_BINS_V01); b++)
    {
      cl_uint sum = 0;
      
      for(cl_uint w = 0; w < (SCAN_HISTOGRAM_BIN_BACKETS); w++)
      {
        cl_uint p = 
          /*WG offset*/
          w*(SCAN_HISTOGRAM_BIN_SIZE_V01*SCAN_HISTOGRAM_TOTAL_BINS_V01) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*SCAN_HISTOGRAM_BIN_SIZE_V01;
          
        dataHistogramGroupEventsTik[p] = 
          cl_uint(abs(SCAN_HISTOGRAM_MAX_COUNT*((double)rand()/((double)RAND_MAX))));
      }
    }
  }
  return SDK_SUCCESS;
#else
  std::cerr << "ERROR, initializeDataForKernelScanHistogramV01: " 
    << "(SCAN_ENABLE_V01 is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int
IntegrationTest::initializeUnsortedEvents
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
  cl_uint *dataUnsortedEventCounts,
  cl_uint *dataUnsortedEventTargets,
  cl_uint *dataUnsortedEventDelays,
  cl_uint *dataUnsortedEventWeights,
  cl_uint *dataHistogram
){
  /*Initialize syn events and neuron counters*/
  for(cl_uint b = 0; b < totalBuffers; b++)
  {
    cl_uint event_total = cl_uint(abs(bufferSize*
      (double)rand()/((double)RAND_MAX)));

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
  return SDK_SUCCESS;
}
/**************************************************************************************************/



int
IntegrationTest::initializeGrouppedEvents
/**************************************************************************************************/
(
  cl_uint enableValues,
  cl_uint totalBuffers,
  cl_uint bufferSize,
  cl_uint destinationBufferSize,
  cl_uint elementSizeWords,
  cl_uint histogramBinSize,
  cl_uint histogramTotalBins,
  cl_uint histogramBitMask,
  cl_uint histogramBitShift,
  cl_uint histogramOutBitShift,
  cl_uint histogramOutTotalGroups,
  cl_uint histogramOutBitMask,
  cl_uint histogramOutSize,
  cl_uint keyOffset,
  cl_uint *dataUnsortedEventCounts,
  cl_uint *dataUnsortedEventTargets,
  cl_uint *dataUnsortedEventDelays,
  cl_uint *dataUnsortedEventWeights,
  cl_uint *dataGroupedEvents,
  cl_uint *dataHistogram,
  cl_uint *dataHistogramOut
){
  int result = SDK_SUCCESS;
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = (histogramTotalBins*histogramBinSize + 1);
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogram, size*sizeof(cl_uint));
  
  /*Init data for verification*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
#endif

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
        result = SDK_FAILURE; 
        break;
      }
      
      /*Calculate offset in the grouped data space*/
      cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
        
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
        result = SDK_FAILURE;
        break;
      }
      /*Increment counter for this bin.*/
      dataHistogramOut[hist_out_ptr]++;
#endif
      
      /*Store event at its group location (grouped by bins)*/
      dataGroupedEvents[dest_offset] = key;
      if(enableValues)
      {
        dataGroupedEvents[destinationBufferSize + dest_offset] = ptr;
      }
      /*Increment ptr for next data item*/
      dataOffsetGroupEventsCopy[bin_offset]++;
    }
    if(result != SDK_SUCCESS){break;}
  }
  free(dataOffsetGroupEventsCopy);
  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::initializeDataForKernelGroupEventsV01(int step, cl_uint keyOffset)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;

#if GROUP_EVENTS_ENABLE_V01
  cl_uint max_offset = 0;
  
  cl_uint shiftFirstStage, shiftNextStage;
  if(step < 0)
  {
    shiftFirstStage = 0, shiftNextStage = 0;
  }
  else
  {
    shiftFirstStage = GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01*step;;
    shiftNextStage = GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01*(step+1);
  }

#if (GROUP_EVENTS_DEBUG_ENABLE)
  memset(dataDebugHostGroupEvents, 0, dataDebugHostGroupEventsSizeBytes);
  memset(dataDebugDeviceGroupEvents, 0, dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  memset(dataErrorGroupEvents, 0, dataErrorGroupEventsSizeBytes);
#endif
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
  memset(dataHistogramGroupEventsTik, 0, dataHistogramGroupEventsTikSizeBytes);
#endif
  memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
  memset(dataHistogramGroupEventsTik, 0, dataHistogramGroupEventsTikSizeBytes);
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
    result = SDK_FAILURE;
  }

  if(result != SDK_FAILURE)
  {
    result = initializeGrouppedEvents
    (
      1,
      GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
      GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
      GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
      GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK,
      shiftFirstStage,
      shiftNextStage,
      GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT,
      dataHistogramGroupEventsTikSize,
      keyOffset,
      dataUnsortedEventCountsTemp,
      dataUnsortedEventTargetsTemp,
      dataUnsortedEventDelaysTemp,
      dataUnsortedEventWeightsTemp,
      dataGroupEventsTik,
      dataHistogramTemp,
      dataHistogramGroupEventsTik
    );
  }

  if(result != SDK_FAILURE)
  {
  runningSum = dataHistogramGroupEventsTik[0];
  dataHistogramGroupEventsTik[0] = 0;
  
  /*Compute offsets*/
  for(cl_uint j = 1; 
    j < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE + 1); j++)
  {
    
    cl_uint d = dataHistogramGroupEventsTik[j];
    dataHistogramGroupEventsTik[j] = runningSum;
    runningSum += d;
  }
  }

  free(dataHistogramTemp);
  free(dataUnsortedEventTargetsTemp);
  free(dataUnsortedEventDelaysTemp);
  free(dataUnsortedEventWeightsTemp);
  free(dataUnsortedEventCountsTemp);
  return result;
#else
  std::cerr << "ERROR, initializeDataForKernelGroupEventsV01: " 
    << "(GROUP_EVENTS_ENABLE_V01 is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::initializeDataForKernelGroupEventsV02_V03
(
  cl_uint step, 
  cl_uint keyOffset
)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;

#if GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  cl_uint shiftFirstStage, shiftNextStage;
  if(step < 0)
  {
    shiftFirstStage = 0, shiftNextStage = 0;
  }
  else
  {
    shiftFirstStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS*step;;
    shiftNextStage = GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT*(step+1);
  }

#if (GROUP_EVENTS_DEBUG_ENABLE)
  memset(dataDebugHostGroupEvents, 0, dataDebugHostGroupEventsSizeBytes);
  memset(dataDebugDeviceGroupEvents, 0, dataDebugDeviceGroupEventsSizeBytes);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  memset(dataErrorGroupEvents, 0, dataErrorGroupEventsSizeBytes);
#endif
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
  memset(dataHistogramGroupEventsTik, 0, dataHistogramGroupEventsTikSizeBytes);
#endif
  memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
  memset(dataUnsortedEventCounts, 0, dataUnsortedEventCountsSizeBytes);
  memset(dataUnsortedEventTargets, 0, dataUnsortedEventTargetsSizeBytes);
  memset(dataUnsortedEventDelays, 0, dataUnsortedEventDelaysSizeBytes);
  memset(dataUnsortedEventWeights, 0, dataUnsortedEventWeightsSizeBytes);
  memset(dataHistogramGroupEventsTik, 0, dataHistogramGroupEventsTikSizeBytes);
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
    dataUnsortedEventCounts,
    dataUnsortedEventTargets,
    dataUnsortedEventDelays,
    dataUnsortedEventWeights,
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
    result = SDK_FAILURE;
  }

  if(result != SDK_FAILURE)
  {
    result = initializeGrouppedEvents
    (
      1,
      GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
      GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
      GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
      GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
      GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK,
      shiftFirstStage,
      shiftNextStage,
      GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
      GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT,
      dataHistogramGroupEventsTikSize,
      keyOffset,
      dataUnsortedEventCounts,
      dataUnsortedEventTargets,
      dataUnsortedEventDelays,
      dataUnsortedEventWeights,
      dataGroupEventsTik,
      dataHistogramTemp,
      dataHistogramGroupEventsTik
    );
  }

  if(result != SDK_FAILURE)
  {
    runningSum = dataHistogramGroupEventsTik[0];
    dataHistogramGroupEventsTik[0] = 0;
    
    /*Compute offsets*/
    for(cl_uint j = 1; 
      j < (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE + 1); j++)
    {
      
      cl_uint d = dataHistogramGroupEventsTik[j];
      dataHistogramGroupEventsTik[j] = runningSum;
      runningSum += d;
    }
  }

  free(dataHistogramTemp);
  return result;
#else
  std::cerr << "ERROR, initializeDataForKernelGroupEventsV02_V03: " 
    << "(GROUP_EVENTS_ENABLE_V02 or GROUP_EVENTS_ENABLE_V03 are not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::initializeSortedEvents
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
  int result = SDK_SUCCESS;
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
        result = SDK_FAILURE;
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
      result = SDK_FAILURE;
    }
  }

  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::initializeDataForKernelMakeEventPtrs
(
  cl_uint mode,
  cl_uint step
)
/**************************************************************************************************/
{
/*
  TODO
  - enhence with test cases:
    - insert key change at the inter WI, inter WF boundaries
*/
#if MAKE_EVENT_PTRS_ENABLE
  int result = SDK_SUCCESS;

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
    eventsPerNeuronDeviation = 30.0;
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
#else
  std::cerr << "ERROR, initializeDataForKernelMakeEventPtrs: " 
    << "(MAKE_EVENT_PTRS_ENABLE is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::initializeEventPointers
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
  int result = SDK_SUCCESS;
  
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
        result = SDK_FAILURE;
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
        result = SDK_FAILURE;
        break;
      }
      if(*((cl_float *)(&sortedEvents[j + totalSortedEvents])) < 
        *((cl_float *)(&sortedEvents[j-1 + totalSortedEvents])))
      {
        std::cerr << "ERROR, initializeEventPointers, detected violation of value sort order " 
          << *((cl_float *)(&sortedEvents[j + totalSortedEvents])) << " < " 
          << *((cl_float *)(&sortedEvents[j-1 + totalSortedEvents])) << std::endl;
        result = SDK_FAILURE;
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



int
IntegrationTest::psInit
(
  cl_uint     totalNeurons,
  cl_uint     injectCurrentUntilStep,
  const char  *neuronVariablesSampleFile
)
/**************************************************************************************************/
{
#if UPDATE_NEURONS_ENABLE_V00
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
    p_connect = 10.0/100.0,
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
  if(neuronVariablesSampleFile != "")
  {
    std::ifstream infile(neuronVariablesSampleFile);
    
    if(!infile.is_open())
    {
      std::cerr << "ERROR, psInit: Failed to open neuronVariablesSampleFile: " << 
        neuronVariablesSampleFile << std::endl;
      return SDK_FAILURE;
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
      return SDK_FAILURE;
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
    if(neuronVariablesSampleFile != "")
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
  
  return SDK_SUCCESS;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::initializeDataForKernelUpdateNeurons
(
  bool          resetEvents,
  bool          resetParameters,
  bool          resetVariables,
  cl_uint       step,
  const char    *neuronVariablesSampleFile
)
/**************************************************************************************************/
{
#if UPDATE_NEURONS_ENABLE_V00
  int result = SDK_SUCCESS;
  
  SET_RANDOM_SEED(srandSeed);
  LOG("initializeDataForKernelUpdateNeurons: set srand seed to " << srandSeed, 0);
  
  memset(dataSpikePackets, 0, dataSpikePacketsSizeBytes);
  
  if(resetVariables || resetParameters)
  {
    result = psInit(UPDATE_NEURONS_TOTAL_NEURONS, UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP, 
      neuronVariablesSampleFile);
    if(result != SDK_SUCCESS){return result;}
  }

  if(resetEvents)
  {
    memset(dataGroupEventsTik, 0, dataGroupEventsTikSizeBytes);
    double gabaRatio = 5.0*(!(step%3));
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
    if(result != SDK_SUCCESS){return result;}
    
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
    if(result != SDK_SUCCESS){return result;}
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
#else
  std::cerr << "ERROR, initializeDataForKernelUpdateNeurons: " 
    << "(UPDATE_NEURONS_ENABLE_V00 is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int
IntegrationTest::registerLocalMemory()
{
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
  {
  cl_uint lmSpikes = sizeof(cl_uint)*EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS; 
  lmSpikesSizeBytes = lmSpikes;
  REGISTER_MEMORY(EXPAND_EVENTS_KERNEL_NAME, MEM_LOCAL, lmSpikes);
  
  cl_uint lmTimeSlotCounters = sizeof(cl_uint)*EXPAND_EVENTS_TIME_SLOTS; 
  lmTimeSlotCountersSizeBytes = lmTimeSlotCounters;
  REGISTER_MEMORY(EXPAND_EVENTS_KERNEL_NAME, MEM_LOCAL, lmTimeSlotCounters);
  
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  cl_uint lmTargetNeuronHistogram = sizeof(cl_uint)*(EXPAND_EVENTS_TIME_SLOTS*
    EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS);
  lmTargetNeuronHistogramSizeBytes = lmTargetNeuronHistogram;
  REGISTER_MEMORY(EXPAND_EVENTS_KERNEL_NAME, MEM_LOCAL, lmTargetNeuronHistogram);
#endif
  }
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
  {
  cl_uint lmScanData = sizeof(cl_uint)*SCAN_WG_SIZE_WI*2; 
  lmScanDataSizeBytes = lmScanData;
  REGISTER_MEMORY(SCAN_KERNEL_NAME, MEM_LOCAL, lmScanData);
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  {
  cl_uint lmSortData = sizeof(cl_uint)*(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI+16); 
  lmSortDataSizeBytes = lmSortData;
  REGISTER_MEMORY(GROUP_EVENTS_KERNEL_NAME, MEM_LOCAL, lmSortData);

  cl_uint lmlocalHistogramReference = sizeof(cl_uint)*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS; 
  lmlocalHistogramReferenceSizeBytes = lmlocalHistogramReference;
  REGISTER_MEMORY(GROUP_EVENTS_KERNEL_NAME, MEM_LOCAL, lmlocalHistogramReference);
  
  cl_uint lmlocalHistogramScratchPad = sizeof(cl_uint)*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*2); 
  lmlocalHistogramScratchPadSizeBytes = lmlocalHistogramScratchPad;
  REGISTER_MEMORY(GROUP_EVENTS_KERNEL_NAME, MEM_LOCAL, lmlocalHistogramScratchPad);
  
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
/**************************************************************************************************/

  return SDK_SUCCESS;
}



int
IntegrationTest::getPlatformStats()
{
  cl_int err;

  // Plaform info
  err = cl::Platform::get(&platforms);

  if (err != CL_SUCCESS) 
  {
    std::cerr << "ERROR: " <<  "cl::Platform::get()" << " (" << err << ")" << std::endl;
    return SDK_FAILURE;
  }

  // Iteratate over platforms
  std::cout << "Number of platforms:\t\t\t\t " 
            << platforms.size() 
            << std::endl;
  for (std::vector<cl::Platform>::iterator i = platforms.begin(); 
       i != platforms.end(); 
       ++i) {
      std::cout << "  Plaform Profile:\t\t\t\t "    
                << (*i).getInfo<CL_PLATFORM_PROFILE>().c_str() 
                << std::endl; 
      std::cout << "  Plaform Version:\t\t\t\t "    
                << (*i).getInfo<CL_PLATFORM_VERSION>().c_str() 
                << std::endl; 
      std::cout << "  Plaform Name:\t\t\t\t\t "     
                << (*i).getInfo<CL_PLATFORM_NAME>().c_str() 
                << std::endl; 
      std::cout << "  Plaform Vendor:\t\t\t\t "   
                << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << std::endl; 
      if ((*i).getInfo<CL_PLATFORM_EXTENSIONS>().size() > 0) {
          std::cout << "  Plaform Extensions:\t\t\t " 
                    << (*i).getInfo<CL_PLATFORM_EXTENSIONS>().c_str() 
                    << std::endl; 
      }
  }

  std::cout << std::endl << std:: endl;
  // Now Iteratate over each platform and its devices
  for (std::vector<cl::Platform>::iterator p = platforms.begin(); 
       p != platforms.end(); 
       ++p) {

      std::cout << "  Plaform Name:\t\t\t\t\t "     
                << (*p).getInfo<CL_PLATFORM_NAME>().c_str() 
                << std::endl; 

      std::vector<cl::Device> devices;
      (*p).getDevices(CL_DEVICE_TYPE_ALL, &devices);
  
      std::cout << "Number of devices:\t\t\t\t " << devices.size() << std::endl;
      for (std::vector<cl::Device>::iterator i = devices.begin(); 
           i != devices.end(); 
           ++i) {
          
          std::cout << "  Device Type:\t\t\t\t\t " ;
          cl_device_type dtype = (*i).getInfo<CL_DEVICE_TYPE>();
          switch (dtype) 
          {
            case CL_DEVICE_TYPE_ACCELERATOR:
                std::cout << "CL_DEVICE_TYPE_ACCRLERATOR" << std::endl;
                break;
            case CL_DEVICE_TYPE_CPU:
                std::cout << "CL_DEVICE_TYPE_CPU" << std::endl;
                break;
            case CL_DEVICE_TYPE_DEFAULT:
                std::cout << "CL_DEVICE_TYPE_DEFAULT" << std::endl;
                break;
            case CL_DEVICE_TYPE_GPU:
                std::cout << "CL_DEVICE_TYPE_GPU" << std::endl;
                break;
          }

          std::cout << "  Device ID:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_VENDOR_ID>() 
                    << std::endl;
          
          std::cout << "  Max compute units:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() 
                    << std::endl;
          
          std::cout << "  Max work items dimensions:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() 
                    << std::endl;
          
          cl::detail::param_traits<cl::detail::cl_device_info,CL_DEVICE_MAX_WORK_ITEM_SIZES>::
            param_type witems = (*i).getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
          for (cl_uint x = 0; x < (*i).getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>(); x++) 
          {
              std::cout << "    Max work items[" 
                        << x << "]:\t\t\t\t " 
                        << witems[x] 
                        << std::endl;
          }

          std::cout << "  Max work group size:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() 
                    << std::endl;
          
          std::cout << "  Preferred vector width char:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR>() 
                    << std::endl;

          std::cout << "  Preferred vector width short:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT>() 
                    << std::endl;
          
          std::cout << "  Preferred vector width int:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT>() 
                    << std::endl;
          
          std::cout << "  Preferred vector width long:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG>() 
                    << std::endl;
          
          std::cout << "  Preferred vector width float:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT>() 
                    << std::endl;
          
          std::cout << "  Preferred vector width double:\t\t " 
                    << (*i).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE>() 
                    << std::endl;
          
          std::cout << "  Max clock frequency:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() 
                    << "Mhz"
                    << std::endl;
          
          std::cout << "  Address bits:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_ADDRESS_BITS>() 
                    << std::endl;        
          
          std::cout << "  Max memeory allocation:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() 
                    << std::endl;        
          
          std::cout << "  Image support:\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_IMAGE_SUPPORT>() ? "Yes" : "No")
                    << std::endl;        
          
          if ((*i).getInfo<CL_DEVICE_IMAGE_SUPPORT>()) {
              std::cout << "  Max number of images read arguments:\t " 
                        << (*i).getInfo<CL_DEVICE_MAX_READ_IMAGE_ARGS>()
                        << std::endl;        

              std::cout << "  Max number of images write arguments:\t " 
                        << (*i).getInfo<CL_DEVICE_MAX_WRITE_IMAGE_ARGS>()
                        << std::endl;        
              
              std::cout << "  Max image 2D width:\t\t\t " 
                        << (*i).getInfo<CL_DEVICE_IMAGE2D_MAX_WIDTH>()
                        << std::endl;        

              std::cout << "  Max image 2D height:\t\t\t " 
                        << (*i).getInfo<CL_DEVICE_IMAGE2D_MAX_HEIGHT>()
                        << std::endl;        
              
              std::cout << "  Max image 3D width:\t\t\t " 
                        << (*i).getInfo<CL_DEVICE_IMAGE3D_MAX_WIDTH>()
                        << std::endl;        

              std::cout << "  Max image 3D height:\t " 
                        << (*i).getInfo<CL_DEVICE_IMAGE3D_MAX_HEIGHT>()
                        << std::endl;        
              
              std::cout << "  Max image 3D depth:\t\t\t " 
                        << (*i).getInfo<CL_DEVICE_IMAGE3D_MAX_DEPTH>()
                        << std::endl;        

              std::cout << "  Max samplers within kernel:\t\t " 
                        << (*i).getInfo<CL_DEVICE_MAX_SAMPLERS>()
                        << std::endl;        
          }

          std::cout << "  Max size of kernel argument:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_PARAMETER_SIZE>()
                    << std::endl;        
          
          std::cout << "  Alignment (bits) of base address:\t\t " 
                    << (*i).getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>()
                    << std::endl;        
          
          std::cout << "  Minimum alignment (bytes) for any datatype:\t " 
                    << (*i).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>()
                    << std::endl;        

          std::cout << "  Single precision floating point capability" << std::endl;
          std::cout << "    Denorms:\t\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() & 
                        CL_FP_DENORM ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Quiet NaNs:\t\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() & 
                        CL_FP_INF_NAN ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Round to nearest even:\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
                        CL_FP_ROUND_TO_NEAREST ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Round to zero:\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
                        CL_FP_ROUND_TO_ZERO ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Round to +ve and infinity:\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
                        CL_FP_ROUND_TO_INF ? "Yes" : "No")
                    << std::endl;
          std::cout << "    IEEE754-2008 fused multiply-add:\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
                        CL_FP_FMA ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Correctly rounded div and sqrt:\t\t " 
                    << ((*i).getInfo<CL_DEVICE_SINGLE_FP_CONFIG>() &  
                        CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT ? "Yes" : "No")
                    << std::endl;
                    

          std::cout << "  Cache type:\t\t\t\t\t " ;
          switch ((*i).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>()) {
          case CL_NONE:
              std::cout << "None" << std::endl;
              break;
          case CL_READ_ONLY_CACHE:
              std::cout << "Read only" << std::endl;
              break;
          case CL_READ_WRITE_CACHE:
              std::cout << "Read/Write" << std::endl;
              break;
          }
          
          std::cout << "  Cache line size:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>()
                    << std::endl;
          
          std::cout << "  Cache size:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>()
                    << std::endl;
          
          std::cout << "  Global memory size:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()
                    << std::endl;
          
          std::cout << "  Constant buffer size:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>()
                    << std::endl;
          
          std::cout << "  Max number of constant args:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_MAX_CONSTANT_ARGS>()
                    << std::endl;

          std::cout << "  Local memory type:\t\t\t\t " ;
          switch ((*i).getInfo<CL_DEVICE_LOCAL_MEM_TYPE>()) {
          case CL_LOCAL:
              std::cout << "Scratchpad" << std::endl;
              break;
          case CL_GLOBAL:
              std::cout << "Global" << std::endl;
              break;
          }
          

          std::cout << "  Local memory size:\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()
                    << std::endl;
          
          std::cout << "  Profiling timer resolution:\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>() 
                    << std::endl;
          
          std::cout << "  Device endianess:\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_ENDIAN_LITTLE>() ? "Little" : "Big") 
                    << std::endl;
          
          std::cout << "  Available:\t\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_AVAILABLE>() ? "Yes" : "No")
                    << std::endl;
   
          std::cout << "  Compiler available:\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_COMPILER_AVAILABLE>() ? "Yes" : "No")
                    << std::endl;
          
          std::cout << "  Execution capabilities:\t\t\t\t " << std::endl;
          std::cout << "    Execute OpenCL kernels:\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_EXECUTION_CAPABILITIES>() & 
                        CL_EXEC_KERNEL ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Execute native function:\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_EXECUTION_CAPABILITIES>() & 
                        CL_EXEC_NATIVE_KERNEL ? "Yes" : "No")
                    << std::endl;
          
          std::cout << "  Queue properties:\t\t\t\t " << std::endl;
          std::cout << "    Out-of-Order:\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_QUEUE_PROPERTIES>() & 
                        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE ? "Yes" : "No")
                    << std::endl;
          std::cout << "    Profiling :\t\t\t\t\t " 
                    << ((*i).getInfo<CL_DEVICE_QUEUE_PROPERTIES>() & 
                        CL_QUEUE_PROFILING_ENABLE ? "Yes" : "No")
                    << std::endl;
          
          
          std::cout << "  Platform ID:\t\t\t\t\t " 
                << (*i).getInfo<CL_DEVICE_PLATFORM>()
                    << std::endl;
          
          std::cout << "  Name:\t\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_NAME>().c_str()
                    << std::endl;
          
          std::cout << "  Vendor:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_VENDOR>().c_str()
                    << std::endl;
          
          std::cout << "  Driver version:\t\t\t\t " 
                    << (*i).getInfo<CL_DRIVER_VERSION>().c_str()
                    << std::endl;
          
          std::cout << "  Profile:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_PROFILE>().c_str()
                    << std::endl;
          
          std::cout << "  Version:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_VERSION>().c_str()
                    << std::endl;

          std::cout << "  Extensions:\t\t\t\t\t " 
                    << (*i).getInfo<CL_DEVICE_EXTENSIONS>().c_str()
                    << std::endl;
      }
      std::cout << std::endl << std::endl;
  }

  return SDK_SUCCESS;
}



int 
IntegrationTest::genBinaryImage()
{
    cl_int err = CL_SUCCESS;
    cl::Program program;

    /*
     * Have a look at the available platforms and pick either
     * the AMD one if available or a reasonable default.
     */
    err = cl::Platform::get(&platforms);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Platform::get() failed."))
    {
        return SDK_FAILURE;
    }

    std::vector<cl::Platform>::iterator i;
    if(platforms.size() > 0)
    {
        for(i = platforms.begin(); i != platforms.end(); ++i)
        {
            if(!strcmp((*i).getInfo<CL_PLATFORM_VENDOR>().c_str(), 
                "Advanced Micro Devices, Inc."))
            {
                break;
            }
        }
    }

    cl_context_properties cps[5] = 
    {
        CL_CONTEXT_PLATFORM, 
        (cl_context_properties)(*i)(), 
        CL_CONTEXT_OFFLINE_DEVICES_AMD,
        (cl_context_properties)1,
        0
    };

    context = cl::Context(CL_DEVICE_TYPE_ALL, cps, NULL, NULL, &err);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Context::Context() failed."))
    {
        return SDK_FAILURE;
    }


    devices = context.getInfo<CL_CONTEXT_DEVICES>();
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Context::getInfo() failed."))
    {
        return SDK_FAILURE;
    }

    std::cout << "Platform :" << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
    size_t numDevices = devices.size();

    if (numDevices == 0) 
    {
        std::cerr << "No device available\n";
        return SDK_FAILURE;
    }

    /* create a CL program using the kernel source */
    streamsdk::SDKFile kernelFile;
    std::string kernelPath = sampleCommon->getPath();
    kernelPath.append("Kernel_ObtainSynapticEvents.cl");
    if(!kernelFile.open(kernelPath.c_str()))
    {
        std::cout << "Failed to load kernel file : " << kernelPath << std::endl;
        return SDK_FAILURE;
    }

    cl::Program::Sources programSource(1, 
        std::make_pair(kernelFile.source().data(), 
        kernelFile.source().size()));
    
    program = cl::Program(context, programSource, &err);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Program::program() failed."))
    {
        return SDK_FAILURE;
    }

    std::string flagsStr = std::string("");

    // Get additional options
    if(isComplierFlagsSpecified())
    {
        streamsdk::SDKFile flagsFile;
        std::string flagsPath = sampleCommon->getPath();
        flagsPath.append(flags.c_str());
        if(!flagsFile.open(flagsPath.c_str()))
        {
            std::cout << "Failed to load flags file: " << flagsPath << std::endl;
            return SDK_FAILURE;
        }
        flagsFile.replaceNewlineWithSpaces();
        const char * flags = flagsFile.source().c_str();
        flagsStr.append(flags);
    }

    if(flagsStr.size() != 0)
        std::cout << "Build Options are : " << flagsStr.c_str() << std::endl;


    err = program.build(devices, flagsStr.c_str());
    sampleCommon->checkVal(err, CL_SUCCESS, "Program::build() failed.");

    std::cout << "Number of devices found : " << numDevices << "\n\n";

    std::vector<size_t> binarySizes = program.getInfo<CL_PROGRAM_BINARY_SIZES>(&err);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Program::getInfo(CL_PROGRAM_BINARY_SIZES) failed."))
    {
        return SDK_FAILURE;
    }

    std::vector<char*> binaries = program.getInfo<CL_PROGRAM_BINARIES>(&err);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Program::getInfo(CL_PROGRAM_BINARIES) failed."))
    {
        return SDK_FAILURE;
    }   

    /* dump out each binary into its own separate file. */
    std::vector<cl::Device>::iterator deviceIndx = devices.begin();
    for(size_t i = 0; i < numDevices; i++, deviceIndx++)
    {
        char fileName[100];
        sprintf(fileName, "%s.%d", dumpBinary.c_str(), (int)i);
        if(binarySizes[i] != 0)
        {
            std::string deviceName = (*deviceIndx).getInfo<CL_DEVICE_NAME>(&err);
            if(!sampleCommon->checkVal(err,
                                       CL_SUCCESS,
                                       "Device::getInfo(CL_DEVICE_NAME) failed."))
            {
                return SDK_FAILURE;
            }

            printf( "%s binary kernel: %s\n", deviceName.c_str(), fileName);
            streamsdk::SDKFile BinaryFile;
            if(!BinaryFile.writeBinaryToFile(fileName, 
                                             binaries[i], 
                                             binarySizes[i]))
            {
                std::cout << "Failed to load kernel file : " << fileName << std::endl;
                return SDK_FAILURE;
            }
        }
        else
        {
            printf("Skipping %s since there is no binary data to write!\n",
                    fileName);
        }
    }

    return SDK_SUCCESS;
}



int 
IntegrationTest::createKernel
(
  cl::Kernel    &kernel,
  std::string   kernelFileName,
  const char    *kernelName,
  std::string   flagsStr,
  size_t        blockSizeX,
  size_t        blockSizeY
){
    cl_int err = CL_SUCCESS;
    cl::Program program;

#if (ERROR_TRY_CATCH_ENABLE)
  try
  {
#endif

    /* create a CL program using the kernel source */

    streamsdk::SDKFile kernelFile;
    std::string kernelPath = sampleCommon->getPath();
    streamsdk::SDKFile defFile;
    std::string defPath = sampleCommon->getPath();

    if(isLoadBinaryEnabled())
    {
        kernelPath.append(loadBinary.c_str());
        if(!kernelFile.readBinaryFromFile(kernelPath.c_str()))
        {
          std::cout << "Failed to load kernel file : " << kernelPath << std::endl;
          return SDK_FAILURE;
        }
        cl::Program::Binaries programBinary(1,std::make_pair(
                                              (const void*)kernelFile.source().data(), 
                                              kernelFile.source().size()));
        
        program = cl::Program(context, device, programBinary, NULL, &err);
        if(!sampleCommon->checkVal(
            err,
            CL_SUCCESS,
            "Program::Program(Binary) failed."))
            return SDK_FAILURE;

    }
    else
    {
      kernelPath.append(kernelFileName);
      if(!kernelFile.open(kernelPath.c_str()))
      {
        std::cout << "Failed to load kernel file : " << kernelPath << std::endl;
        return SDK_FAILURE;
      }
      
      defPath.append("Definitions.h");
      if(!defFile.open(defPath.c_str()))
      {
        std::cout << "Failed to load kernel file : " << defPath << std::endl;
        return SDK_FAILURE;
      }

      /*To nearest 4B chunk*/
      size_t sourceCodeSize = ((defFile.source().size() + 
        kernelFile.source().size())/sizeof(int)+1)*sizeof(int);
      /*size_t sourceCodeSize = defFile.source().size() + kernelFile.source().size();*/
      std::cout << "Source code size (" << kernelName << "): " << ((float)sourceCodeSize)/1024.0 
        << " KB" << std::endl;
      char *sourceCode = (char *) calloc(sourceCodeSize, sizeof(char));
      
      sourceCode[0] = '\0';
      strcat ( sourceCode, defFile.source().data() );
      strcat ( sourceCode, kernelFile.source().data() );
      //std::cout << "Source code:\n" << sourceCode << std::endl;
      cl::Program::Sources programSource(1, std::make_pair(sourceCode, sourceCodeSize));

      program = cl::Program(context, programSource, &err);
      if(!sampleCommon->checkVal(
          err,
          CL_SUCCESS,
          "Program::Program(Source) failed."))
          return SDK_FAILURE;
    }

    streamsdk::SDKFile flagsFile;
    std::string flagsPath = OCL_COMPILER_OPTIONS_FILE_NAME;
    if(!flagsFile.open(flagsPath.c_str()))
    {
      std::cout << "Failed to load flags file: " << flagsPath << std::endl;
      return SDK_FAILURE;
    }
    flagsFile.replaceNewlineWithSpaces();
    const char *flags = flagsFile.source().c_str();
    flagsStr.append(flags);

    std::cout << "Building for devices: ";
    for (std::vector<cl::Device>::iterator d = device.begin(); d != device.end(); d++) 
    {
      std::cout << (*d).getInfo<CL_DEVICE_NAME>() << ",";
    }
    std::cout <<std::endl;
    
    if(flagsStr.size() != 0)
        std::cout << "Build Options (" << kernelName << "): " << flagsStr.c_str() << std::endl;

    err = program.build(device, flagsStr.c_str());
    if(err != CL_SUCCESS)
    {
      if(err == CL_BUILD_PROGRAM_FAILURE)
      {
          std::string str = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device.front());

          std::cout << " \n\t\t\tBUILD LOG\n";
          std::cout << " ************************************************\n";
          std::cout << str << std::endl;
          std::cout << " ************************************************\n";
      }
    }

    if(!sampleCommon->checkVal(
        err,
        CL_SUCCESS,
        "Program::build() failed."))
        return SDK_FAILURE;
    
    /* Create kernel */
    kernel = cl::Kernel(program, kernelName, &err);
    if(!sampleCommon->checkVal(err, CL_SUCCESS, "Kernel::Kernel() failed for"))
    {
        return SDK_FAILURE;
    }

    /* Check group size against group size returned by kernel */
    size_t kernelWorkGroupSize = 
      kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device.front(), &err);
    if(!sampleCommon->checkVal(err, CL_SUCCESS, "Kernel::getWorkGroupInfo()  failed."))
    {
      return SDK_FAILURE;
    }

    if((blockSizeX * blockSizeY) > kernelWorkGroupSize)
    {
      std::cout << "Out of Resources!" << std::endl;
      std::cout << "Group Size specified : "
                << blockSizeX * blockSizeY << std::endl;
      std::cout << "Max Group Size supported on the kernel : "
                << kernelWorkGroupSize << std::endl;
                
      return SDK_FAILURE;
    }
#if (ERROR_TRY_CATCH_ENABLE)
  }
  catch (cl::Error err) 
  {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    return SDK_FAILURE;
  }
#endif

  return SDK_SUCCESS;
}



/*Find target device from the list of target devices*/
bool 
IntegrationTest::findTargetDevice
(
  vector<cl::Device>                  platformDevices,
  const char                          *targetDevices,
  std::vector<cl::Device>::iterator   *d
){
  bool found = false;
  char *str = (char *) calloc(0xFFFF, sizeof(char));
  strcpy(str, targetDevices);
  char *pch;
  pch = strtok (str, ",");
  
  while (pch != NULL)
  {
    for(*d = platformDevices.begin(); *d != platformDevices.end(); ++(*d)) 
    {
      std::string deviceName = (*(*d)).getInfo<CL_DEVICE_NAME>();
      
      if(strcmp(deviceName.c_str(), pch) == 0)
      {
        found = true; break;
      }
    }
    if(found){break;}
    pch = strtok (NULL, ",");
  }
  
  free(str);
  return found;
}



int 
IntegrationTest::setupCL()
{
    cl_int err = CL_SUCCESS;
    cl_device_type dType;
    cl_uint myDeviceId;
    bool found = false;

#if (ERROR_TRY_CATCH_ENABLE)
  try
  {
#endif
  
    err = cl::Platform::get(&platforms);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Platform::get() failed."))
    {
        return SDK_FAILURE;
    }

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
      std::cout << "Unable to find target platform with vendor " << 
        (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
      return SDK_FAILURE;
    }
    
    (*i).getDevices(CL_DEVICE_TYPE_ALL, &devices);
 
    std::vector<cl::Device>::iterator d;
    found = findTargetDevice(devices, TARGET_DEVICE_NAME, &d);
    
    if(!found)
    {
      std::cout << "Unable to find target devices " << TARGET_DEVICE_NAME 
        << " on platform from vendor "  << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
      return SDK_FAILURE;
    }
    
    dType = (*d).getInfo<CL_DEVICE_TYPE>();
    myDeviceId = (*d).getInfo<CL_DEVICE_VENDOR_ID>();

    cl_context_properties cps[3] = 
    { 
        CL_CONTEXT_PLATFORM, 
        (cl_context_properties)(*i)(),
        0 
    };

    context = cl::Context(dType, cps, NULL, NULL, &err);
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Context::Context() failed."))
    {
        return SDK_FAILURE;
    }

    devices = context.getInfo<CL_CONTEXT_DEVICES>();
    if(!sampleCommon->checkVal(err, 
        CL_SUCCESS,
        "Context::getInfo() failed."))
    {
        return SDK_FAILURE;
    }
    
    int deviceCount = (int)devices.size();

    found = findTargetDevice(devices, TARGET_DEVICE_NAME, &d);
    
    if(!found)
    {
      std::cout << "Unable to find target devices " << TARGET_DEVICE_NAME 
        << " within created context on platform from vendor "  
        << (*i).getInfo<CL_PLATFORM_VENDOR>().c_str() << "\n";
      return SDK_FAILURE;
    }
    
    device.push_back(*d);
    
    /*Validate mem allocations against device-specific restrictions for all registered kernels*/
    {
      bool print = 1;
      set<std::string>::iterator kernel_name;
      cl_uint memBaseAllignBytes = (*d).getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>()/8;
      cl_uint minDataTypeAlignSize = (*d).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>();
      cl_ulong memMaxAllocactionSize = (*d).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
      cl_ulong maxGlobalMemSize = (*d).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
      cl_ulong maxConstantMemSize = (*d).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>();
      cl_ulong maxLocalMemSize = (*d).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();

      print ? std::cout << "Memory Allocations:" << std::endl, true : false;
      
      for
      (
        kernel_name = kernelStats.kernelNames.begin();
        kernel_name != kernelStats.kernelNames.end(); 
        ++kernel_name
      ){
        map<std::string, cl_uint>::iterator m;
        
        print ? std::cout << " Kernel " << *kernel_name << ":" << std::endl, true : false;

        /*Validate GM allocations*/
        print ? std::cout << "  Global Memory (global scope):" << std::endl, true : false;
        map<std::string, cl_uint> gmSizes = kernelStats.gmSizes[*kernel_name];
        cl_ulong gmAllSizeBytes = 0, gmMaxSizeBytes = 0;
        for
        (
          m = gmSizes.begin(); 
          m != gmSizes.end(); 
          ++m
        ){
          std::string name = (m->first);
          cl_uint size = (m->second);
          
          size = (size/minDataTypeAlignSize + 1)*minDataTypeAlignSize;
          gmAllSizeBytes += size;
          if(gmMaxSizeBytes < size){gmMaxSizeBytes = size;}
          
          print ? std::cout << "   " << name << ": " << ((float)size)/(1024.0*1024.0) << " MB" 
            << std::endl, true : false;
          
          if(size > memMaxAllocactionSize)
          {
            std::cout << "Memory object in kernel " << (*kernel_name) 
            << " with size identifier " << name << " and size " 
            << ((float)size)/(1024.0*1024.0) << " MB exceeds CL_DEVICE_MAX_MEM_ALLOC_SIZE, " 
            << ((float)memMaxAllocactionSize)/(1024.0*1024.0) << "\n";
            if(*kernel_name == KERNEL_ALL)
            {
              LOG("Max GM:" << (((float)size)/(1024.0*1024.0)), 1);
            }
            return SDK_FAILURE;
          }
        }
        print ? std::cout << "   TOTAL: " << ((float)gmAllSizeBytes)/(1024.0*1024.0) << " MB" 
          << std::endl, true : false;
          
        if(*kernel_name == KERNEL_ALL)
        {
          LOG("Total GM:" << (((float)gmAllSizeBytes)/(1024.0*1024.0)), 1);
          LOG("Max GM:" << (((float)gmMaxSizeBytes)/(1024.0*1024.0)), 1);
        }

        if(gmAllSizeBytes > maxGlobalMemSize)
        {
          std::cout << "Total allocation for global memory in kernel " << (*kernel_name) 
          << ", " << ((float)gmAllSizeBytes)/(1024.0*1024.0) 
          << " MB, exceeds CL_DEVICE_GLOBAL_MEM_SIZE, " 
          << ((float)maxGlobalMemSize)/(1024.0*1024.0) << "\n";
          return SDK_FAILURE;
        }
        
        /*Validate CM allocations*/
        print ? std::cout << "  Constant Memory (global scope):" << std::endl, true : false;
        map<std::string, cl_uint> cmSizes = kernelStats.cmSizes[*kernel_name];
        cl_ulong cmAllSizeBytes = 0;
        for
        (
          m = cmSizes.begin(); 
          m != cmSizes.end(); 
          ++m
        ){
          std::string name = (m->first);
          cl_uint size = (m->second);
          
          size = (size/minDataTypeAlignSize + 1)*minDataTypeAlignSize;
          cmAllSizeBytes += size;
          
          print ? std::cout << "   " << name << ": " << ((float)size)/(1024.0) << " KB" 
            << std::endl, true : false;
            
        }
        print ? std::cout << "   TOTAL: " << ((float)cmAllSizeBytes)/(1024.0) << " KB" 
          << std::endl, true : false;
        if(cmAllSizeBytes > maxConstantMemSize)
        {
          std::cout << "Total allocation for constant memory in kernel " << (*kernel_name) 
          << ", " << ((float)cmAllSizeBytes)/(1024.0) 
          << " KB, exceeds CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, " 
          << ((float)maxConstantMemSize)/(1024.0) << " KB\n";
          return SDK_FAILURE;
        }
        
        /*Validate LM allocations*/
        print ? std::cout << "  Local Memory (WG scope):" << std::endl, true : false;
        map<std::string, cl_uint> lmSizes = kernelStats.lmSizes[*kernel_name];
        cl_ulong lmAllSizeBytes = 0;
        for
        (
          m = lmSizes.begin(); 
          m != lmSizes.end(); 
          ++m
        ){
          std::string name = (m->first);
          cl_uint size = (m->second);
          
          size = (size/minDataTypeAlignSize + 1)*minDataTypeAlignSize;
          lmAllSizeBytes += size;
          
          print ? std::cout << "   " << name << ": " << ((float)size)/(1024.0) << " KB" 
            << std::endl, true : false;
            
        }
        print ? std::cout << "   TOTAL: " << ((float)lmAllSizeBytes)/(1024.0) << " KB" 
          << std::endl, true : false;
        if(lmAllSizeBytes > maxLocalMemSize)
        {
          std::cout << "Total allocation for local memory in kernel " << (*kernel_name) 
          << ", " << ((float)lmAllSizeBytes)/(1024.0) 
          << " KB, exceeds CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, " 
          << ((float)maxLocalMemSize)/(1024.0) << " KB\n";
          /*return SDK_FAILURE;*/
        }
      }
    }

    
    std::cout << "Creating command queue for device " << (*d).getInfo<CL_DEVICE_NAME>() 
      << std::endl;
    commandQueue = cl::CommandQueue(context, *d, 0, &err);
    if(!sampleCommon->checkVal(err, CL_SUCCESS, "CommandQueue::CommandQueue() failed."))
    {
        return SDK_FAILURE;
    }

    /*Create and initialize memory objects and kernels*/
    
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
  CREATE_BUFFER(CL_MEM_READ_WRITE, dataHistogramGroupEventsTikBuffer, 
    dataHistogramGroupEventsTikSizeBytes);
#endif

#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
  CREATE_BUFFER(CL_MEM_READ_WRITE, dataUnsortedEventCountsBuffer, dataUnsortedEventCountsSizeBytes);
  CREATE_BUFFER(CL_MEM_READ_WRITE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes);
  CREATE_BUFFER(CL_MEM_READ_WRITE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes);
  CREATE_BUFFER(CL_MEM_READ_WRITE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes);
#endif

#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
  CREATE_BUFFER(CL_MEM_READ_WRITE, dataHistogramBuffer, dataHistogramSizeBytes);
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

#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
    {
#if (SCAN_DEBUG_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataScanErrorBuffer, dataScanErrorSizeBytes);
#endif
    }
#endif

#if EXPAND_EVENTS_ENABLE || UPDATE_NEURONS_ENABLE_V00
    {
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes);
    }
#endif
    
#if EXPAND_EVENTS_ENABLE
    {
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataExpandEventsDebugHostBuffer, 
      dataExpandEventsDebugHostSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataExpandEventsDebugDeviceBuffer, 
      dataExpandEventsDebugDeviceSizeBytes);
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    CREATE_BUFFER(CL_MEM_READ_WRITE, dataExpandEventsErrorBuffer, dataExpandEventsErrorSizeBytes);
#endif
    CREATE_BUFFER(CL_MEM_READ_ONLY, dataSynapseTargetsBuffer, dataSynapseTargetsSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_ONLY, dataSynapseDelaysBuffer, dataSynapseDelaysSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_ONLY, dataSynapseWeightsBuffer, dataSynapseWeightsSizeBytes);
    CREATE_BUFFER(CL_MEM_READ_ONLY, dataSynapsePointerBuffer, dataSynapsePointerSizeBytes);
    
    createKernel
    (
      kernelExpandEvents,
      EXPAND_EVENTS_KERNEL_FILE_NAME,
      EXPAND_EVENTS_KERNEL_NAME,
      "",
      blockSizeX_KernelExpandEvents,
      blockSizeY_KernelExpandEvents
    );
    }
#endif

#if SCAN_ENABLE_V00
    {
    createKernel
    (
      kernelScanHistogramV00,
      SCAN_KERNEL_FILE_NAME,
      SCAN_KERNEL_NAME,
      "-D SCAN_DEVICE_V00",
      blockSizeX_kernelScanHistogramV00,
      blockSizeY_kernelScanHistogramV00
    );
    }
#endif

#if GROUP_EVENTS_ENABLE_V00
    {
    createKernel
    (
      kernelGroupEventsV00,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D GROUP_EVENTS_DEVICE_V00",
      blockSizeX_kernelGroupEventsV00,
      blockSizeY_kernelGroupEventsV00
    );
    }
#endif

#if SCAN_ENABLE_V01
    {
    createKernel
    (
      kernelScanHistogramV01,
      SCAN_KERNEL_FILE_NAME,
      SCAN_KERNEL_NAME,
      "-D SCAN_DEVICE_V01",
      blockSizeX_kernelScanHistogramV01,
      blockSizeY_kernelScanHistogramV01
    );
    }
#endif

#if GROUP_EVENTS_ENABLE_V01 
    {
    createKernel
    (
      kernelGroupEventsV01,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D GROUP_EVENTS_DEVICE_V01",
      blockSizeX_kernelGroupEventsV01,
      blockSizeY_kernelGroupEventsV01
    );
    }
#endif

#if GROUP_EVENTS_ENABLE_V02
    {
    createKernel
    (
      kernelGroupEventsV02,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D GROUP_EVENTS_DEVICE_V02",
      blockSizeX_kernelGroupEventsV02,
      blockSizeY_kernelGroupEventsV02
    );
    }
#endif

#if GROUP_EVENTS_ENABLE_V03
    {
    createKernel
    (
      kernelGroupEventsV03,
      GROUP_EVENTS_KERNEL_FILE_NAME,
      GROUP_EVENTS_KERNEL_NAME,
      "-D GROUP_EVENTS_DEVICE_V03",
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
      kernelMakeEventPtrs,
      MAKE_EVENT_PTRS_KERNEL_FILE_NAME,
      MAKE_EVENT_PTRS_KERNEL_NAME,
      "",
      blockSizeX_kernelMakeEventPtrs,
      blockSizeY_kernelMakeEventPtrs
    );
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
    createKernel
    (
      kernelGlueEventPtrs,
      GLUE_EVENT_PTRS_KERNEL_FILE_NAME,
      GLUE_EVENT_PTRS_KERNEL_NAME,
      "",
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
      kernelUpdateNeuronsV00,
      UPDATE_NEURONS_KERNEL_FILE_NAME,
      UPDATE_NEURONS_KERNEL_NAME,
#if COMPILER_FLAGS_OPTIMIZE_ENABLE
      "-D UPDATE_NEURONS_DEVICE_V00 -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros "
      "-cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math",
#else
      /*"-D UPDATE_NEURONS_DEVICE_V00 -cl-fp32-correctly-rounded-divide-sqrt -cl-opt-disable",*/
      /*"-D UPDATE_NEURONS_DEVICE_V00",*/
      "-D UPDATE_NEURONS_DEVICE_V00",
#endif
      blockSizeX_kernelUpdateNeuronsV00,
      blockSizeY_kernelUpdateNeuronsV00
    );
    
    createKernel
    (
      kernelUpdateSpikedNeuronsV00,
      UPDATE_NEURONS_KERNEL_FILE_NAME,
      UPDATE_NEURONS_SPIKED_KERNEL_NAME,
#if COMPILER_FLAGS_OPTIMIZE_ENABLE
      "-D UPDATE_NEURONS_DEVICE_V00 -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros "
      "-cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math",
#else
      /*"-D UPDATE_NEURONS_DEVICE_V00 -cl-fp32-correctly-rounded-divide-sqrt -cl-opt-disable",*/
      /*"-D UPDATE_NEURONS_DEVICE_V00 -cl-fp32-correctly-rounded-divide-sqrt",*/
      "-D UPDATE_NEURONS_DEVICE_V00 -cl-fp32-correctly-rounded-divide-sqrt",
#endif
      blockSizeX_kernelUpdateNeuronsV00,
      blockSizeY_kernelUpdateNeuronsV00
    );
    }
#endif

#if (ERROR_TRY_CATCH_ENABLE)
  }
  catch (cl::Error err) 
  {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    return SDK_FAILURE;
  }
#endif

  return SDK_SUCCESS;
}



double
IntegrationTest::timeStampNs()
{
#ifdef WIN32
  QueryPerformanceCounter((LARGE_INTEGER *)&current);
  return ((double)current / (double)performanceFrequency);
#else
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return (double) tp.tv_sec * (1000ULL * 1000ULL * 1000ULL) + (double) tp.tv_nsec;
#endif
}



int 
IntegrationTest::runCLKernels()
{
  cl_int status;
  cl_int eventStatus = CL_QUEUED;
  cl::Event writeEvt = NULL;

#if (ERROR_TRY_CATCH_ENABLE)
  try
  {
#endif

/*
  Initialize network
*/
#if EXPAND_EVENTS_ENABLE
    if(initializeExpandEventsData() != SDK_SUCCESS)
    {
      std::cout << "Failed initializeExpandEventsData" << std::endl;
      return SDK_FAILURE;
    }
#if PREINITIALIZE_NETWORK_STATE
    if(initializeSpikeData(PREINITIALIZE_NETWORK_MIN_SPIKE_PERCENT, 
      PREINITIALIZE_NETWORK_MAX_SPIKE_PERCENT) != SDK_SUCCESS)
#else
    if(initializeSpikeData(EXPAND_EVENTS_MIN_MAX_SPIKE_PERCENT) != SDK_SUCCESS)
#endif
    {
      std::cout << "Failed initializeSpikeData" << std::endl;
      return SDK_FAILURE;
    }
#endif
#if PREINITIALIZE_NETWORK_STATE && GROUP_EVENTS_ENABLE_V00 && EXPAND_EVENTS_ENABLE
    {
      int res = initializeEventBuffers
      (
        GROUP_EVENTS_TIME_SLOTS,
        GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
        GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE,
        GROUP_EVENTS_TOTAL_NEURONS,
        PREINITIALIZE_NETWORK_TIME_SLOT_DELTA,
        GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00,
        GROUP_EVENTS_HISTOGRAM_BIN_MASK,
        GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
        GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
        GROUP_EVENTS_MIN_DELAY,
        1.0,
        PREINITIALIZE_NETWORK_TIME_SLOT_DELTA_DEVIATION,
        PREINITIALIZE_NETWORK_PERCENT_INHIBITORY,
        dataUnsortedEventCounts,
        dataUnsortedEventTargets,
        dataUnsortedEventDelays,
        dataUnsortedEventWeights,
        dataHistogram
      );
      
      if(res != SDK_SUCCESS)
      {
        std::cout << "Failed initializeEventBuffers" << std::endl;
        return SDK_FAILURE;
      }
      
      memcpy(dataUnsortedEventCountsVerify, dataUnsortedEventCounts, 
        dataUnsortedEventCountsSizeBytes);
      memcpy(dataUnsortedEventsTargetsVerify, dataUnsortedEventTargets, 
        dataUnsortedEventTargetsSizeBytes);
      memcpy(dataUnsortedEventsDelaysVerify, dataUnsortedEventDelays, 
        dataUnsortedEventDelaysSizeBytes);
      memcpy(dataUnsortedEventsWeightsVerify, dataUnsortedEventWeights, 
        dataUnsortedEventWeightsSizeBytes);
      memcpy(dataHistogramVerify, dataHistogram, 
        dataHistogramSizeBytes);
    }
#endif
#if MAKE_EVENT_PTRS_ENABLE
    if(initializeDataForKernelMakeEventPtrs(0, 0) != SDK_SUCCESS)
    {
      std::cout << "Failed initializeDataForKernelMakeEventPtrs" << std::endl; 
    }
#endif
#if UPDATE_NEURONS_ENABLE_V00
#if PREINITIALIZE_NETWORK_STATE
    if(initializeDataForKernelUpdateNeurons(0, 1, 1, 0, 
      PREINITIALIZE_NETWORK_NEURON_VARIABLES_SAMPLE_FILE) != SDK_SUCCESS)
#else
    if(initializeDataForKernelUpdateNeurons(0, 1, 1, 0, "") != SDK_SUCCESS)
#endif
    {
      std::cerr << "ERROR, runCLKernels: Failed initializeDataForKernelUpdateNeurons" << std::endl;
      return SDK_FAILURE;
    }
#if (LOG_MODEL_VARIABLES)
    traceFile = new std::ofstream(LOG_MODEL_VARIABLES_FILE_NAME);
    dataToTraceFile = new std::stringstream("", std::ios::out | std::ios::app);

    if(!(*traceFile).is_open())
    {
      std::cerr << "ERROR, runCLKernels: Unable to open trace file";
      return SDK_FAILURE;
    }

    *dataToTraceFile << startTimeStamp;
    *dataToTraceFile << LOG_MODEL_VARIABLES_FILE_HEADER << std::endl;
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
        GROUP_EVENTS_TOTAL_NEURONS,
        REFERENCE_EVENT_QUEUE_SIZE,
        dataUnsortedEventCounts,
        dataUnsortedEventTargets,
        dataUnsortedEventWeights,
        dataUnsortedEventDelays,
        nrn_ps
      );
      
      if(res != SDK_SUCCESS)
      {
        std::cout << "Failed injectUnsortedEvents" << std::endl;
        return SDK_FAILURE;
      }
    }
#endif
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataHistogramGroupEventsTikBuffer, 
    dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif

#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataUnsortedEventCountsBuffer, dataUnsortedEventCountsSizeBytes, 
    dataUnsortedEventCounts);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
    dataUnsortedEventTargets);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
    dataUnsortedEventDelays);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
    dataUnsortedEventWeights);
#endif

#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
#endif

#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataGroupEventsTikBuffer, 
    dataGroupEventsTikSizeBytes, dataGroupEventsTik);
  }
#endif

#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
  {
#if (GROUP_EVENTS_DEBUG_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
    dataDebugHostGroupEvents);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataDebugDeviceGroupEventsBuffer, 
    dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
    dataErrorGroupEvents);
#endif
  }
#endif

#if EXPAND_EVENTS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
    dataSpikePackets);
#endif

#if EXPAND_EVENTS_ENABLE
  {
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataExpandEventsDebugHostBuffer, 
    dataExpandEventsDebugHostSizeBytes, dataExpandEventsDebugHost);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataExpandEventsDebugDeviceBuffer, 
    dataExpandEventsDebugDeviceSizeBytes, dataExpandEventsDebugDevice);
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataExpandEventsErrorBuffer, dataExpandEventsErrorSizeBytes, 
    dataExpandEventsError);
#endif
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataSynapseTargetsBuffer, dataSynapseTargetsSizeBytes, 
    dataSynapseTargets);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataSynapseDelaysBuffer, dataSynapseDelaysSizeBytes, 
    dataSynapseDelays);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataSynapseWeightsBuffer, dataSynapseWeightsSizeBytes, 
    dataSynapseWeights);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataSynapsePointerBuffer, dataSynapsePointerSizeBytes, 
    dataSynapsePointer);
  }
#endif

#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
  {
#if (SCAN_DEBUG_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
    dataScanDebugHost);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
    dataScanDebugDevice);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanErrorBuffer, dataScanErrorSizeBytes, dataScanError);
#endif
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataHistogramGroupEventsTokBuffer, 
    dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
  }
#endif

#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataGroupEventsTokBuffer, 
    dataGroupEventsTokSizeBytes, dataGroupEventsTok);
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  {
  ENQUEUE_WRITE_BUFFER(CL_FALSE, dataMakeEventPtrsStructBuffer, 
    dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
  }
#endif

#if UPDATE_NEURONS_ENABLE_V00
  {
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
  ENQUEUE_WRITE_BUFFER(CL_FALSE, psToleranceBuffer, psToleranceSizeBytes, psTolerance);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUpdateNeuronsErrorBuffer, dataUpdateNeuronsErrorSizeBytes, 
      dataUpdateNeuronsError);
#endif
  ENQUEUE_WRITE_BUFFER(CL_FALSE, constantCoefficientsBuffer, constantCoefficientsSizeBytes, 
    constantCoefficients);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, modelVariablesBuffer, modelVariablesSizeBytes, 
    modelVariables);
  ENQUEUE_WRITE_BUFFER(CL_FALSE, modelParametersBuffer, modelParametersSizeBytes, 
    modelParameters);
  }
#endif

#if MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataMakeEventPtrsErrorBuffer, dataMakeEventPtrsErrorSizeBytes, 
      dataMakeEventPtrsError);
#endif
#endif

  /* Set arguments to the kernels */
  
#if EXPAND_EVENTS_ENABLE
  cl_uint argNumExpandEvents = 0;
  {
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelExpandEvents, dataExpandEventsDebugHostBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataExpandEventsDebugDeviceBuffer, argNumExpandEvents++);
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelExpandEvents, dataExpandEventsErrorBuffer, argNumExpandEvents++);
#endif
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  SET_KERNEL_ARG(kernelExpandEvents, dataHistogramBuffer, argNumExpandEvents++);
#endif
  SET_KERNEL_ARG(kernelExpandEvents, dataSpikePacketsBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataUnsortedEventCountsBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataUnsortedEventTargetsBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataUnsortedEventDelaysBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataUnsortedEventWeightsBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataSynapseTargetsBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataSynapseDelaysBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataSynapseWeightsBuffer, argNumExpandEvents++);
  SET_KERNEL_ARG(kernelExpandEvents, dataSynapsePointerBuffer, argNumExpandEvents++);
  }
#endif

#if SCAN_ENABLE_V00
  cl_uint argNumScan00 = 0;
  {
#if (SCAN_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelScanHistogramV00, dataScanDebugHostBuffer, argNumScan00++);
  SET_KERNEL_ARG(kernelScanHistogramV00, dataScanDebugDeviceBuffer, argNumScan00++);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelScanHistogramV00, dataScanErrorBuffer, argNumScan00++);
#endif
  SET_KERNEL_ARG(kernelScanHistogramV00, dataHistogramBuffer, argNumScan00++);
  }
#endif

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
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  SET_KERNEL_ARG(kernelGroupEventsV00, dataHistogramGroupEventsTikBuffer, argNumGrupEventsV00++);
#endif
  SET_KERNEL_ARG(kernelGroupEventsV00, dataUnsortedEventCountsBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, dataUnsortedEventDelaysBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, dataGroupEventsTikBuffer, argNumGrupEventsV00++);
  SET_KERNEL_ARG(kernelGroupEventsV00, dataHistogramBuffer, argNumGrupEventsV00++);
  }
#endif

#if SCAN_ENABLE_V01
  cl_uint argNumScan01 = 0;
  {
#if (SCAN_DEBUG_ENABLE)
  SET_KERNEL_ARG(kernelScanHistogramV01, dataScanDebugHostBuffer, argNumScan01++);
  SET_KERNEL_ARG(kernelScanHistogramV01, dataScanDebugDeviceBuffer, argNumScan01++);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  SET_KERNEL_ARG(kernelScanHistogramV01, dataScanErrorBuffer, argNumScan01++);
#endif
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
  SET_KERNEL_ARG(kernelGroupEventsV02, dataUnsortedEventTargetsBuffer, argNumGrupEventsV02++);
  SET_KERNEL_ARG(kernelGroupEventsV02, dataUnsortedEventDelaysBuffer, argNumGrupEventsV02++);
  SET_KERNEL_ARG(kernelGroupEventsV02, dataUnsortedEventWeightsBuffer, argNumGrupEventsV02++);
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
  SET_KERNEL_ARG(kernelGroupEventsV03, dataUnsortedEventTargetsBuffer, argNumGrupEventsV03++);
  SET_KERNEL_ARG(kernelGroupEventsV03, dataUnsortedEventDelaysBuffer, argNumGrupEventsV03++);
  SET_KERNEL_ARG(kernelGroupEventsV03, dataUnsortedEventWeightsBuffer, argNumGrupEventsV03++);
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
  SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataSpikePacketsBuffer, argNumUpdateNeuronsV00++);
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
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataSpikePacketsBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataMakeEventPtrsStructBuffer, 
    argNumUpdateSpikedNeuronsV00++);
  }
#endif

  /* 
  * Enqueue a kernel run call.
  */
  
  cl::Event ndrEvt;
  
#if EXPAND_EVENTS_ENABLE
  cl::NDRange globalThreadsExpandEvents(EXPAND_EVENTS_WG_SIZE_WI*EXPAND_EVENTS_GRID_SIZE_WG);
  cl::NDRange localThreadsExpandEvents(blockSizeX_KernelExpandEvents, blockSizeY_KernelExpandEvents);
#endif

#if SCAN_ENABLE_V00
  cl::NDRange globalThreadsScanV00(SCAN_WG_SIZE_WI*SCAN_GRID_SIZE_WG);
  cl::NDRange localThreadsScanV00(blockSizeX_kernelScanHistogramV00, 
    blockSizeY_kernelScanHistogramV00);
#endif

#if GROUP_EVENTS_ENABLE_V00
  cl::NDRange globalThreadsGroupEventsV00(GROUP_EVENTS_WG_SIZE_WI*GROUP_EVENTS_GRID_SIZE_WG);
  cl::NDRange localThreadsGroupEventsV00(blockSizeX_kernelGroupEventsV00, 
    blockSizeY_kernelGroupEventsV00);
#endif

#if SCAN_ENABLE_V01
  cl::NDRange globalThreadsScanV01(SCAN_WG_SIZE_WI*SCAN_GRID_SIZE_WG);
  cl::NDRange localThreadsScanV01(blockSizeX_kernelScanHistogramV01, 
    blockSizeY_kernelScanHistogramV01);
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

  bool verified = true;
#if (PROFILING_MODE == 1 || PROFILING_MODE == 2) && (START_PROFILING_AT_STEP > 0)
  double startAppTime = 0, endAppTime = 0;
#endif
#if SIMULATION_SNAPSHOT
  bool takeSimSnapshot = false;
#endif

  /*Iterate through steps*/
  std::cout << "\nStarting execution" << std::endl;
  for(currentTimeStep = 0; currentTimeStep < SIMULATION_TIME_STEPS; currentTimeStep++)
  {
#if SIMULATION_SNAPSHOT
    if(currentTimeStep == SIMULATION_TIME_STEPS-1){takeSimSnapshot = true;}
    else{takeSimSnapshot = false;}
#endif
#if (SIMULATION_MODE == 0 || ERROR_TRACK_ENABLE)
    currentTimeSlot = currentTimeStep%EVENT_TIME_SLOTS;
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
#if PROFILING_MODE == 1 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep == START_PROFILING_AT_STEP)
    {
      std::cout << "\nStarting device timing at step " << currentTimeStep << std::endl;
      startAppTime = timeStampNs();
    }
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
    {
    cl_uint expandEventsTimeStep = currentTimeStep;
/*Unit test initialization*/
#if !(UPDATE_NEURONS_ENABLE_V00)
    /*Reset the data with new values every resetAtSlot steps for unit test for better represenation*/
    expandEventsTimeStep = expandEventsTimeStep%EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT;
    if(!expandEventsTimeStep)
    {
    if(initializeExpandEventsData() != SDK_SUCCESS)
    {
      std::cout << "Failed initializeExpandEventsData" << std::endl; 
      verified = false; 
      break;
    }
    if(initializeSpikeData(EXPAND_EVENTS_MIN_MAX_SPIKE_PERCENT) != SDK_SUCCESS)
    {
      std::cout << "Failed initializeSpikeData" << std::endl;
      verified = false; 
      break;
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
      dataSpikePackets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSynapseTargetsBuffer, dataSynapseTargetsSizeBytes, 
      dataSynapseTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSynapseDelaysBuffer, dataSynapseDelaysSizeBytes, 
      dataSynapseDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSynapseWeightsBuffer, dataSynapseWeightsSizeBytes, 
      dataSynapseWeights);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSynapsePointerBuffer, dataSynapsePointerSizeBytes, 
      dataSynapsePointer);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventCountsBuffer, dataUnsortedEventCountsSizeBytes, 
      dataUnsortedEventCounts);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
#if EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
#endif
    }
    else
    {
      int result = SDK_SUCCESS;
      
      if(expandEventsTimeStep == 1)
      {
        result = initializeSpikeData(0.0, 0.0);
      }
      else if(expandEventsTimeStep == EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT - 1)
      {
        result = initializeSpikeData(100.0, 100.0);
      }
      else
      {
        result = initializeSpikeData(0.0, 
          (double)(100*expandEventsTimeStep/EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT));
      }
      
      if(result != SDK_SUCCESS)
      {
        std::cout << "Failed initializeSpikeData" << std::endl; 
        verified = false; 
        break;
      }
      ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
        dataSpikePackets);
    }
/*End unit test initialization*/
#elif OVERWRITE_SPIKES_UNTILL_STEP > 0
  if(currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP)
  {
    if(initializeSpikeData(OVERWRITE_SPIKES_MIN_MAX_PERCENT) != SDK_SUCCESS)
    {
      std::cout << "Failed initializeSpikeData" << std::endl; 
      verified = false; 
      break;
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
      dataSpikePackets);
    if(currentTimeStep == OVERWRITE_SPIKES_UNTILL_STEP-1)
    {
      std::cout << "\nCompleted injecting spikes at step " << currentTimeStep << std::endl;
    }
  }
#endif
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataExpandEventsDebugHostBuffer, 
      dataExpandEventsDebugHostSizeBytes, dataExpandEventsDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataExpandEventsDebugDeviceBuffer, 
      dataExpandEventsDebugDeviceSizeBytes, dataExpandEventsDebugDevice);
#endif
#if EXPAND_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
      dataSpikePackets);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSynapseTargetsBuffer, dataSynapseTargetsSizeBytes, 
      dataSynapseTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSynapseDelaysBuffer, dataSynapseDelaysSizeBytes, 
      dataSynapseDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSynapseWeightsBuffer, dataSynapseWeightsSizeBytes, 
      dataSynapseWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSynapsePointerBuffer, dataSynapsePointerSizeBytes, 
      dataSynapsePointer);
#endif
    SET_KERNEL_ARG(kernelExpandEvents, expandEventsTimeStep, argNumExpandEvents);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelExpandEvents, cl::NullRange, 
      globalThreadsExpandEvents, localThreadsExpandEvents, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelExpandEvents, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (EXPAND_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataExpandEventsDebugHostBuffer, 
      dataExpandEventsDebugHostSizeBytes, dataExpandEventsDebugHost);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataExpandEventsDebugDeviceBuffer, 
      dataExpandEventsDebugDeviceSizeBytes, dataExpandEventsDebugDevice);
#endif
#if EXPAND_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventCountsBuffer, dataUnsortedEventCountsSizeBytes, 
      dataUnsortedEventCounts);
#elif SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventCountsBuffer, dataUnsortedEventCountsSizeBytes, 
        dataUnsortedEventCounts);
    }
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataExpandEventsErrorBuffer, dataExpandEventsErrorSizeBytes, 
        dataExpandEventsError);
      if(dataExpandEventsError[0] != 0)
      {
        std::cout << "EXPAND_EVENTS_ERROR_TRACK_ENABLE: received error code from the device: " 
          << dataExpandEventsError[0] << std::endl;
        verified = false; break;
      }
    }
#endif
#if EXPAND_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
#if EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
#endif
    if(verifyKernelExpandEvents(expandEventsTimeStep, 
      (!expandEventsTimeStep)&(!PREINITIALIZE_NETWORK_STATE)) != SDK_SUCCESS)
    {
      std::cout << "Failed verifyKernelExpandEvents" << std::endl; 
      verified = false; 
      break;
    }
#endif
    }
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V00
    {
/*Unit test initialization*/
#if !(EXPAND_EVENTS_ENABLE)
    if(initializeDataForKernelScanHistogramV00(currentTimeSlot, 0) != SDK_SUCCESS)
    {
      std::cout 
      << "Failed initializeDataForKernelScanHistogramV00" << std::endl; 
      verified = false; 
      break;
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
#endif
/*End unit test initialization*/
#if (SCAN_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
      dataScanDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
      dataScanDebugDevice);
#endif
#if SCAN_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
    if(initializeDataForKernelScanHistogramV00(currentTimeSlot, 1) != SDK_SUCCESS)
    {
      std::cout 
      << "Failed initializeDataForKernelScanHistogramV00" << std::endl; 
      verified = false; 
      break;
    }
#endif
    SET_KERNEL_ARG(kernelScanHistogramV00, currentTimeStep, argNumScan00);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelScanHistogramV00, cl::NullRange, 
      globalThreadsScanV00, localThreadsScanV00, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelScanHistogramV00, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (SCAN_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
      dataScanDebugHost);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
      dataScanDebugDevice);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataScanErrorBuffer, dataScanErrorSizeBytes, dataScanError);
      if(dataScanError[0] != 0)
      {
        std::cout << "SCAN_ERROR_TRACK_ENABLE, kernelScanHistogramV00: "
          << "Received error code from device: " << dataScanError[0] << std::endl;
        verified = false; break;
      }
    }
#endif
#if SCAN_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
    if(verifyKernelScanHistogramV00(currentTimeSlot) != SDK_SUCCESS)
    {
      std::cout << "Failed verifyKernelScanHistogramV00" << std::endl; verified = false; 
      break;
    }
#endif
    }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00
    {
/*Unit test initialization*/
#if !(EXPAND_EVENTS_ENABLE && SCAN_ENABLE_V00)
    {
      double deltaDev = 50.0; double perecentInh = 5.0;
      cl_uint detla = cl_uint(((GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE*deltaDev)/100)/
        GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS);
      if(initializeDataForKernelGroupEventsV00(detla, deltaDev, perecentInh)!= SDK_SUCCESS)
        {std::cout << "Failed initializeDataForKernelGroupEventsV00" 
        << std::endl; verified = false; break;}
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventCountsBuffer, dataUnsortedEventCountsSizeBytes, 
      dataUnsortedEventCounts);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
#endif 
/*End unit test initialization*/
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if GROUP_EVENTS_VERIFY_ENABLE || SORT_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramBuffer, dataHistogramSizeBytes, dataHistogram);
#if SORT_VERIFY_ENABLE && \
    (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 && \
     GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE)
    if(captureUnsortedEvents(dataUnsortedEventCounts, dataUnsortedEventTargets, 
       dataUnsortedEventDelays, dataUnsortedEventWeights) != SDK_SUCCESS)
    {
      std::cout << "Failed captureUnsortedEvents" << std::endl; verified = false; break;
    }
#endif
#endif
    SET_KERNEL_ARG(kernelGroupEventsV00, currentTimeStep, argNumGrupEventsV00);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelGroupEventsV00, cl::NullRange, 
      globalThreadsGroupEventsV00, localThreadsGroupEventsV00, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelGroupEventsV00, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::cout << "GROUP_EVENTS_ERROR_TRACK_ENABLE, verifyKernelGroupEventsV00: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        verified = false; 
        break;
      }
    }
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    if(verifyKernelGroupEventsV00(GROUP_EVENTS_SOURCE_KEY_OFFSET_V00) != SDK_SUCCESS)
      {std::cout << "Failed verifyKernelGroupEventsV00" 
      << std::endl; verified = false; break;}
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
#if SCAN_ENABLE_V01
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    if(initializeDataForKernelScanHistogramV01() != SDK_SUCCESS){std::cout 
      << "Failed initializeDataForKernelScanHistogramV01" << std::endl; verified = false; break;}
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif 
/*End unit test initialization*/
#if (SCAN_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
      dataScanDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
      dataScanDebugDevice);
#endif
#if SCAN_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    SET_KERNEL_ARG(kernelScanHistogramV01, dataHistogramGroupEventsTikBuffer, argNumScan01);
    /*TODO: get rid of unnecessary argument with preprocessor*/
    SET_KERNEL_ARG(kernelScanHistogramV01, 0, argNumScan01+1);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelScanHistogramV01, cl::NullRange, 
      globalThreadsScanV01, localThreadsScanV01, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelScanHistogramV01, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (SCAN_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_FALSE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
      dataScanDebugHost);
    ENQUEUE_READ_BUFFER(CL_FALSE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
      dataScanDebugDevice);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_FALSE, dataScanErrorBuffer, dataScanErrorSizeBytes, dataScanError);
      if(dataScanError[0] != 0)
      {
        std::cout << "SCAN_ERROR_TRACK_ENABLE, kernelScanHistogramV01: "
          << "Received error code from device: " << dataScanError[0] << std::endl;
        verified = false; break;
      }
    }
#endif
#if SCAN_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsVerifySizeBytes, dataHistogramGroupEventsVerify);
    if(verifyKernelScanHistogramV01() != SDK_SUCCESS){std::cout 
      << "Failed kernelScanHistogramV01" 
      << std::endl; verified = false; break;}
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
    if(initializeDataForKernelGroupEventsV01(sortStep-1, 1) != SDK_SUCCESS){std::cout 
      << "Failed initializeDataForKernelGroupEventsV01" << std::endl; verified = false; break;}
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif  
/*End unit test initialization*/
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    SET_KERNEL_ARG(kernelGroupEventsV01, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV01++);
#endif
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTikBuffer, argNumGrupEventsV01);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTokBuffer, argNumGrupEventsV01+1);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataHistogramGroupEventsTikBuffer, argNumGrupEventsV01+2);
    SET_KERNEL_ARG(kernelGroupEventsV01, sortStep, argNumGrupEventsV01+3);
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    argNumGrupEventsV01--;
#endif
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelGroupEventsV01, cl::NullRange, 
      globalThreadsGroupEventsV01, localThreadsGroupEventsV01, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelGroupEventsV01, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::cout << "GROUP_EVENTS_ERROR_TRACK_ENABLE, verifyKernelGroupEventsV01: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        verified = false; 
        break;
      }
    }
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV01(sortStep) != SDK_SUCCESS){std::cout 
      << "Failed verifyKernelGroupEventsV01" << std::endl; verified = false; break;}
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
    if(initializeDataForKernelGroupEventsV02_V03(sortStep-1, 1) != SDK_SUCCESS)
    {
      std::cout 
      << "Failed initializeDataForKernelGroupEventsV02_V03" << std::endl; 
      verified = false; 
      break;
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif  
/*End unit test initialization*/
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    SET_KERNEL_ARG(kernelGroupEventsV02, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV02++);
#endif
    SET_KERNEL_ARG(kernelGroupEventsV02, dataGroupEventsTikBuffer, argNumGrupEventsV02);
    SET_KERNEL_ARG(kernelGroupEventsV02, dataGroupEventsTokBuffer, argNumGrupEventsV02+1);
    SET_KERNEL_ARG(kernelGroupEventsV02, dataHistogramGroupEventsTikBuffer, argNumGrupEventsV02+2);
    SET_KERNEL_ARG(kernelGroupEventsV02, sortStep, argNumGrupEventsV02+3);
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    argNumGrupEventsV02--;
#endif
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelGroupEventsV02, cl::NullRange, 
      globalThreadsGroupEventsV02, localThreadsGroupEventsV02, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelGroupEventsV02, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::cout << "GROUP_EVENTS_ERROR_TRACK_ENABLE, kernelGroupEventsV02: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        verified = false; 
        break;
      }
    }
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV02(sortStep, GROUP_EVENTS_REPLACEMENT_KEY_OFFSET_V02) != SDK_SUCCESS)
      {std::cout << "Failed verifyKernelGroupEventsV02" << std::endl; verified = false; break;}
#endif
    }
#endif
/**************************************************************************************************/
    }/*if*/
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02
    /*Swap buffers for next iteration*/
    swap2(dataHistogramGroupEventsTikBuffer, dataHistogramGroupEventsTokBuffer);
    swap2(dataGroupEventsTikBuffer, dataGroupEventsTokBuffer);
#endif
    }/*for*/
    if(!verified)  {break;}
/**************************************************************************************************/
#if (GROUP_EVENTS_TOTAL_NEURON_BITS && EXPAND_EVENTS_TOTAL_NEURON_BITS) &&\
    (GROUP_EVENTS_TOTAL_NEURON_BITS != EXPAND_EVENTS_TOTAL_NEURON_BITS)
  #error (GROUP_EVENTS_TOTAL_NEURON_BITS != EXPAND_EVENTS_TOTAL_NEURON_BITS)
#endif

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
#if SCAN_ENABLE_V01
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && \
      GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT && GROUP_EVENTS_ENABLE_V02)
    if(initializeDataForKernelScanHistogramV01() != SDK_SUCCESS){std::cout 
      << "Failed initializeDataForKernelScanHistogramV01" << std::endl; verified = false; break;}
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif 
/*End unit test initialization*/
#if (SCAN_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
      dataScanDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_FALSE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
      dataScanDebugDevice);
#endif
#if SCAN_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    SET_KERNEL_ARG(kernelScanHistogramV01, dataHistogramGroupEventsTikBuffer, argNumScan01);
    /*TODO: get rid of unnecessary argument with preprocessor*/
    SET_KERNEL_ARG(kernelScanHistogramV01, 0, argNumScan01+1);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelScanHistogramV01, cl::NullRange, 
      globalThreadsScanV01, localThreadsScanV01, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelScanHistogramV01, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (SCAN_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_FALSE, dataScanDebugHostBuffer, dataScanDebugHostSizeBytes, 
      dataScanDebugHost);
    ENQUEUE_READ_BUFFER(CL_FALSE, dataScanDebugDeviceBuffer, dataScanDebugDeviceSizeBytes, 
      dataScanDebugDevice);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_FALSE, dataScanErrorBuffer, dataScanErrorSizeBytes, dataScanError);
      if(dataScanError[0] != 0)
      {
        std::cout << "SCAN_ERROR_TRACK_ENABLE, kernelScanHistogramV01: "
          << "Received error code from device: " << dataScanError[0] << std::endl;
        verified = false; break;
      }
    }
#endif
#if SCAN_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsVerifySizeBytes, dataHistogramGroupEventsVerify);
    if(verifyKernelScanHistogramV01() != SDK_SUCCESS){std::cout 
      << "Failed kernelScanHistogramV01" 
      << std::endl; verified = false; break;}
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
      GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT && GROUP_EVENTS_ENABLE_V02)
    if(initializeDataForKernelGroupEventsV01(sortStep-1, 0) != SDK_SUCCESS){std::cout 
      << "Failed initializeDataForKernelGroupEventsV01" << std::endl; verified = false; break;}
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif  
/*End unit test initialization*/
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    SET_KERNEL_ARG(kernelGroupEventsV01, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV01++);
#endif
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTikBuffer, argNumGrupEventsV01);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataGroupEventsTokBuffer, argNumGrupEventsV01+1);
    SET_KERNEL_ARG(kernelGroupEventsV01, dataHistogramGroupEventsTikBuffer, argNumGrupEventsV01+2);
    SET_KERNEL_ARG(kernelGroupEventsV01, sortStep, argNumGrupEventsV01+3);
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    argNumGrupEventsV01--;
#endif
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelGroupEventsV01, cl::NullRange, 
      globalThreadsGroupEventsV01, localThreadsGroupEventsV01, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelGroupEventsV01, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::cout << "GROUP_EVENTS_ERROR_TRACK_ENABLE, kernelGroupEventsV01: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        verified = false; 
        break;
      }
    }
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV01(sortStep) != SDK_SUCCESS){std::cout 
      << "Failed verifyKernelGroupEventsV01" << std::endl; verified = false; break;}
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
      GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT && GROUP_EVENTS_ENABLE_V02)
    if(initializeDataForKernelGroupEventsV02_V03(sortStep-1, 0) != SDK_SUCCESS)
    {
      std::cout 
      << "Failed initializeDataForKernelGroupEventsV02_V03" << std::endl; 
      verified = false; 
      break;
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
#if GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#endif
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif  
/*End unit test initialization*/
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventTargetsBuffer, dataUnsortedEventTargetsSizeBytes, 
      dataUnsortedEventTargets);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventDelaysBuffer, dataUnsortedEventDelaysSizeBytes, 
      dataUnsortedEventDelays);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
      dataUnsortedEventWeights);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTikBuffer, 
      dataHistogramGroupEventsTikSizeBytes, dataHistogramGroupEventsTik);
#elif SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataUnsortedEventWeightsBuffer, dataUnsortedEventWeightsSizeBytes, 
        dataUnsortedEventWeights);
    }
#endif
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    SET_KERNEL_ARG(kernelGroupEventsV03, dataHistogramGroupEventsTokBuffer, argNumGrupEventsV03++);
#endif
    SET_KERNEL_ARG(kernelGroupEventsV03, dataGroupEventsTikBuffer, argNumGrupEventsV03);
    SET_KERNEL_ARG(kernelGroupEventsV03, dataGroupEventsTokBuffer, argNumGrupEventsV03+1);
    SET_KERNEL_ARG(kernelGroupEventsV03, dataHistogramGroupEventsTikBuffer, argNumGrupEventsV03+2);
    SET_KERNEL_ARG(kernelGroupEventsV03, sortStep, argNumGrupEventsV03+3);
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    argNumGrupEventsV03--;
#endif
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelGroupEventsV03, cl::NullRange, 
      globalThreadsGroupEventsV03, localThreadsGroupEventsV03, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelGroupEventsV03, (endAppTime-startAppTime), 1.0)
    }
#endif
#if (GROUP_EVENTS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugHostGroupEventsBuffer, dataDebugHostGroupEventsSizeBytes, 
      dataDebugHostGroupEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataDebugDeviceGroupEventsBuffer, 
      dataDebugDeviceGroupEventsSizeBytes, dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataErrorGroupEventsBuffer, dataErrorGroupEventsSizeBytes, 
        dataErrorGroupEvents);
      if(dataErrorGroupEvents[0] != 0)
      {
        std::cout << "GROUP_EVENTS_ERROR_TRACK_ENABLE, verifyKernelGroupEventsV03: "
          << "Received error code from device: " << dataErrorGroupEvents[0] << std::endl;
        verified = false; 
        break;
      }
    }
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
#endif
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTokBuffer, 
      dataGroupEventsTokSizeBytes, dataGroupEventsTok);
    if(verifyKernelGroupEventsV03(sortStep) != SDK_SUCCESS)
    {
      std::cout << "Failed verifyKernelGroupEventsV03, sort step " << sortStep << std::endl; 
      verified = false; 
      break;
    }
#endif
    }
#endif
/**************************************************************************************************/
    }/*if*/
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V03
    /*Swap buffers for next iteration*/
    swap2(dataHistogramGroupEventsTikBuffer, dataHistogramGroupEventsTokBuffer);
    swap2(dataGroupEventsTikBuffer, dataGroupEventsTokBuffer);
#endif
    }/*for*/
    if(!verified)  {break;}
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02 &&\
      GROUP_EVENTS_ENABLE_V03 && SCAN_ENABLE_V01)
    /*Use milder mode if UPDATE_NEURONS_ENABLE_V00*/
    cl_uint mode = (currentTimeStep!=0)*(1 + !UPDATE_NEURONS_ENABLE_V00);
    if(initializeDataForKernelMakeEventPtrs(mode, currentTimeStep-1) != SDK_SUCCESS)
    {
      std::cout 
      << "Failed initializeDataForKernelMakeEventPtrs" << std::endl; 
      verified = false; 
      break;
    }
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
#endif
/*End unit test initialization*/
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataMakeEventPtrsDebugHostBuffer, 
      dataMakeEventPtrsDebugHostSizeBytes, dataMakeEventPtrsDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataMakeEventPtrsDebugDeviceBuffer, 
      dataMakeEventPtrsDebugDeviceSizeBytes, dataMakeEventPtrsDebugDevice);
#endif
#if MAKE_EVENT_PTRS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataHistogramGroupEventsTokBuffer, 
      dataHistogramGroupEventsTokSizeBytes, dataHistogramGroupEventsTok);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
#endif
    SET_KERNEL_ARG(kernelMakeEventPtrs, dataHistogramGroupEventsTokBuffer, argNumMakeEventPtrs);
    SET_KERNEL_ARG(kernelMakeEventPtrs, dataGroupEventsTikBuffer, argNumMakeEventPtrs+1);
    SET_KERNEL_ARG(kernelMakeEventPtrs, dataMakeEventPtrsStructBuffer, argNumMakeEventPtrs+2); /*TODO: instead of a new buffer it could be a reuse*/
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelMakeEventPtrs, cl::NullRange, 
      globalThreadsMakeEventPtrs, localThreadsMakeEventPtrs, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelMakeEventPtrs, (endAppTime-startAppTime), 1.0)
    }
#endif
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelGlueEventPtrs, cl::NullRange, 
      globalThreadsGlueEventPtrs, localThreadsGlueEventPtrs, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelGlueEventPtrs, (endAppTime-startAppTime), 1.0)
    }
#endif
#endif
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsDebugHostBuffer, 
      dataMakeEventPtrsDebugHostSizeBytes, dataMakeEventPtrsDebugHost);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsDebugDeviceBuffer, 
      dataMakeEventPtrsDebugDeviceSizeBytes, dataMakeEventPtrsDebugDevice);
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsErrorBuffer, dataMakeEventPtrsErrorSizeBytes, 
        dataMakeEventPtrsError);
      if(dataMakeEventPtrsError[0] != 0)
      {
        std::cout << "MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE: Received error code from device: " 
          << dataMakeEventPtrsError[0] << std::endl;
        verified = false; 
        break;
      }
    }
#endif
#if MAKE_EVENT_PTRS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsStructBuffer, 
      dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
    if(verifyKernelMakeEventPtrs() != SDK_SUCCESS)
    {
      std::cout << "Failed verifyKernelMakeEventPtrs" << std::endl; 
      verified = false; 
      break;
    }
#endif
    }
#endif
/**************************************************************************************************/
#if (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 && \
     GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE)
#if SORT_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
      dataGroupEventsTikSizeBytes, dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsStructBuffer, 
      dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);
    if(verifySortedEvents(dataGroupEventsTik, dataMakeEventPtrsStruct, 2) != SDK_SUCCESS){std::cout 
      << "Failed verifySortedEvents" << std::endl; verified = false; break;}
#elif SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, 
        dataGroupEventsTikSizeBytes, dataGroupEventsTik);
      ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsStructBuffer, 
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
    if(initializeDataForKernelUpdateNeurons(!(MAKE_EVENT_PTRS_ENABLE), 
      resetVarsAndParams, resetVarsAndParams, currentTimeStep, "") != SDK_SUCCESS)
    {
      std::cout 
      << "Failed initializeDataForKernelUpdateNeurons" << std::endl; 
      verified = false; 
      break;
    }
    if(resetVarsAndParams)
    {
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
      ENQUEUE_WRITE_BUFFER(CL_TRUE, psToleranceBuffer, psToleranceSizeBytes, 
        psTolerance);
#endif
      ENQUEUE_WRITE_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
        dataSpikePackets);
      ENQUEUE_WRITE_BUFFER(CL_TRUE, constantCoefficientsBuffer, constantCoefficientsSizeBytes, 
        constantCoefficients);
      ENQUEUE_WRITE_BUFFER(CL_TRUE, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_WRITE_BUFFER(CL_TRUE, modelParametersBuffer, modelParametersSizeBytes, 
        modelParameters);
    }
#if !(MAKE_EVENT_PTRS_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
      dataMakeEventPtrsStruct);
# else
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
      dataMakeEventPtrsStruct);
#endif
#endif
/*END unit test initialization*/

/*Debugging*/
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUpdateNeuronsDebugHostBuffer, 
      dataUpdateNeuronsDebugHostSizeBytes, dataUpdateNeuronsDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUpdateNeuronsDebugDeviceBuffer, 
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
    (PROFILING_MODE == 1 && START_PROFILING_AT_STEP > 0) || SIMULATION_SNAPSHOT
#if UPDATE_NEURONS_VERIFY_ENABLE || (KERNEL_ENDSTEP_VERIFY_EVERY_STEPS > 0)
    if(!((currentTimeStep+1)%KERNEL_ENDSTEP_VERIFY_EVERY_STEPS) || (currentTimeStep == 0)
#else
    if(false
#endif
#if OVERWRITE_SPIKES_UNTILL_STEP
    || (currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP)
#endif
#if INJECT_CURRENT_UNTILL_STEP
    || (currentTimeStep < INJECT_CURRENT_UNTILL_STEP)
#endif
#if SIMULATION_SNAPSHOT
    || (takeSimSnapshot)
#endif
    ){
      ENQUEUE_READ_BUFFER(CL_TRUE, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
        dataMakeEventPtrsStruct);
    }
#endif
#endif
/*END Get dataMakeEventPtrsStruct*/

/*Set loop-dependent arguments for update kernel and equeue it*/
    SET_KERNEL_ARG(kernelUpdateNeuronsV00, dataGroupEventsTikBuffer, argNumUpdateNeuronsV00);
    SET_KERNEL_ARG(kernelUpdateNeuronsV00, currentTimeStep, argNumUpdateNeuronsV00+1);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelUpdateNeuronsV00, cl::NullRange, 
      globalThreadsUpdateNeuronsV00, localThreadsUpdateNeuronsV00, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelUpdateNeuronsV00, (endAppTime-startAppTime), 1.0)
    }
#endif
/*END enqueue update kernel*/

/*Set loop-dependent arguments for update kernel for spiked neurons and equeue it*/
    SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, dataGroupEventsTikBuffer, 
      argNumUpdateSpikedNeuronsV00);
    SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, currentTimeStep, 
      argNumUpdateSpikedNeuronsV00+1);
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP){startAppTime = timeStampNs();}
#endif
    status = commandQueue.enqueueNDRangeKernel(kernelUpdateSpikedNeuronsV00, cl::NullRange, 
      globalThreadsUpdateNeuronsV00, localThreadsUpdateNeuronsV00, NULL, &ndrEvt);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "CommandQueue::enqueueNDRangeKernel() failed."))
    {
      return SDK_FAILURE;
    }
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
    if(currentTimeStep >= START_PROFILING_AT_STEP)
    {
      status = ndrEvt.wait();
      endAppTime = timeStampNs();
      if(!sampleCommon->checkVal(status, CL_SUCCESS, "commandQueue, cl:Event.wait() failed."))
      {
        return SDK_FAILURE;
      }
      REGISTER_TIME(kernelUpdateSpikedNeuronsV00, (endAppTime-startAppTime), 1.0)
    }
#endif
/*END enqueue update kernel for spiked neurons*/

/*Debugging: buffer exchange*/
#if (UPDATE_NEURONS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUpdateNeuronsDebugHostBuffer, 
      dataUpdateNeuronsDebugHostSizeBytes, dataUpdateNeuronsDebugHost);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataUpdateNeuronsDebugDeviceBuffer, 
      dataUpdateNeuronsDebugDeviceSizeBytes, dataUpdateNeuronsDebugDevice);
#endif
/*END debugging*/

/*Error tracking: read error mask*/
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
    if(!((currentTimeStep+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))
    {
      bool ignoreSolverFailuresDevice = IGNORE_SOLVER_EXCEPTIONS;
      ENQUEUE_READ_BUFFER(CL_TRUE, dataUpdateNeuronsErrorBuffer, dataUpdateNeuronsErrorSizeBytes, 
        dataUpdateNeuronsError);
      if(dataUpdateNeuronsError[0] != 0)
      {
        std::cout << "UPDATE_NEURONS_ERROR_TRACK_ENABLE: received error code from the device: ";
        PRINT_HEX(4, dataUpdateNeuronsError[0]); std::cout << std::endl;
        
        if((!ignoreSolverFailuresDevice && dataUpdateNeuronsError[0] != 0) || 
           (ignoreSolverFailuresDevice && 
           ((UPDATE_NEURONS_ERROR_NON_SOLVER_FAILURE_MASK&dataUpdateNeuronsError[0]) != 0))
        ){
          verified = false; break;
        }
        else
        {
          std::cout << "UPDATE_NEURONS_ERROR_TRACK_ENABLE: ignoring solver failures\n";
          memset(dataUpdateNeuronsError, 0, UPDATE_NEURONS_ERROR_BUFFER_SIZE_WORDS*sizeof(cl_uint));
          ENQUEUE_WRITE_BUFFER(CL_TRUE, dataUpdateNeuronsErrorBuffer, 
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
#if OVERWRITE_SPIKES_UNTILL_STEP
    verify |= (currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP);
#endif
#if INJECT_CURRENT_UNTILL_STEP
    verify |= (currentTimeStep < INJECT_CURRENT_UNTILL_STEP);
#endif
    if(verify)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);
    }
#if OVERWRITE_SPIKES_UNTILL_STEP
    if(currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP)
    {
      cl_uint size = (UPDATE_NEURONS_SPIKE_PACKETS_V00*UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS);
      cl_uint *dataSpikePacketsToInject = (cl_uint *)calloc(size, sizeof(cl_uint));
      memcpy(dataSpikePacketsToInject, dataSpikePackets, (size*sizeof(cl_uint)));
  
      ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
        dataSpikePackets);
        
      if(verifyKernelUpdateNeurons
        (
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          dataSynapsePointer,
          dataSynapseTargets,
          dataSynapseDelays,
          dataSynapseWeights,
          dataSpikePacketsToInject,
          modelVariables,
          dataSpikePackets
        ) != SDK_SUCCESS
      ){
        std::cout << "Failed verifyKernelUpdateNeurons" << std::endl; 
        verified = false; 
      }
      free(dataSpikePacketsToInject);
      if(verified == false){break;}
    }
    else
#endif
    {
      if(verify)
      {
        ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
          dataSpikePackets);
      }
      if(verifyKernelUpdateNeurons
          (
            verify,
            currentTimeStep,
            dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
            dataMakeEventPtrsStruct,
            dataGroupEventsTik,
            dataSynapsePointer,
            dataSynapseTargets,
            dataSynapseDelays,
            dataSynapseWeights,
            NULL,
            modelVariables,
            dataSpikePackets
          ) != SDK_SUCCESS
      ){
        std::cout << "Failed verifyKernelUpdateNeurons" << std::endl; 
        verified = false; 
        break;
      }
    }
    }
/*END: Verification without profiling*/

/*Profiling mode 1: verify host-device until steps reach OVERWRITE_SPIKES_UNTILL_STEP (while
  spikes or currents are injected)*/
#elif (PROFILING_MODE == 1 && START_PROFILING_AT_STEP > 0)
#if OVERWRITE_SPIKES_UNTILL_STEP
    if(currentTimeStep < OVERWRITE_SPIKES_UNTILL_STEP)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);
      
      cl_uint size = (UPDATE_NEURONS_SPIKE_PACKETS_V00*UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS);
      cl_uint *dataSpikePacketsToInject = (cl_uint *)calloc(size, sizeof(cl_uint));
      memcpy(dataSpikePacketsToInject, dataSpikePackets, (size*sizeof(cl_uint)));
  
      ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
        dataSpikePackets);
        
      if(verifyKernelUpdateNeurons
        (
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          dataSynapsePointer,
          dataSynapseTargets,
          dataSynapseDelays,
          dataSynapseWeights,
          dataSpikePacketsToInject,
          modelVariables,
          dataSpikePackets
        ) != SDK_SUCCESS)
      {
        std::cout << "Failed verifyKernelUpdateNeurons" << std::endl; 
        verified = false; 
      }
      free(dataSpikePacketsToInject);
      if(verified == false){break;}
    }
#elif INJECT_CURRENT_UNTILL_STEP
    if(currentTimeStep < INJECT_CURRENT_UNTILL_STEP)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, modelVariablesBuffer, modelVariablesSizeBytes, 
        modelVariables);
      ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
        dataGroupEventsTik);
      ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
        dataSpikePackets);
      if(verifyKernelUpdateNeurons
        (
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          dataSynapsePointer,
          dataSynapseTargets,
          dataSynapseDelays,
          dataSynapseWeights,
          NULL,
          modelVariables,
          dataSpikePackets
        ) != SDK_SUCCESS)
      {
        std::cout << "Failed verifyKernelUpdateNeurons" << std::endl; 
        verified = false; 
        break;
      }
    }
#endif
#endif
/*END verification during inject of spikes or currents*/

/*Unit test verification*/
#else
    ENQUEUE_READ_BUFFER(CL_TRUE, modelVariablesBuffer, modelVariablesSizeBytes, 
      modelVariables);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
      dataSpikePackets);
    if(verifyKernelUpdateNeurons
        (
          true,
          currentTimeStep,
          dataGroupEventsTikSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          dataGroupEventsTik,
          NULL,
          NULL,
          NULL,
          NULL,
          NULL,
          modelVariables,
          dataSpikePackets
        ) != SDK_SUCCESS
    ){
      std::cout << "Failed verifyKernelUpdateNeurons" << std::endl; 
      verified = false; 
      break;
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
    ENQUEUE_READ_BUFFER(CL_TRUE, dataSpikePacketsBuffer, dataSpikePacketsSizeBytes, 
      dataSpikePackets);
    ENQUEUE_READ_BUFFER(CL_TRUE, modelVariablesBuffer, modelVariablesSizeBytes, 
      modelVariables);
    ENQUEUE_READ_BUFFER(CL_TRUE, dataGroupEventsTikBuffer, dataGroupEventsTikSizeBytes, 
      dataGroupEventsTik);
        
    takeSimulationSnapshot
    (
      currentTimeStep,
      1000,
      dataUnsortedEventCounts,
      dataUnsortedEventWeights,
      dataMakeEventPtrsStruct,
      dataGroupEventsTik,
      dataSpikePackets,
      modelVariables
    );
  }
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V03
#if ((GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS%2) == 0)
    /*Restore original buffer pointer for the next step*/
    swap2(dataHistogramGroupEventsTikBuffer, dataHistogramGroupEventsTokBuffer);
    swap2(dataGroupEventsTikBuffer, dataGroupEventsTokBuffer);
#endif
#endif
/**************************************************************************************************/
  }/*for SIMULATION_TIME_STEPS*/
  
  std::cout << "\nDispatched all kernels" << std::endl;
  status = commandQueue.flush();
  if(!sampleCommon->checkVal(status, CL_SUCCESS, "cl::CommandQueue.flush failed."))
  {
    return SDK_FAILURE;
  }
  std::cout << "\nEnqueued all kernels" << std::endl;
  std::cout << "\nWaiting for the command queue to complete..." << std::endl;
  eventStatus = CL_QUEUED;
  while(eventStatus != CL_COMPLETE)
  {
    status = ndrEvt.getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS, &eventStatus);
    if(!sampleCommon->checkVal(status, CL_SUCCESS, 
      "commandQueue, cl:Event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS) failed."))
        return SDK_FAILURE;
  }
  
#if PROFILING_MODE == 1 && START_PROFILING_AT_STEP > 0
  status = commandQueue.finish();
  endAppTime = timeStampNs();
  if(!sampleCommon->checkVal(status, CL_SUCCESS, "cl::CommandQueue.finish failed."))
  {
    return SDK_FAILURE;
  }
  cl_uint profiledSteps = SIMULATION_TIME_STEPS - START_PROFILING_AT_STEP;
  double appTimeSec = (endAppTime - startAppTime);
  std::cout << "\nDevice timing: " << profiledSteps << " steps, " << appTimeSec << " sec, " 
    << appTimeSec/profiledSteps << " sec/step" << std::endl;
  
  /*Host equivalent computation*/
#if (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 &&\
    GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE &&\
    EXPAND_EVENTS_ENABLE)
  {
    int result = SDK_SUCCESS;
    double startAppTime = 0, endAppTime = 0;
    
    std::cout << "\nStarting host timing at step " << START_PROFILING_AT_STEP << std::endl;
    
    startAppTime = timeStampNs();
    for(cl_uint step = START_PROFILING_AT_STEP; step < SIMULATION_TIME_STEPS; step++)
    {
      result = propagateSpikes
      (
        UPDATE_NEURONS_TOTAL_NEURONS,
        REFERENCE_EVENT_QUEUE_SIZE,
        nrn_ps,
        ne,
        te_ps,
        dataSynapsePointer,
        dataSynapseTargets,
        dataSynapseDelays,
        dataSynapseWeights
      );
      if(result != SDK_SUCCESS){break;}

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
      if(result != SDK_SUCCESS && !IGNORE_SOLVER_EXCEPTIONS){break;}
    }
    endAppTime = timeStampNs();
    
    cl_uint profiledSteps = SIMULATION_TIME_STEPS - START_PROFILING_AT_STEP;
    double appTimeSec = (endAppTime - startAppTime);
    std::cout << "\nHost timing: " << profiledSteps << " steps, " << appTimeSec << " sec, " 
      << appTimeSec/profiledSteps << " sec/step" << std::endl;
  }
#endif
#elif PROFILING_MODE == 2 && START_PROFILING_AT_STEP > 0
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
    LOG("Kernel " << *kernel_name << " Total Time:" << execTime["Time"], 1);
    LOG("Kernel " << *kernel_name << " Average Time:" << averageTime, 1);
    LOG("Kernel " << *kernel_name << " Execution Count:" << (cl_uint)execTime["Count"], 1);
  }
  std::cout << "Total time: " << totalExecTime << std::endl;
  LOG("Total Time:" << totalExecTime, 1);
  std::cout.setf(std::ios::scientific,std::ios::floatfield);
#endif

  if(!verified)  {return SDK_FAILURE;}

#if (ERROR_TRY_CATCH_ENABLE)
  }
  catch (cl::Error err) 
  {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    return SDK_FAILURE;
  }
#endif

  return SDK_SUCCESS;
}



int 
IntegrationTest::initialize()
{
  /*Call base class Initialize to get default configuration*/
  if(!this->SDKSample::initialize())
      return SDK_FAILURE;
  
  srandSeed = 0; 
  srandCounter = 0;
  
  time_t rawtime; time ( &rawtime ); startTimeStamp = std::string(ctime (&rawtime));
  
#ifdef WIN32
  QueryPerformanceFrequency((LARGE_INTEGER *)&performanceFrequency);
#endif
      
#if (LOG_SIMULATION)
  simulationLogFile = new std::ofstream(LOG_SIMULATION_FILE_NAME, 
    std::ofstream::out | std::ofstream::app);
  dataToSimulationLogFile = new std::stringstream("", std::ios::out | std::ios::app);

  if(!(*simulationLogFile).is_open())
  {
    std::cerr << "ERROR, initialize: Unable to open log file";
    return SDK_FAILURE;
  }
  
  *dataToSimulationLogFile << "-------\n" << startTimeStamp;
#endif

#if (LOG_REPORT)
    reportLogFile = new std::ofstream(LOG_REPORT_FILE_NAME);
    dataToReportLogFile = new std::stringstream("", std::ios::out | std::ios::app);

    if(!(*reportLogFile).is_open())
    {
      std::cerr << "ERROR, initialize: Unable to open report file";
      return SDK_FAILURE;
    }
#endif

#if (SIMULATION_SNAPSHOT)
    snapshotLogFile = new std::ofstream(LOG_SNAPSHOT_FILE_NAME);
    dataToSnapshotLogFile = new std::stringstream("", std::ios::out | std::ios::app);

    if(!(*snapshotLogFile).is_open())
    {
      std::cerr << "ERROR, initialize: Unable to open snapshot log file";
      return SDK_FAILURE;
    }
#endif

  return SDK_SUCCESS;
}



int 
IntegrationTest::setup()
{
  /* Allocate host memory*/
  if(allocateHostData() != SDK_SUCCESS)
      return SDK_FAILURE;
      
  /* Init host memory
  if(initializeHostData() != SDK_SUCCESS)
      return SDK_FAILURE;*/
  
  /* Register local memory for stats */
  if(registerLocalMemory() != SDK_SUCCESS)
      return SDK_FAILURE;

  /* create and initialize timers */
  int timer = sampleCommon->createTimer();
  sampleCommon->resetTimer(timer);
  sampleCommon->startTimer(timer);

  if(setupCL() != SDK_SUCCESS)
      return SDK_FAILURE;

  sampleCommon->stopTimer(timer);
  /* Compute setup time */
  setupTime = (double)(sampleCommon->readTimer(timer));

  return SDK_SUCCESS;
}



int 
IntegrationTest::run()
{
  /*Warm up*/
  if(warmupcount > 0)
  {
    std::cout << "Warming up kernel for " << warmupcount << " iterations" <<std::endl;
    
    for(int i = 0; i < warmupcount; i++)
    {
      /* Set kernel arguments and run kernel */
      if(runCLKernels() != SDK_SUCCESS)
          return SDK_FAILURE;
    }
  }

  std::cout << "Executing test for " << iterations << 
      " iterations" <<std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  /* create and initialize timers */
  int timer = sampleCommon->createTimer();
  sampleCommon->resetTimer(timer);
  sampleCommon->startTimer(timer);

  for(int i = 0; i < iterations; i++)
  {
    /* Set kernel arguments and run kernel */
    if(runCLKernels() != SDK_SUCCESS)
        return SDK_FAILURE;
  }

  sampleCommon->stopTimer(timer);    
  /* Compute kernel time */
  kernelTime = (double)(sampleCommon->readTimer(timer)) / iterations;
  
  LOG("Result:PASS", 1);
  return SDK_SUCCESS;
}



int 
IntegrationTest::verifyResults()
{
  return SDK_SUCCESS;
}



int 
IntegrationTest::verifyKernelExpandEvents
(
  cl_uint timeStep,
  bool reset
)
/**************************************************************************************************/
{
  bool result = SDK_SUCCESS;
  
#if EXPAND_EVENTS_ENABLE
  {
  /*Reset data (useful for unit test)*/
  if(reset)
  {
    memset(dataUnsortedEventCountsVerify, 0, dataUnsortedEventCountsVerifySizeBytes);
    memset(dataUnsortedEventsTargetsVerify, 0, dataUnsortedEventsTargetsVerifySizeBytes);
    memset(dataUnsortedEventsDelaysVerify, 0, dataUnsortedEventsDelaysVerifySizeBytes);
    memset(dataUnsortedEventsWeightsVerify, 0, dataUnsortedEventsWeightsVerifySizeBytes);
    memset(dataHistogramVerify, 0, dataHistogramVerifySizeBytes);
  }
  
  /*Init verification data*/
  cl_uint maxSynapticBuffer = 0;
  cl_uint timeSlot = timeStep%EXPAND_EVENTS_TIME_SLOTS;

  /*Reset histogram data in previous step*/
  for(cl_uint i = 0; i < (EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*EXPAND_EVENTS_GRID_SIZE_WG+1); i++)
  {
    cl_uint offset = 
    /*time slot*/
    (((EXPAND_EVENTS_TIME_SLOTS + (timeSlot-1))%EXPAND_EVENTS_TIME_SLOTS)*
    (EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*EXPAND_EVENTS_GRID_SIZE_WG+1) + 
    /*data*/
    i);
    dataHistogramVerify[offset] = 0;
  }
  /*Reset event counter in previous step*/
  for(cl_uint p = 0; p < EXPAND_EVENTS_SPIKE_PACKETS; p++)
  {
    cl_uint offset = 
      /*WG*/
      EXPAND_EVENTS_TIME_SLOTS * p/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG + 
      /*time slot*/
      ((EXPAND_EVENTS_TIME_SLOTS + (timeSlot-1))%EXPAND_EVENTS_TIME_SLOTS);
      
    dataUnsortedEventCountsVerify[offset] = 0;
  }
  /*Iterate through spike packets*/
  for(cl_uint packet = 0; packet < EXPAND_EVENTS_SPIKE_PACKETS; packet++)
  {
    cl_uint packet_index = packet * EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS;
    cl_uint total_spikes = dataSpikePackets[packet_index];

    /*Iterate through spikes in a current packet*/
    for(cl_uint i = 0; i < total_spikes; i++)
    {
      cl_uint spiked_neuron = dataSpikePackets[packet_index + 
        EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * i];
      cl_float spike_time = *((cl_float *)(&dataSpikePackets[packet_index + 
        EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS * i + 1]));
        
      cl_uint synapsePointer = dataSynapsePointer[spiked_neuron];
      cl_uint synapse_count = dataSynapsePointer[spiked_neuron + 1] - synapsePointer;
      
      /*Iterate through synapses of spiked neuron*/
      for(cl_uint j = 0; j < synapse_count; j++)
      {
        cl_uint synapseOffset = (synapsePointer + j);
        /*target neuron*/
        cl_uint target_neuron = dataSynapseTargets[synapseOffset];
        /*weight*/
        cl_float weight = dataSynapseWeights[synapseOffset];
        /*delay*/
        cl_float delay = dataSynapseDelays[synapseOffset];
        /*Add delay to spike time, decrement*/
        cl_float event_time = spike_time + delay; event_time = event_time - 1.0f;
        /*Make it relative to its time slot*/
        cl_float event_time_binned = event_time - (int)event_time;
        /*Events with 0.0 time bounce back to the previous time slot*/
        if(event_time_binned == 0.0f){event_time_binned = 1.0f;}
        cl_uint bin_correction = (int)event_time_binned;
        /*Calculate time slot*/
        cl_uint time_slot = 
          ((timeSlot + (int)event_time - bin_correction)%EXPAND_EVENTS_TIME_SLOTS);
        
        /*Obtain offsets for storing synaptic event*/
        cl_uint event_local_ptr = dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * 
          packet/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG + time_slot];
        
        /*Catch overflows and record max overflow*/
        if(event_local_ptr >= EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
        {
          if(maxSynapticBuffer < event_local_ptr)
          {
            maxSynapticBuffer = event_local_ptr;
          }
          event_local_ptr = EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE-1;
        }

        /*Increment event counter*/
        dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * 
          packet/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG + time_slot]++;
        
        cl_uint event_global_ptr = 
          /*Event data buffers*/
          (packet/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG) * EXPAND_EVENTS_TIME_SLOTS * 
          EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
          /*Current event data buffer*/
          time_slot * EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
          /*Current event*/
          event_local_ptr;

        /*Store the event*/
        dataUnsortedEventsTargetsVerify[event_global_ptr] = target_neuron; 
        *((cl_float *)(&dataUnsortedEventsDelaysVerify[event_global_ptr])) = event_time_binned; 
        *((cl_float *)(&dataUnsortedEventsWeightsVerify[event_global_ptr])) = weight;
        
#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
        /*Compute histogram key for target neuron based on MSBs*/
        cl_uint event_time_uint = *((cl_uint *)(&event_time_binned));
        cl_uint bin = (event_time_uint>>EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT) & 
          EXPAND_EVENTS_HISTOGRAM_BIN_MASK;
          
        /*Offset is based on time slot, bin, WG*/
        cl_uint offset = 
          /*WG offset*/
          (packet/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG) + 
          /*time slot + bin with EXPAND_EVENTS_GRID_SIZE_WG as a pitch*/
          (time_slot*(EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*EXPAND_EVENTS_GRID_SIZE_WG+1) + 
          bin*EXPAND_EVENTS_GRID_SIZE_WG);
          
        /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
        dataHistogramVerify[offset]++;
#endif
      }
    }
  }
  
  /*Fail if overflow*/
  if(maxSynapticBuffer != 0)
  {
    /*Correct time bin counters in case if simulation still proceeds*/
    for(cl_uint packet = 0; packet < EXPAND_EVENTS_SPIKE_PACKETS; packet++)
    {
      for(cl_uint time_slot = 0; time_slot < EXPAND_EVENTS_TIME_SLOTS; time_slot++)
      {
        if(dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * 
          packet/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG + time_slot] >= 
          EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
        {
          dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * 
            packet/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG + time_slot] = 
            EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE-1;
        }
      }
    }
    
    std::cout << "Event buffer overflow. " << 
      "Need to increase EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE above " 
      << maxSynapticBuffer << std::endl;
    return SDK_FAILURE;
  }
  
  cl_uint error_event_totals = 0;

  /*Iterate through synaptic counters*/
  for(cl_uint i = 0; i < EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS; i++)
  {
    if(dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * i + timeSlot] != 
      dataUnsortedEventCounts[EXPAND_EVENTS_TIME_SLOTS * i + timeSlot])
    {
      error_event_totals++;
      /**/
      std::cout << i <<"->(" << dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * i + 
        timeSlot] << "," << dataUnsortedEventCounts[EXPAND_EVENTS_TIME_SLOTS * i + 
        timeSlot] << "),";
      
    }
  }
  
  if(error_event_totals)
  {
    std::cout << "Failed to match synaptic event time slot counters " << error_event_totals 
      << " times." << std::endl;
    return SDK_FAILURE;
  }

  /*Iterate through spike packets*/
  for(cl_uint p = 0; p < EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS; p++)
  {
    cl_uint checksum_event_target_neuron_host = 0, checksum_event_target_neuron_device = 0;
    cl_uint checksum_event_delay_host = 0, checksum_event_delay_device = 0;
    cl_uint checksum_event_weight_host = 0, checksum_event_weight_device = 0;
    cl_uint event_total = dataUnsortedEventCountsVerify[EXPAND_EVENTS_TIME_SLOTS * p + timeSlot];
    
    /*Iterate through synapses of spiked neuron*/
    for(cl_uint j = 0; j < event_total; j++)
    {
      cl_uint event_global_ptr = 
        /*Event data buffers*/
        p * EXPAND_EVENTS_TIME_SLOTS * 
        EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event data buffer*/
        timeSlot * EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE +
        /*Current event*/
        j;
        
      /*Compute checksums*/
      //checksum_event_target_neuron_host ^= dataUnsortedEventsTargetsVerify[event_global_ptr];
      CHECKSUM01(checksum_event_target_neuron_host, dataUnsortedEventsTargetsVerify[event_global_ptr]);
      //checksum_event_target_neuron_device ^= dataUnsortedEventTargets[event_global_ptr];
      CHECKSUM01(checksum_event_target_neuron_device, dataUnsortedEventTargets[event_global_ptr]);

      //checksum_event_delay_host ^= dataUnsortedEventsDelaysVerify[event_global_ptr];
      CHECKSUM01(checksum_event_delay_host, dataUnsortedEventsDelaysVerify[event_global_ptr]);
      //checksum_event_delay_device ^= dataUnsortedEventDelays[event_global_ptr];
      CHECKSUM01(checksum_event_delay_device, dataUnsortedEventDelays[event_global_ptr]);
      cl_float timeA = *((cl_float *)(&dataUnsortedEventDelays[event_global_ptr]));
      cl_float timeV = *((cl_float *)(&dataUnsortedEventsDelaysVerify[event_global_ptr]));
      if(timeA < 0.0f || timeA > (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY))
      {
        std::cout << "Found an event in actual data with time outside of valid range: value " 
          << timeA << ", range " << 0.0f
          << " - " << (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY) << std::endl;
        return SDK_FAILURE;
      }
      if(timeV< 0.0f || timeV > (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY))
      {
        std::cout << "Found an event in verification data with delay outside of valid range: value " 
          << timeV << ", range " << 0.0f
          << " - " << (float)(EXPAND_EVENTS_MAX_DELAY-EXPAND_EVENTS_MIN_DELAY) << std::endl;
        return SDK_FAILURE;
      }
      
      //checksum_event_weight_host ^= dataUnsortedEventsWeightsVerify[event_global_ptr];
      CHECKSUM01(checksum_event_weight_host, dataUnsortedEventsWeightsVerify[event_global_ptr]);
      //checksum_event_weight_device ^= dataUnsortedEventWeights[event_global_ptr];
      CHECKSUM01(checksum_event_weight_device, dataUnsortedEventWeights[event_global_ptr]);
    }
    
    if(
      (checksum_event_target_neuron_host != checksum_event_target_neuron_device) ||
      (checksum_event_delay_host != checksum_event_delay_device) ||
      (checksum_event_weight_host != checksum_event_weight_device)
    ){
      std::cout << "Failed to verify sunaptic data checksum(s), time slot " 
      << timeSlot << ": " << std::endl;
      std::cout << "\tTarget neuron: host = " << checksum_event_target_neuron_host 
      << " vs device = " << 
        checksum_event_target_neuron_device << std::endl;
      std::cout << "\tDelay: host = " << checksum_event_delay_host << " vs device = " << 
        checksum_event_delay_device << std::endl;
      std::cout << "\tWeight host = " << checksum_event_weight_host << " vs device = " << 
        checksum_event_weight_device << std::endl;
        
      return SDK_FAILURE;
    }
  }

#if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
  /* init scan and verify data */
  for(cl_uint j = 0; j < EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS; j++)
  {
    for(cl_uint k = 0; k < EXPAND_EVENTS_GRID_SIZE_WG; k++)
    {
      /*Offset is based on time slot, bin, WG*/
      cl_uint offset = 
      /*Time slot space = (time slot #) x (bins per WG) x WGs*/
      timeSlot*(EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*EXPAND_EVENTS_GRID_SIZE_WG+1) + 
      /*Bin totals are stored aligned in GM (bin 0 from all WGs, bin 1 from all WGs etc*/
      j*EXPAND_EVENTS_GRID_SIZE_WG + k;
      /*Increment counter for a target neuron failing into slot/bin/WG specified by the offset.*/
      error_event_totals += (dataHistogram[offset] != dataHistogramVerify[offset]);
    }
  }
  
  if(error_event_totals)
  {
    std::cout << "Failed to match histogram " << error_event_totals 
      << " times." << std::endl;
    return SDK_FAILURE;
  }
#endif
  }
#endif

  return result;
}



int 
IntegrationTest::verifyKernelScanHistogramV00(cl_uint timeSlot)
/**************************************************************************************************/
{
#if SCAN_ENABLE_V00
  cl_uint errorEventTotals = 0; 
  cl_uint offsetTimeSlot = timeSlot*(SCAN_HISTOGRAM_TOTAL_BINS_V00*SCAN_HISTOGRAM_BIN_SIZE_V00 + 1);

  for(cl_uint j = 1; j <= SCAN_HISTOGRAM_TOTAL_BINS_V00*SCAN_HISTOGRAM_BIN_SIZE_V00; j++)
  {
    cl_uint offset = offsetTimeSlot + j;
    if(dataHistogramVerify[offset] != dataHistogram[offset])
    {
      errorEventTotals++;
      /*
      std::cout << timeSlot << "->(" << dataHistogramVerify[offset] << "," 
        << dataHistogram[offset] << "), ";
      */
    }
  }

  if(errorEventTotals)
  {
    std::cout << "verifyKernelScanHistogramV00: failed to match verification data " 
      << errorEventTotals << " times." << std::endl;
    return SDK_FAILURE;
  }
#endif

  return SDK_SUCCESS;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelGroupEventsV00(cl_uint keyOffset)
/**************************************************************************************************/
{
  bool result = SDK_SUCCESS;
#if GROUP_EVENTS_ENABLE_V00
  {
  cl_uint print = 0;
  cl_uint offset = currentTimeSlot*
    (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE + 1);
  
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_TIME_SLOTS*(GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE+1);
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogram, dataHistogramSizeBytes);
  
  /*Init data for verification*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  for(cl_uint i = 0; i < dataHistogramGroupEventsVerifySize; i++)
  {
    dataHistogramGroupEventsVerify[i] = 0;
  }
  
  cl_uint total_synaptic_events = 
    dataHistogram[offset+GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
  
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
#endif
  for(cl_uint i = 0; i < dataGroupEventsTikVerifySize; i++)
  {
    dataGroupEventsTikVerify[i] = 0;
  }
  
  /*Group events for verification*/
  for(cl_uint b = 0; b < GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS; b++)
  {
    cl_uint event_total = dataUnsortedEventCounts[GROUP_EVENTS_TIME_SLOTS*b + currentTimeSlot];
    
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
      cl_uint key = dataUnsortedEventDelays[ptr];
      
      /*Compute offset key for target neuron*/
      cl_uint bin = (key>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00)&
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
      if(dataOffsetGroupEventsCopy[bin_offset] >= dataHistogram[bin_offset+1])
      {
        std::cout << "verifyKernelGroupEventsV00: Destination event bin pointer overlaps "
          << "with next pointer for bin " << 
          bin << ", time slot " << currentTimeSlot << ", buffer " << b << std::endl;
        result = SDK_FAILURE; 
        break;
      }
      
      /*Calculate offset in the grouped data space*/
      cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
        
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
        result = SDK_FAILURE;
        break;
      }
      /*Increment counter for this bin.*/
      dataHistogramGroupEventsVerify[hist_out_ptr]++;
#endif
      
      /*Store event at its group location (grouped by bins)*/
      dataGroupEventsTikVerify[dest_offset] = key;
#if GROUP_EVENTS_VALUES_MODE_V00
      dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE+dest_offset] = ptr;
#endif
      /*Increment ptr for next data item*/
      dataOffsetGroupEventsCopy[bin_offset]++;
    }
    if(result != SDK_SUCCESS){break;}
  }
  
  free(dataOffsetGroupEventsCopy);
  
  print ? std::cout << "Time slot " << currentTimeSlot << ": " << std::endl, true : false;

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
    cl_uint start = dataHistogram[offset + j*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    cl_uint end = dataHistogram[offset + (j+1)*GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    
    print ? std::cout << "Start-End: " << start << "-" << end << std::endl, true : false;
    
    /*Verify correct bin and checksum in that bin for current time slot*/
    for(cl_uint p = start; p < end; p++)
    {
      cl_uint v_key = (dataGroupEventsTikVerify[p]>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00) & 
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
      cl_uint a_key = (dataGroupEventsTik[p]>>GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00) & 
        GROUP_EVENTS_HISTOGRAM_BIN_MASK;
      
      verify_error_count_target_neuron += (v_key != j);
      actual_error_count_target_neuron += (a_key != j);
      CHECKSUM01(verify_checksum_target_neuron, dataGroupEventsTikVerify[p]);
      //verify_checksum_target_neuron += dataGroupEventsTikVerify[p];
      CHECKSUM01(actual_checksum_target_neuron, dataGroupEventsTik[p]);
      //actual_checksum_target_neuron += dataGroupEventsTik[p];
#if GROUP_EVENTS_VALUES_MODE_V00
      CHECKSUM01(verify_checksum_value_01, dataGroupEventsTikVerify[p + 
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]);
      //verify_checksum_value_01 += dataGroupEventsTikVerify[p + 
      //  GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE];
      CHECKSUM01(actual_checksum_value_01, dataGroupEventsTik[p + 
        GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE]);
      //actual_checksum_value_01 += dataGroupEventsTik[p + 
      //  GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE];
#endif
    }
    
    if(verify_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to validate correct bin for " 
        << verify_error_count_target_neuron 
        << " keys out of " << (end-start) << " in verification data in bin " << j 
        << ", time slot " << currentTimeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
    
    if(actual_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to validate correct bin for " 
        << actual_error_count_target_neuron 
        << " keys out of " << (end-start) << " in actual data in bin " << j 
        << ", time slot " << currentTimeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
    
    if(verify_checksum_target_neuron != actual_checksum_target_neuron)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to match neuron checksum in bin " 
        << j << ", time slot " 
        << currentTimeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
#if GROUP_EVENTS_VALUES_MODE_V00
    if(verify_checksum_value_01 != actual_checksum_value_01)
    {
      std::cout << "verifyKernelGroupEventsV00: Failed to match value 01 checksum in bin " 
        << j << ", time slot " 
        << currentTimeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
#endif
    /*
    print ? std::cout << verify_checksum_target_neuron << ", ", true : false;
    */
  }
  
  print ? std::cout << std::endl, true : false;
  
  /*Verify event data*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
  print ? std::cout << "Output histogram:" << std::endl, true : false;
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
        //verify_checksum_histogram_out += dataHistogramGroupEventsVerify[p];
        CHECKSUM01(actual_checksum_histogram_out, dataHistogramGroupEventsTik[p]);
        //actual_checksum_histogram_out += dataHistogramGroupEventsTik[p];
        print ? std::cout << p << "," << dataHistogramGroupEventsVerify[p] << ","
          << dataHistogramGroupEventsTik[p] << std::endl, true : false;
      }
      /*Verify*/
      if(verify_checksum_histogram_out != actual_checksum_histogram_out)
      {
        std::cout << "verifyKernelGroupEventsV00: Failed to match output histogram checksum in bin " 
          << b << ", WG " << w << ", time slot " << currentTimeSlot << std::endl;
        result = SDK_FAILURE;
        break;
      }
    }
    if(result != SDK_SUCCESS){break;}
  }
#endif
  }
#endif
  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelScanHistogramV01()
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;
  
#if SCAN_ENABLE_V01
  /* allocate memory for histogram */
  cl_uint *temp = (cl_uint *)calloc(SCAN_HISTOGRAM_BIN_SIZE_V01*SCAN_HISTOGRAM_TOTAL_BINS_V01, 
    sizeof(cl_uint));
  
  /*Compute histogram*/
  for(cl_uint j = 0; j < (SCAN_HISTOGRAM_BIN_SIZE_V01); j++)
  {
    for(cl_uint b = 0; b < (SCAN_HISTOGRAM_TOTAL_BINS_V01); b++)
    {
      cl_uint sum = 0;
      
      for(cl_uint w = 0; w < (SCAN_HISTOGRAM_BIN_BACKETS); w++)
      {
        cl_uint p = 
          /*WG offset*/
          w*(SCAN_HISTOGRAM_BIN_SIZE_V01*SCAN_HISTOGRAM_TOTAL_BINS_V01) + 
          /*WG offset for histogram out*/
          j +
          /*bin*/
          b*SCAN_HISTOGRAM_BIN_SIZE_V01;
          
        sum += dataHistogramGroupEventsTik[p];
      }
      temp[j + b*SCAN_HISTOGRAM_BIN_SIZE_V01] = sum;
    }
  }
  
  cl_uint runningSum = temp[0];
  temp[0] = 0;
  
  /*Compute offsets*/
  for(cl_uint j = 1; j < (SCAN_HISTOGRAM_TOTAL_BINS_V01*SCAN_HISTOGRAM_BIN_SIZE_V01 + 1); j++)
  {
    cl_uint d = temp[j];
    temp[j] = runningSum;
    runningSum += d;
  }
  
  unsigned long long error_count_histogram_out = 0;

  for(cl_uint j = 1; j < (SCAN_HISTOGRAM_TOTAL_BINS_V01*SCAN_HISTOGRAM_BIN_SIZE_V01 + 1); j++)
  {
    error_count_histogram_out += 
      (dataHistogramGroupEventsVerify[j] != temp[j]);
  }

  /*Verify*/
  if(error_count_histogram_out)
  {
    std::cout << "Failed to match 01 histogram " << error_count_histogram_out << " times"
    << std::endl;
    result = SDK_FAILURE;
  }
  
  free(temp);
#endif

  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelGroupEvents
(
  cl_uint verifyKeyBins,
  cl_uint timeSlot,
  cl_uint step,
  cl_uint value1CarryEnable,
  cl_uint value2CarryEnable,
  cl_uint stepShiftEnable,
  cl_uint elementSizeWords,
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
  int result = SDK_SUCCESS;
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
      //verify_checksum_target_neuron += dataGroupedEventsVerify[p];
      CHECKSUM01(actual_checksum_target_neuron, dataGroupedEvents[p]);
      //actual_checksum_target_neuron += dataGroupedEvents[p];
      if(value1CarryEnable)
      {
        CHECKSUM01(verify_checksum_value_01, dataGroupedEventsVerify[destinationBufferSize + p]);
        //verify_checksum_value_01 += dataGroupedEventsVerify[destinationBufferSize + p];
        CHECKSUM01(actual_checksum_value_01, dataGroupedEvents[destinationBufferSize + p]);
        //actual_checksum_value_01 += dataGroupedEvents[destinationBufferSize + p];
      }
      if(value2CarryEnable)
      {
        CHECKSUM01(verify_checksum_value_02, dataGroupedEventsVerify[2*destinationBufferSize + p]);
        //verify_checksum_value_02 += dataGroupedEventsVerify[2*destinationBufferSize + p];
        CHECKSUM01(actual_checksum_value_02, dataGroupedEvents[2*destinationBufferSize + p]);
        //actual_checksum_value_02 += dataGroupedEvents[2*destinationBufferSize + p];
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
      result = SDK_FAILURE;
      break;
    }
    
    if(actual_error_count_target_neuron)
    {
      std::cout << "verifyKernelGroupEvents: Failed to validate correct bin for " 
        << actual_error_count_target_neuron 
        << " keys out of " << (end-start) << " in actual data in bin " << j 
        << ", time slot " << timeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
    }
    
    if(verify_checksum_target_neuron != actual_checksum_target_neuron)
    {
      std::cout << "verifyKernelGroupEvents: Failed to match neuron checksum in bin " 
        << j << ", time slot " << timeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
    if(value1CarryEnable && (verify_checksum_value_01 != actual_checksum_value_01))
    {
      std::cout << "verifyKernelGroupEvents: Failed to match value 01 checksum in bin " 
        << j << ", time slot " 
        << timeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
    if(value2CarryEnable && (verify_checksum_value_02 != actual_checksum_value_02))
    {
      std::cout << "verifyKernelGroupEvents: Failed to match value 02 checksum in bin " 
        << j << ", time slot " 
        << timeSlot << std::endl;
      result = SDK_FAILURE;
      break;
    }
    /*
    print ? std::cout << verify_checksum_target_neuron << ", ", true : false;
    */
  }
  
  if(result != SDK_SUCCESS){return result;}
  
  print ? std::cout << std::endl, true : false;
  
  /*Verify event data*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
        //verify_checksum_histogram_out += dataHistogramOutVerify[p];
        CHECKSUM01(actual_checksum_histogram_out, dataHistogramOut[p]);
        //actual_checksum_histogram_out += dataHistogramOut[p];
      }
      /*Verify*/
      if(verify_checksum_histogram_out != actual_checksum_histogram_out)
      {
        std::cout << "verifyKernelGroupEvents: Failed to match output histogram checksum in bin "
          << b << ", WG " << w << ", time slot " << currentTimeSlot << std::endl;
        result = SDK_FAILURE;
        break;
      }
    }
    if(result != SDK_SUCCESS){break;}
  }
  }
#endif

  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelGroupEventsV01(cl_uint step)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;
  
#if GROUP_EVENTS_ENABLE_V01
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogramGroupEventsTik, 
    dataHistogramGroupEventsTikSizeBytes);
  
  cl_uint event_total = dataHistogramGroupEventsTik[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    
  /*Init data for verification*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
#endif
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
    if(dataOffsetGroupEventsCopy[bin_offset] >= dataHistogramGroupEventsTik[bin_offset+1])
    {
      std::cout << "verifyKernelGroupEventsV01: Destination event bin pointer overlaps "
      << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot << ", WG " 
      << wg_id << ", pointer " << e << ", sort step " << step << std::endl;
      result = SDK_FAILURE; 
      break;
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
      
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
      result = SDK_FAILURE;
      break;
    }
    /*Increment counter for this bin.*/
    dataHistogramGroupEventsVerify[hist_out_ptr]++;
#endif
    
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
  if(result != SDK_SUCCESS){return result;}

  result = verifyKernelGroupEvents
  (
    1,
    currentTimeSlot,
    step,
    GROUP_EVENTS_VALUES_MODE_V01,
    0,
    GROUP_EVENTS_ENABLE_STEP_SHIFT_V01,
    GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_GRID_SIZE_WG,
    dataHistogramGroupEventsTik,
    dataHistogramGroupEventsTok,
    dataHistogramGroupEventsVerify,
    dataGroupEventsTok,
    dataGroupEventsTikVerify
  );
  
#endif
  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelGroupEventsV02(cl_uint step, cl_uint keyOffset)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;
  
#if GROUP_EVENTS_ENABLE_V02
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogramGroupEventsTik, 
    dataHistogramGroupEventsTikSizeBytes);
    
  cl_uint event_total = dataHistogramGroupEventsTik[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE];
    
  /*Init data for verification*/
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
#endif
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
    cl_uint new_key = dataUnsortedEventTargets[value+keyOffset];
    
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
    if(dataOffsetGroupEventsCopy[bin_offset] >= dataHistogramGroupEventsTik[bin_offset+1])
    {
      std::cout << "verifyKernelGroupEventsV02: Destination event bin pointer overlaps "
      << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot << ", WG " 
      << wg_id << ", pointer " << e << std::endl;
      result = SDK_FAILURE; 
      break;
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
      
#if (GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT)
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
      result = SDK_FAILURE;
      break;
    }
    /*Increment counter for this bin.*/
    dataHistogramGroupEventsVerify[hist_out_ptr]++;
#endif
    
    /*Store event at its group location (grouped by bins)*/
    dataGroupEventsTikVerify[dest_offset] = new_key;
    dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = value;

    /*Increment ptr for next data item*/
    dataOffsetGroupEventsCopy[bin_offset]++;
  }
  
  free(dataOffsetGroupEventsCopy);
  if(result != SDK_SUCCESS){return result;}

  result = verifyKernelGroupEvents
  (
    0,
    currentTimeSlot,
    step,
    GROUP_EVENTS_VALUES_MODE_V02,
    0,
    GROUP_EVENTS_ENABLE_STEP_SHIFT_V02,
    GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_GRID_SIZE_WG,
    dataHistogramGroupEventsTik,
    dataHistogramGroupEventsTok,
    dataHistogramGroupEventsVerify,
    dataGroupEventsTok,
    dataGroupEventsTikVerify
  );
  
#endif

  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelGroupEventsV03(cl_uint step)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;
  
#if GROUP_EVENTS_ENABLE_V03
  /* allocate memory for offset copies used in verification data generation*/
  cl_uint size = GROUP_EVENTS_GRID_SIZE_WG*GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE*
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
  cl_uint *dataOffsetGroupEventsCopy = (cl_uint *)calloc(size, sizeof(cl_uint));
  memcpy(dataOffsetGroupEventsCopy, dataHistogramGroupEventsTik, 
    dataHistogramGroupEventsTikSizeBytes);
    
  cl_uint event_total = dataHistogramGroupEventsTik[GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*
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
    if(dataOffsetGroupEventsCopy[bin_offset] >= dataHistogramGroupEventsTik[bin_offset+1])
    {
      std::cout << "verifyKernelGroupEventsV03: Destination event bin pointer overlaps "
      << "with next pointer for bin " << bin << ", time slot " << currentTimeSlot << ", WG " 
      << wg_id << ", pointer " << e << std::endl;
      result = SDK_FAILURE; 
      break;
    }
    
    /*Calculate offset in the grouped data space*/
    cl_uint dest_offset = dataOffsetGroupEventsCopy[bin_offset];
    
    /*Store event at its group location (grouped by bins)*/
    dataGroupEventsTikVerify[dest_offset] = key;
    dataGroupEventsTikVerify[GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = 
      dataUnsortedEventDelays[valuePtr];
    dataGroupEventsTikVerify[2*GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE + dest_offset] = 
      dataUnsortedEventWeights[valuePtr];
    
    /*Increment ptr for next data item*/
    dataOffsetGroupEventsCopy[bin_offset]++;
  }
  
  free(dataOffsetGroupEventsCopy);
  if(result != SDK_SUCCESS){return result;}

  result = verifyKernelGroupEvents
  (
    1,
    currentTimeSlot,
    step,
    GROUP_EVENTS_VALUES_MODE_V03,
    1,
    GROUP_EVENTS_ENABLE_STEP_SHIFT_V03,
    GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
    GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS,
    GROUP_EVENTS_HISTOGRAM_BIN_MASK,
    GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03,
    GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT,
    GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE,
    GROUP_EVENTS_GRID_SIZE_WG,
    dataHistogramGroupEventsTik,
    dataHistogramGroupEventsTok,
    NULL,
    dataGroupEventsTok,
    dataGroupEventsTikVerify
  );
  
#endif

  return result;
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelMakeEventPtrs()
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;
  
#if MAKE_EVENT_PTRS_ENABLE
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

  if(result == SDK_SUCCESS)
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
          result = SDK_FAILURE;
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
          result = SDK_FAILURE;
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
            result = SDK_FAILURE;
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
      result = SDK_FAILURE;
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
  
  if(!error && (result != SDK_FAILURE))
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
            result = SDK_FAILURE;
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
    result = SDK_FAILURE;
  }
  
  /*Verify total count match*/
  if(result != SDK_FAILURE && totalCountsVerifyStruct != totalCountsVerifyOriginal)
  {
    std::cerr << "ERROR, verifyKernelMakeEventPtrs, sum of pointer counters in all structs "
      << "doesn't match total in original data: " << totalCountsVerifyStruct << " vs " 
      << totalCountsVerifyOriginal << std::endl;
    result = SDK_FAILURE;
  }

  free(wfDataIter);
#endif
#endif

  return result;
}
/**************************************************************************************************/



/*
  Inject sorted synaptic events into event queue
*/
int 
IntegrationTest::injectSortedEvents
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
  int result = SDK_SUCCESS;

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
      result = SDK_FAILURE;
      break;
    }
    
    if(verify && nrn[nId].n_in != 0)
    {
      std::cerr << "ERROR, injectSortedEvents, detected another entry for  " 
        << nId << " with count " << nrn[nId].n_in << std::endl;
      result = SDK_FAILURE;
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
        result = SDK_FAILURE;
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
          result = SDK_FAILURE;
          break;
        }
        
        if(t > 1.0 || t < 0.0)
        {
          std::cerr << "ERROR, injectSortedEvents, detected event time outside of bounds for neuron "
            << nId << ": " << t << " is not within (0.0, 1.0)" << std::endl;
          result = SDK_FAILURE;
          break;
        }
      }
      
      nrn[nId].in_t[j] = t;
      nrn[nId].in_w[j] = w;
    }
    if(result != SDK_SUCCESS){break;}
  }

  return result;
}
/**************************************************************************************************/



/*
  Inject unsorted synaptic events into event queue
*/
int 
IntegrationTest::injectUnsortedEvents
(
  cl_uint       timeSlots,
  cl_uint       buffers,
  cl_uint       bufferSize,
  cl_uint       totalNeurons,
  cl_uint       eventQueueSize,
  cl_uint       *dataUnsortedEventCounts,
  cl_uint       *dataUnsortedEventTargets,
  cl_uint       *dataUnsortedEventWeights,
  cl_uint       *dataUnsortedEventDelays,
  neuron_iz_ps  *nrn
)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;

  for(cl_uint s = 0; s < timeSlots; s++)
  {
    for(cl_uint b = 0; b < buffers; b++)
    {
      int total = dataUnsortedEventCounts[timeSlots*b + s];
      
      cl_uint ptr = 
        /*Event data buffers*/
        b * timeSlots * 
        bufferSize +
        /*Current event data buffer*/
        s * bufferSize;

      for(int e = 0; e < total; e++)
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
          result = SDK_FAILURE;
          break;
        }
      }
      if(result != SDK_SUCCESS){break;}
    }
    if(result != SDK_SUCCESS){break;}
  }
  
#if 0//FLIXIBLE_DELAYS_ENABLE
  if(result == SDK_SUCCESS)
  {
    /*Decrement synaptic events: */
    for(unsigned int i = 0; i < totalNeurons; i++)
    {
      for(int j = 0; j < nrn[i].n_in; j++)
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
  Propagation of spike events to synaptic events for PS method
*/
int 
IntegrationTest::propagateSpikes
(
  unsigned int  totalNeurons,
  unsigned int  eventQueueSize,
  neuron_iz_ps  *nrn,
  int           *ne,
  DATA_TYPE     *te,
  unsigned int  *synapsePointer,
  unsigned int  *synapseTargets,
  DATA_TYPE     *synapseDelays,
  DATA_TYPE     *synapseWeights
)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;

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
      unsigned int ptrStart = synapsePointer[i];
      unsigned int ptrEnd = synapsePointer[i+1];

      /*Iterate over target neurons of this source neuron*/
      for(unsigned int s = ptrStart; s < ptrEnd; s++)
      {
        /*Get target neuron and its event count*/
        unsigned int k = synapseTargets[s];
        unsigned int n_in = nrn[k].n_in;
        DATA_TYPE t_event = 0.0;
        
        if(n_in<eventQueueSize)
        {
	        unsigned int j=n_in; 
#if FLIXIBLE_DELAYS_ENABLE
          t_event = te[i] + synapseDelays[s];
#else
          t_event = te[i];
#endif
          /*Use insertion sort to maintain ordered synaptic events*/
          while ((j > 0) && (nrn[k].in_t[j-1] > t_event))
          {
	          nrn[k].in_t[j] = nrn[k].in_t[j-1]; /*shift*/
	          nrn[k].in_w[j] = nrn[k].in_w[j-1]; 
            j--;
	        }
          
          nrn[k].in_t[j] = t_event;
          nrn[k].in_w[j] = synapseWeights[s];
          nrn[k].n_in++;
        }
        else
        {
          std::cerr << "ERROR, propagateSpikes, detected event queue overflow for nID " 
            << k << ": " << n_in << " >= " << eventQueueSize << std::endl;
          result = SDK_FAILURE;
          break;
        }
      }
      if(result != SDK_SUCCESS){break;}
    }
  }
  
#if FLIXIBLE_DELAYS_ENABLE
  if(result == SDK_SUCCESS)
  {
    /*Decrement synaptic events: */
    for(unsigned int i = 0; i < totalNeurons; i++)
    {
      for(int j = 0; j < nrn[i].n_in; j++)
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
IntegrationTest::verifyEvents
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
  int result = SDK_SUCCESS;

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
        result = SDK_FAILURE;
        break;
      }
      
      unsigned int nD = sortedEvents[entryAddress + e];
      DATA_TYPE tD = *((DATA_TYPE *)(&sortedEvents[entryAddress + e + sortedEventsSize]));
      DATA_TYPE wD = *((DATA_TYPE *)(&sortedEvents[entryAddress + e + 2*sortedEventsSize]));
      
      /*Verify neuron IDs*/
      if(nH != nD)
      {
        std::cerr << "ERROR, verifyEvents, neuron ID mismatch: " << nH << " != " << nD << std::endl;
        result = SDK_FAILURE;
        break;
      }
      /*Verify event time*/
      if(tH != tD)
      {
        std::cerr << "ERROR, verifyEvents, event time mismatch for neuron " << nH << ": " 
          << tH << " != " << tD << " ("; PRINT_HEX(4, tH); std::cerr << " != "; PRINT_HEX(4, tD);
          std::cerr << ")" << std::endl;
        result = SDK_FAILURE;
        break;
      }
      if(result == SDK_FAILURE){break;}
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
              result = SDK_FAILURE;
              break;
            }
          }

          free(weightCheck);
          
          if(result != SDK_FAILURE)
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
          result = SDK_FAILURE;
          break;
        }
      }
    }
    if(result == SDK_FAILURE){break;}
    
    if(result == SDK_SUCCESS && e != entryCount)
    {
      std::cerr << "ERROR, verifyEvents, event count mismatch 2 for neuron " << nH << ": " 
        << e << " != " << entryCount << std::endl;
      result = SDK_FAILURE;
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
IntegrationTest::stepIzPs
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
  int *fcount,
  unsigned long long int *icount, 
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
	int result = SDK_SUCCESS;
  
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
    result = SDK_FAILURE;
  }
  
	if(ps_order >= ps_order_limit)
  {
#if PRINT_stepIzPs
    std::cerr << "WARNING, stepIzPs, PS order limit overflow: " << ps_order << std::endl;
#endif
    result = SDK_FAILURE;
  }

#if STATISTICS_ENABLE
  (*icount)++;
	*mu += ((DATA_TYPE)ps_order - *mu)/(DATA_TYPE)*icount;
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
      result = SDK_FAILURE;
    }

		if(dt_part>dt_full || dt_part<0)
    {
      dt_part=dt_full/2;
#if PRINT_stepIzPs
      std::cerr << "WARNING, stepIzPs, detected NR divergence for neuron ID: " << nrn_ind 
        << std::endl;
#endif
      result = SDK_FAILURE;
    }
    
#if STATISTICS_ENABLE
		/*Record spike and schedule events*/
    (*fcount)++;
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
      result = SDK_FAILURE;
    }
    
    if(ps_order >= ps_order_limit)
    {
#if PRINT_stepIzPs
      std::cerr << "WARNING, stepIzPs, PS order limit overflow: " << ps_order << std::endl;
#endif
      result = SDK_FAILURE;
    }

#if STATISTICS_ENABLE
    (*icount)++; 
		*mu += ((DATA_TYPE)ps_order- *mu)/(DATA_TYPE)*icount;
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



int 
IntegrationTest::updateStep
(
  bool    ignoreFailures,
  cl_uint injectCurrentUntilStep,
  cl_uint currentTimeStep,
  cl_uint totalNeurons,
  cl_uint psOrderLimit,
  cl_uint nrOrderLimit,
  double  nrTolerance
)
/**************************************************************************************************/
{
#define PRINT_updateStep                                                           STATISTICS_ENABLE
  int result = SDK_SUCCESS;
  int nrn_ind, substep;

#if STATISTICS_ENABLE
  int fcount = 0, max_order_ps = 0, max_events = 0;
  unsigned long long int icount = 0;
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
  
  if((injectCurrentUntilStep > 0) && (currentTimeStep == injectCurrentUntilStep))
  {
    for(nrnp_ps=nrn_ps;nrnp_ps<nrnx_ps;nrnp_ps++)nrnp_ps->I=0;
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
      substep = 0;
      DATA_TYPE event = nrnp_ps->in_t[substep];

      while((event <= 1.0f) && (substep<(nrnp_ps->n_in)))
#else
      for(substep=0;substep<(nrnp_ps->n_in);substep++)
#endif
      {
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
            &fcount,
            &icount,
            &mu_order_ps,
            &max_order_ps,
#endif
            psOrderLimit,
            nrOrderLimit,
            (DATA_TYPE)nrTolerance,
            FLIXIBLE_DELAYS_ENABLE
          );
          if(result != SDK_SUCCESS && !ignoreFailures){break;}

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
      if(result != SDK_SUCCESS && !ignoreFailures){break;}
      
#if(FLIXIBLE_DELAYS_ENABLE)
      /*Shift future synaptic events: */
      nrnp_ps->n_in -= substep;
      
      if(nrnp_ps->n_in > 0)
      {
        int ptr1 = 0, ptr2 = substep; 
        substep = nrnp_ps->n_in;
        
        while (substep > 0)
        {
          nrnp_ps->in_t[ptr1] = nrnp_ps->in_t[ptr2];
          nrnp_ps->in_w[ptr1] = nrnp_ps->in_w[ptr2];
          ptr1++; ptr2++; substep--;
        }
      }
      else if(nrnp_ps->n_in < 0)
      {
        std::cerr << "ERROR updateStep: negative overflow in synaptic event number." << std::endl;
        nrnp_ps->n_in = 0;
        result = SDK_FAILURE;
        if(!ignoreFailures){break;}
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
        &fcount,
        &icount,
        &mu_order_ps,
        &max_order_ps,
#endif
        psOrderLimit,
        nrOrderLimit,
        (DATA_TYPE)nrTolerance,
        FLIXIBLE_DELAYS_ENABLE
      );
      if(result != SDK_SUCCESS && !ignoreFailures){break;}
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
        &fcount,
        &icount,
        &mu_order_ps,
        &max_order_ps,
#endif
        psOrderLimit,
        nrOrderLimit,
        (DATA_TYPE)nrTolerance,
        FLIXIBLE_DELAYS_ENABLE
      );
      if(result != SDK_SUCCESS && !ignoreFailures){break;}
    }
  }/*For neurons*/
  
#if PRINT_updateStep && STATISTICS_ENABLE
  std::cout << "\n Spikes: " 
            << "\n  Total: " << fcount 
            << "\n  Average: " << ((double)fcount/(double)totalNeurons)
            << "\n Synaptic events: " 
            << "\n  Max per neuron: " << max_events 
            << "\n ps_step calls: " << icount 
            << "\n Average PS order: " << mu_order_ps 
            << "\n Max PS order: " << max_order_ps 
            << std::endl;
#endif

  free(yold); free(ynew);
	for(int i=0;i<5;i++){free(yp[i]);} free(yp);

  return result;
#undef PRINT_updateStep
}
/**************************************************************************************************/



int 
IntegrationTest::verifyKernelUpdateNeurons
(
  bool          verify,
  cl_uint       step,
  cl_uint       sortedEventsSize,
  unsigned int  *pointerStruct,
  unsigned int  *sortedEvents,
  unsigned int  *synapsePointer,
  unsigned int  *synapseTargets,
  DATA_TYPE     *synapseDelays,
  DATA_TYPE     *synapseWeights,
  unsigned int  *spikePackets,
  DATA_TYPE     *modelVariables,
  unsigned int  *dataSpikePackets
)
/**************************************************************************************************/
{
  int result = SDK_SUCCESS;
#if UPDATE_NEURONS_ENABLE_V00
  
  bool breakOnFailure = 1;
  bool ignoreSolverFailuresHost = IGNORE_SOLVER_EXCEPTIONS;

  if((synapsePointer != NULL) && (synapseTargets != NULL) && (synapseDelays != NULL) &&
    (synapseWeights != NULL))
  {
    if(spikePackets != NULL)
    {
      memset(ne, 0, UPDATE_NEURONS_TOTAL_NEURONS*sizeof(int));
      
      /*Iterate through spike packets*/
      for(cl_uint packet = 0; packet < UPDATE_NEURONS_SPIKE_PACKETS_V00; packet++)
      {
        cl_uint packet_index = packet * UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS;
        cl_uint total_spikes = spikePackets[packet_index];

        /*Iterate through spikes in a current packet and inject spikes*/
        for(cl_uint i = 0; i < total_spikes; i++)
        {
          cl_uint spiked_neuron = spikePackets[packet_index + 
            UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + 
            UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS * i];
          cl_float spike_time = *((cl_float *)(&spikePackets[packet_index + 
            UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + 
            UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS * i + 1]));
            
          if(ne[spiked_neuron] == 1)
          {
            std::cerr << "ERROR, verifyKernelUpdateNeurons, duplicate entry detected for neuron ID " 
              << spiked_neuron << " while injecting its spike: 1)" << spike_time << ", 2)" 
              << te_ps[spiked_neuron] << std::endl;
            return SDK_FAILURE;
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
      synapsePointer,
      synapseTargets,
      synapseDelays,
      synapseWeights
    );
    if(result != SDK_SUCCESS){return result;}
    
#if NETWORK_VERIFY_ENABLE
    if((pointerStruct != NULL) && (sortedEvents != NULL) && FLIXIBLE_DELAYS_ENABLE && verify)
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
      if(result != SDK_SUCCESS){return result;}
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
    if(result != SDK_SUCCESS){return result;}
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

  if(result != SDK_SUCCESS && !ignoreSolverFailuresHost){return result;}

#if NETWORK_VERIFY_ENABLE
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
      result = SDK_FAILURE;
    }
    
    underTestType1 = modelVariables[UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].u;
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable u for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = SDK_FAILURE;
    }
    
    underTestType1 = modelVariables[2*UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].g_ampa;
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable g_ampa for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = SDK_FAILURE;
    }
    
    underTestType1 = modelVariables[3*UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].g_gaba;
    if(underTestType1 != reference)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable g_gaba for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = SDK_FAILURE;
    }
    
#if (LOG_MODEL_VARIABLES)
    if(i == LOG_MODEL_VARIABLES_NEURON_ID)
    {
      *dataToTraceFile LOG_MODEL_VARIABLES_FILE_BODY(LOG_MODEL_VARIABLES_NEURON_ID);
    }
#endif
    /*
    cl_uint underTestType2 = pointerStruct[i*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE + 1];
    if(underTestType2 != 0)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, event count was not reset for neuron " << i
        << ": " << underTestType2 << std::endl;
      result = SDK_FAILURE;
    }
    */
    if(result != SDK_SUCCESS && breakOnFailure){break;}
  }

  /*Verify spikes*/
  if(result != SDK_FAILURE)
  {
    char *spikeCheck = (char *)calloc(UPDATE_NEURONS_TOTAL_NEURONS, sizeof(char));
    
    /*All verified spikes have equivalents in the reference data*/
    for(cl_uint packet = 0; packet < UPDATE_NEURONS_SPIKE_PACKETS_V00; packet++)
    {
      cl_uint packet_index = packet * UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS;
      cl_uint total_spikes = dataSpikePackets[packet_index];

      /*Iterate through spikes in a current packet*/
      for(cl_uint i = 0; i < total_spikes; i++)
      {
        cl_uint spiked_neuron = dataSpikePackets[packet_index + 
          UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS * i];
        cl_float spike_time = *((cl_float *)(&dataSpikePackets[packet_index + 
          UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS * i 
          + 1]));
        if(te_ps[spiked_neuron] != spike_time)
        {
          std::cerr << "ERROR, verifyKernelUpdateNeurons, spike time mismatch for neuron " 
            << spiked_neuron << ", packet " << packet << ": " << te_ps[spiked_neuron] << "!=" 
            << spike_time << std::endl;
          result = SDK_FAILURE;
          if(breakOnFailure){break;}
        }
        else
        {
          spikeCheck[spiked_neuron] = 1;
        }
      }
      if(breakOnFailure && (result == SDK_FAILURE)){break;}
    }
    
    /*There are no reference spikes, which are not present in the verified data*/
    if(result != SDK_FAILURE)
    {
      for(cl_uint n = 0; n < UPDATE_NEURONS_TOTAL_NEURONS; n++)
      {
        if(te_ps[n] != 0.0 && spikeCheck[n] == 0)
        {
          std::cerr << "ERROR, verifyKernelUpdateNeurons, a spike from neuron " << n 
            << " and spike time " << te_ps[n] << " is absent"<< std::endl;
          result = SDK_FAILURE;
          if(breakOnFailure){break;}
        }
      }
    }
    
    free(spikeCheck);
  }
  }
#endif
#endif

  return result;
}
/**************************************************************************************************/



//TODO: it could be useful to replace kernel-specific defines with generic ones in this method:
int 
IntegrationTest::verifySortedEvents
(
  cl_uint *sortedEvents, 
  cl_uint *pointerStruct, 
  cl_uint level
)
/**************************************************************************************************/
{
#if SORT_VERIFY_ENABLE
#define PRINT_VERIFY_SORTED_EVENTS  0
  int result = SDK_SUCCESS;
  cl_uint keyOffset = 0;
  cl_uint val1Offset = (keyOffset+1)%GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  cl_uint val2Offset = (keyOffset+2)%GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
  char *sortedEventsCheck;
  
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
      result = SDK_FAILURE;
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
        result = SDK_FAILURE;
        break;
      }
      
      if(time1 < time2)
      {
        std::cout << "verifySortedEvents " << currentTimeStep 
          << ": failed to verify sort order for value element " << i << ": " << time1 << "<" 
          << time2 << std::endl;
        result = SDK_FAILURE;
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
  if(result != SDK_FAILURE)
  {
    std::cout << "verifySortedEvents " << currentTimeStep << ": verified sort order" << std::endl;
  }
#endif

  if(result != SDK_FAILURE && level > 0)
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
        result = SDK_FAILURE;
        break;
      }
#endif

      /*Find the entry*/
      result = SDK_FAILURE;
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
              result = SDK_SUCCESS;
              sortedEventsCheck[a] = 1;
              break;
            }
          }
        }
      }

      if(result == SDK_FAILURE)
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
    if(result != SDK_FAILURE)
    {
      std::cout 
        << "\nverifySortedEvents" << currentTimeStep << ": verified mapping of unsorted "
        << "data in sorted data" << std::endl;
    }
#endif
  }
  
  /*Check for unrecognized keys*/
  if(result != SDK_FAILURE && level > 1)
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
        //result = SDK_FAILURE;
        //break;
      }
    }
#if PRINT_VERIFY_SORTED_EVENTS
    if(result != SDK_FAILURE)
    {
      std::cout << "verifySortedEvents " << currentTimeStep << ": checked for unrecognized keys" 
        << std::endl;
    }
#endif
  }
  
#if PRINT_VERIFY_SORTED_EVENTS
  std::cout << "verifySortedEvents " << currentTimeStep << ": finished" << std::endl;
#endif

  free(dataUnsortedEventsSnapShot); dataUnsortedEventsSnapShot = NULL;
  if(level > 0)
  {
    free(sortedEventsCheck);
  }
  return result;
#else
  std::cerr << "ERROR, verifySortedEvents: " 
    << "(SORT_VERIFY_ENABLE is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



int 
IntegrationTest::captureUnsortedEvents
(
  cl_uint *unsortedEventCounts,
  cl_uint *unsortedEventTargets,
  cl_uint *unsortedEventDelays,
  cl_uint *unsortedEventWeights
)
/**************************************************************************************************/
{
#if SORT_VERIFY_ENABLE
  cl_uint ptr_store = 0;
  cl_uint event_total = 0;
  int result = SDK_SUCCESS;
  
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
        result = SDK_FAILURE;
        break;
      }
      
      dataUnsortedEventsSnapShot[ptr_store] = neuron;
      dataUnsortedEventsSnapShot[ptr_store+1] = time;
      dataUnsortedEventsSnapShot[ptr_store+2] = weight;
      ptr_store += GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
    }
  }
  
  if(result != SDK_SUCCESS)
  {
    free(dataUnsortedEventsSnapShot);
  }
  
  return result;
#else
  std::cerr << "ERROR, captureUnsortedEvents: " 
    << "(SORT_VERIFY_ENABLE is not true)" << std::endl;
  return SDK_FAILURE;
#endif
}
/**************************************************************************************************/



#if SIMULATION_SNAPSHOT
int 
IntegrationTest::takeSimulationSnapshot
(
  cl_uint   step,
  cl_uint   sampleSizeNeurons,
  cl_uint   *dataUnsortedEventCounts,
  cl_uint   *dataUnsortedEventWeights,
  cl_uint   *dataMakeEventPtrsStruct,
  cl_uint   *dataGroupEventsTik,
  cl_uint   *dataSpikePackets,
  cl_float  *modelVariables
)
/**************************************************************************************************/
{
  SET_RANDOM_SEED(srandSeed);
  LOG("takeSimulationSnapshot: set srand seed to " << srandSeed, 0);
  (*dataToSnapshotLogFile).str("");
  *dataToSnapshotLogFile << "\n\nSNAPSHOT AT STEP: " << step << "\n\n";
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE && GROUP_EVENTS_ENABLE_V03
  {
    *dataToSnapshotLogFile << "\n\nEVENT BUFFERS\n\n";
    *dataToSnapshotLogFile << "Time Slot,Parameter Name,Parameter Value" << std::endl;
      
    for(cl_uint s = 0; s < GROUP_EVENTS_TIME_SLOTS; s++)
    {
      double totalEventsMean = 0, totalEventsVariance = 0, percentInh = 0, percentExc = 0;
      int totalEvents = 0, totalEventsMax = 0, 
        totalEventsMin = GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE;
      
      for(cl_uint b = 0; b < GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS; b++)
      {
        int total = dataUnsortedEventCounts[GROUP_EVENTS_TIME_SLOTS*b + s];
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
          cl_float w = *((cl_float *)(&dataUnsortedEventWeights[ptr + e]));
          if(w < 0){percentInh++;}
          else{percentExc++;}
        }
      }
      totalEventsVariance /= GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS;
      totalEventsVariance = sqrt(totalEventsVariance);
      percentInh = 100.0*percentInh/totalEvents;
      percentExc = 100.0*percentExc/totalEvents;
      
      *dataToSnapshotLogFile
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

  *dataToSnapshotLogFile << "\n\nEVENTS\n\n";
  *dataToSnapshotLogFile << "Parameter Name,Parameter Value" << std::endl;

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
  
  *dataToSnapshotLogFile << "Total Neurons," << UPDATE_NEURONS_TOTAL_NEURONS << std::endl
    << "Total Inhibitory Events," << totalEventsInh << std::endl
    << "Total Excitatory Events," << totalEventsExc << std::endl
    << "Events Per Neuron Max," << eventsPerNeuronMax << std::endl
    << "Events Per Neuron Min," << eventsPerNeuronMin << std::endl;
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE && UPDATE_NEURONS_ENABLE_V00
  *dataToSnapshotLogFile << "\n\nSPIKES\n\n";
  *dataToSnapshotLogFile << "Parameter Name,Parameter Value" << std::endl;
  
  cl_uint totalSpikes = 0, spikesPerPacketMax = 0, 
    spikesPerPacketMin = UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE;
  cl_float spikeTimeMax = 0.0f, spikeTimeMin = 100.0f;
  
  /*Iterate through spike packets*/
  for(cl_uint packet = 0; packet < UPDATE_NEURONS_SPIKE_PACKETS_V00; packet++)
  {
    cl_uint packetIndex = packet * UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS;
    cl_uint spikes = dataSpikePackets[packetIndex];
    
    totalSpikes += spikes;
    if(spikesPerPacketMax < spikes){spikesPerPacketMax = spikes;}
    if(spikesPerPacketMin > spikes){spikesPerPacketMin = spikes;}

    /*Iterate through spikes in a current packet*/
    for(cl_uint i = 0; i < spikes; i++)
    {
      cl_uint spiked_neuron = dataSpikePackets[packetIndex + 
        UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + 
        UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS * i];
      cl_float spike_time = *((cl_float *)(&dataSpikePackets[packetIndex + 
        UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + 
        UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS * i + 1]));

      if(spikeTimeMax < spike_time){spikeTimeMax = spike_time;}
      if(spikeTimeMin > spike_time){spikeTimeMin = spike_time;}
    }
  }
  
  *dataToSnapshotLogFile << "Total Spikes," << totalSpikes << std::endl
    << "Packet Spikes Max," << spikesPerPacketMax << std::endl
    << "Packet Spikes Min," << spikesPerPacketMin << std::endl
    << "Spike Time Max," << spikeTimeMax << std::endl
    << "Spike Time Min," << spikeTimeMin << std::endl;

#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
  *dataToSnapshotLogFile << "\n\nNEURON MODEL VARIABLES\n\n";
  *dataToSnapshotLogFile << "Neuron ID,v,u,g_ampa,g_gaba" << std::endl;
  cl_uint window = UPDATE_NEURONS_TOTAL_NEURONS/sampleSizeNeurons;
  
  for(cl_uint i = 0; i < sampleSizeNeurons; i++)
  {
    cl_uint nId = i*window + cl_uint(abs((window-1)*((double)rand()/((double)RAND_MAX))));
          
    *dataToSnapshotLogFile << nId << "," << modelVariables[nId] << "," 
      << modelVariables[UPDATE_NEURONS_TOTAL_NEURONS+nId] << "," 
      << modelVariables[2*UPDATE_NEURONS_TOTAL_NEURONS+nId] << "," 
      << modelVariables[3*UPDATE_NEURONS_TOTAL_NEURONS+nId] << std::endl;
  }
#endif
/**************************************************************************************************/

  *snapshotLogFile << (*dataToSnapshotLogFile).str();
  (*dataToSnapshotLogFile).str("");
  
  return SDK_SUCCESS;
}
#endif
/**************************************************************************************************/



void 
IntegrationTest::psClean()
/**************************************************************************************************/
{
  for(int i=0;i<4;i++){free(co[i]);} free(co); free(nrn_ps); free(te_ps); free(ne);
}
/**************************************************************************************************/



int 
IntegrationTest::cleanup()
{
/**************************************************************************************************/
#if ((GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) || SCAN_ENABLE_V01 ||\
     (GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT) ||\
     (GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT))
  if(dataHistogramGroupEventsTik)
      free(dataHistogramGroupEventsTik);
  if(dataHistogramGroupEventsVerify)
      free(dataHistogramGroupEventsVerify);
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V02 || \
    GROUP_EVENTS_ENABLE_V03
  if(dataUnsortedEventCounts)
      free(dataUnsortedEventCounts);
  if(dataUnsortedEventTargets)
      free(dataUnsortedEventTargets);
  if(dataUnsortedEventDelays)
      free(dataUnsortedEventDelays);
  if(dataUnsortedEventWeights)
      free(dataUnsortedEventWeights);
#endif
/**************************************************************************************************/
#if ((EXPAND_EVENTS_ENABLE && EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM) || SCAN_ENABLE_V00\
     || GROUP_EVENTS_ENABLE_V00)
  if(dataHistogram)
      free(dataHistogram);
  if(dataHistogramVerify)
      free(dataHistogramVerify);
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03 || MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  if(dataGroupEventsTik)
      free(dataGroupEventsTik);
  if(dataGroupEventsTikVerify)
      free(dataGroupEventsTikVerify);
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
#if (GROUP_EVENTS_DEBUG_ENABLE)
  if(dataDebugHostGroupEvents)
      free(dataDebugHostGroupEvents);
  if(dataDebugDeviceGroupEvents)
      free(dataDebugDeviceGroupEvents);
#endif
#if (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  if(dataErrorGroupEvents)
      free(dataErrorGroupEvents);
#endif
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  if(dataSpikePackets) 
      free(dataSpikePackets);
#endif
/**************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
  if(dataSynapsePointer)
      free(dataSynapsePointer);
  if(dataSynapseTargets)
      free(dataSynapseTargets);
  if(dataSynapseDelays)
      free(dataSynapseDelays);
  if(dataSynapseWeights)
      free(dataSynapseWeights);
  if(dataUnsortedEventCountsVerify)
      free(dataUnsortedEventCountsVerify);
  if(dataUnsortedEventsTargetsVerify)
      free(dataUnsortedEventsTargetsVerify);
  if(dataUnsortedEventsDelaysVerify)
      free(dataUnsortedEventsDelaysVerify);
  if(dataUnsortedEventsWeightsVerify)
      free(dataUnsortedEventsWeightsVerify);
#if (EXPAND_EVENTS_DEBUG_ENABLE)
  if(dataExpandEventsDebugHost)
      free(dataExpandEventsDebugHost);
  if(dataExpandEventsDebugDevice)
      free(dataExpandEventsDebugDevice);
#endif
#if (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  if(dataExpandEventsError)
      free(dataExpandEventsError);
#endif
#endif
/**************************************************************************************************/
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
#if (SCAN_DEBUG_ENABLE)
  if(dataScanDebugHost)
      free(dataScanDebugHost);
  if(dataScanDebugDevice)
      free(dataScanDebugDevice);
#endif
#if (SCAN_ERROR_TRACK_ENABLE)
  if(dataScanError)
      free(dataScanError);
#endif
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03 ||\
    MAKE_EVENT_PTRS_ENABLE
  if(dataHistogramGroupEventsTok)
      free(dataHistogramGroupEventsTok);
#endif
/**************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  if(dataGroupEventsTok)
      free(dataGroupEventsTok);
  if(dataGroupEventsTokVerify)
      free(dataGroupEventsTokVerify);
#endif
/**************************************************************************************************/
#if SORT_VERIFY_ENABLE
  if(dataUnsortedEventsSnapShot)
      free(dataUnsortedEventsSnapShot);
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  if(dataMakeEventPtrsStruct)
      free(dataMakeEventPtrsStruct);
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  if(dataMakeEventPtrsDebugHost)
      free(dataMakeEventPtrsDebugHost);
  if(dataMakeEventPtrsDebugDevice)
      free(dataMakeEventPtrsDebugDevice);
#endif
#if (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  if(dataMakeEventPtrsError)
      free(dataMakeEventPtrsError);
#endif
#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
#if (UPDATE_NEURONS_DEBUG_ENABLE)
  if(dataUpdateNeuronsDebugHost)
      free(dataUpdateNeuronsDebugHost);
  if(dataUpdateNeuronsDebugDevice)
      free(dataUpdateNeuronsDebugDevice);
#endif
#if (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  if(dataUpdateNeuronsError)
      free(dataUpdateNeuronsError);
#endif
#if (UPDATE_NEURONS_TOLERANCE_MODE > 1)
  if(psTolerance)
      free(psTolerance);
#endif
  if(modelParameters)
      free(modelParameters);
  if(modelVariables)
      free(modelVariables);
  if(constantCoefficients)
      free(constantCoefficients);
#if (LOG_MODEL_VARIABLES)
  *traceFile << (*dataToTraceFile).str();
  (*traceFile).close();
#endif
  psClean();
#endif
#if (LOG_SIMULATION)
  *simulationLogFile << (*dataToSimulationLogFile).str();
  (*simulationLogFile).close();
#endif
#if (LOG_REPORT)
  *reportLogFile << (*dataToReportLogFile).str();
  (*reportLogFile).close();
#endif
#if (SIMULATION_SNAPSHOT)
  (*snapshotLogFile).close();
#endif
/**************************************************************************************************/

  return SDK_SUCCESS;
}
/**************************************************************************************************/



void 
IntegrationTest::printStats()
{
  std::string strArray[4] = 
  {
      "Time(sec)", 
      "[Transfer+Kernel]Time(sec)"
  };
  std::string stats[2];

  totalTime = setupTime + kernelTime;

  stats[0] = sampleCommon->toString(totalTime, std::dec);
  stats[1] = sampleCommon->toString(kernelTime, std::dec);

  this->SDKSample::printStats(strArray, stats, 4);
}



int 
main(int argc, char *argv[])
{
  IntegrationTest test("OpenCL IntegrationTest");

  if(test.initialize() !=  SDK_SUCCESS)
    return SDK_FAILURE;

  if(!test.parseCommandLine(argc, argv))
    return SDK_FAILURE;

  if(test.getPlatformStats() !=  SDK_SUCCESS)
    {test.cleanup(); std::cout << "\nRESULT:FAIL\n"; return SDK_FAILURE;}
  if(test.setup() !=  SDK_SUCCESS)
    {test.cleanup(); std::cout << "\nRESULT:FAIL\n"; return SDK_FAILURE;}
  if(test.run() !=  SDK_SUCCESS)
    {test.cleanup(); std::cout << "\nRESULT:FAIL\n"; return SDK_FAILURE;}
  std::cout << "\nRESULT:PASS\n";
  test.cleanup();
  test.printStats();

  return SDK_SUCCESS;
}
