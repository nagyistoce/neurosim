
/*Compiler options:

-I dir
-D name=definition
-g
-O0
-fbin-source -fbin-llvmir -fbin-amdil -fbin-exe
-fno-bin-source -fno-bin-llvmir -fno-bin-amdil -fno-bin-exe

export AMD_OCL_BUILD_OPTIONS_APPEND="-g -O0"
export AMD_OCL_BUILD_OPTIONS="-g -O0"
*/

#ifndef INC_DEFINITIONS_H
#define INC_DEFINITIONS_H

/*Somer parameters here should be passed rather then hard-coded:*/

/*Connectivity density*/
#ifndef P_CONNECT
  #define P_CONNECT               2
#endif

/*Connectivity ratio (% of excitatory): */
#ifndef R_CONNECT
  #define R_CONNECT               0.94
#endif

/*IZ model variable tolerance select (0-15): */
#define IZ_TOL                    7

/*If ZERO_TOL ==  1 it overwrites IZ_TOL with zero tolerance*/
#ifndef ZERO_TOL
  #define ZERO_TOL                1
#endif

/*Newton-Raphson error tolerance: */
#ifndef NR_TOL
  #define NR_TOL                  (1e-8f)
#endif

/*Verification error tolerance: */
#define VERIFY_TOL                (1e-7f)

/*If VERIFY_ZERO_TOL ==  1 it overwrites VERIFY_TOL with zero tolerance*/
#ifndef VERIFY_ZERO_TOL
  #define VERIFY_ZERO_TOL         1 
#endif

/*Math macros:*/
#define TESTPOW2(x)           (x && !( (x-1) & x ))
#define FLOAT_UNFLIP(f)       (f ^ (((f >> 31) - 1) | 0x80000000))
#define FLOAT_FLIP(f)         (f ^ (-int(f >> 31) | 0x80000000))
#define DIV_POW2(i, j)        (i >> j)
#define INT_CEIL(i,j)         (i/j + ((i%j) > 0))
#define ADD(i, j)             (i+j)
#define MUL(i, j)             (i*j)
#define MAD(i, j, k)          (i*j+k)
#define DIV(i, j)             (i/j)
#define MUL_UINT24(i, j)      (i*j)

/*Binary conversion*/
#define BIN_0000 0
#define BIN_0001 1
#define BIN_0010 2
#define BIN_0011 3
#define BIN_0100 4
#define BIN_0101 5
#define BIN_0110 6
#define BIN_0111 7
#define BIN_1000 8
#define BIN_1001 9
#define BIN_1010 a
#define BIN_1011 b
#define BIN_1100 c
#define BIN_1101 d
#define BIN_1110 e
#define BIN_1111 f

#define BIN_8_HEXIFY(b1,b2) (0x ## b1 ## b2)
#define BIN_8_RELAY(b1,b2) BIN_8_HEXIFY(b1, b2)
#define BIN_8(b1,b2) BIN_8_RELAY(BIN_ ## b1, BIN_ ## b2)

#define BIN_16_HEXIFY(b1,b2,b3,b4) (0x ## b1 ## b2 ## b3 ## b4)
#define BIN_16_RELAY(b1,b2,b3,b4) BIN_16_HEXIFY(b1, b2, b3, b4)
#define BIN_16(b1,b2,b3,b4) BIN_16_RELAY(BIN_##b1, BIN_##b2, BIN_##b3, BIN_##b4)

/* Print macros: */
#define PRINTF              printf
#define PRINTERR(...)       fprintf(stderr,__VA_ARGS__)
#define PRINT_HEX_TO_FILE(filePtr, byteSize, data)\
  {\
    unsigned int i;\
    for(i = (sizeof(byteSize)-1)*8; i > 0; i-=8)\
    {\
      fprintf (filePtr,"%02X ",(unsigned char)(data >> i));\
    }\
    fprintf (filePtr,"%02X",(unsigned char)(data >> i));\
  }
#define PRINT_HEX(byteSize, data)\
  {\
    long dataLong = *((long *)(&(data)));\
    unsigned int i;\
    for(i = (sizeof(byteSize)-1)*8; i > 0; i-=8)\
    {\
      PRINTF("%02X ",(unsigned char)(dataLong >> i));\
    }\
    PRINTF("%02X",(unsigned char)(dataLong >> i));\
  }

/*Checksum*/
/*For finding primes refer to http://www.archimedes-lab.org/primOmatic.html*/
#define CHECKSUM01(checksum, data)\
  {\
    checksum += ((1705662821 + data)%2147483659);\
  }
  
/*Memory size and name registration for statistics*/
#define REGISTER_MEMORY(kernel_name, mem_type, mem_name)\
  {\
    switch (mem_type)\
    {\
      case MEM_CONSTANT:\
      {\
        map<std::string, cl_uint> memSizes = memStats.cmSizes[kernel_name];\
        memSizes[#mem_name] = mem_name ##SizeBytes;\
        memStats.cmSizes[kernel_name] = memSizes;\
      }\
        break;\
      case MEM_GLOBAL:\
      {\
        map<std::string, cl_uint> memSizes = memStats.gmSizes[kernel_name];\
        memSizes[#mem_name] = mem_name ##SizeBytes;\
        memStats.gmSizes[kernel_name] = memSizes;\
      }\
        break;\
      case MEM_LOCAL:\
      {\
        map<std::string, cl_uint> memSizes = memStats.lmSizes[kernel_name];\
        memSizes[#mem_name] = mem_name ##SizeBytes;\
        memStats.lmSizes[kernel_name] = memSizes;\
      }\
        break;\
      default:\
        std::cout << "REGISTER_MEMORY: unsupported memory type\n" << std::endl;\
        return SDK_FAILURE;\
    }\
    set<std::string> kernels = memStats.kernelNames;\
    kernels.insert(kernel_name);\
    memStats.kernelNames = kernels;\
  }
/*Memory types*/
#define MEM_CONSTANT                                          0
#define MEM_GLOBAL                                            1
#define MEM_LOCAL                                             2

/*Calloc host memory*/
#define CALLOC(name, type, size)\
  {\
    name ##Size = size;\
    name ##SizeBytes = (name ##Size) * sizeof(type);\
    name = (type *)calloc(name ##Size, sizeof(type));\
    if(name == NULL)\
    {\
      sampleCommon->error("Failed to allocate memory for (" #name ")");\
      return SDK_FAILURE;\
    }\
  }
  
/*Reference host memory to another memory: name1 is referenced to name2 assuming that both have 
  the same charackteristics*/
#define REFERENCE(name1, name2)\
  {\
    name1 ##Size = name2 ##Size;\
    name1 ##SizeBytes = name2 ##SizeBytes;\
    name1 = name2;\
    if(name1 == NULL)\
    {\
      sampleCommon->error("Failed to allocate memory for (" #name1 ")");\
      return SDK_FAILURE;\
    }\
  }
  
/*Allocate buffer*/
#define CREATE_BUFFER(flags, buffer, size)\
  {\
    cl_int err = CL_SUCCESS;\
    if(isAmdPlatform())\
    {\
      buffer = cl::Buffer(context, flags, size, NULL, &err);\
      if(!sampleCommon->checkVal(err, CL_SUCCESS, "Buffer::Buffer() failed. (" #buffer ")"))\
      {\
        return SDK_FAILURE;\
      }\
    }\
    else\
    {\
      std::cout << "Unable to allocate buffer " << #buffer << " on non-AMD platform\n";\
      return SDK_FAILURE;\
    }\
  }
  
/*Enqueu write buffer*/
#define ENQUEUE_WRITE_BUFFER(block, buffer, size, data)\
  {\
    cl_int status;\
    cl::Event writeEvt = NULL;\
    cl_int eventStatus = CL_QUEUED;\
    status = commandQueue.enqueueWriteBuffer(buffer, block, 0, size, data, NULL, &writeEvt);\
    if(!sampleCommon->checkVal(status, CL_SUCCESS,\
      "CommandQueue::enqueueWriteBuffer() failed. (" #buffer ")"))\
    {\
      return SDK_FAILURE;\
    }\
    status = commandQueue.flush();\
    if(!sampleCommon->checkVal(status, CL_SUCCESS,\
      "cl::CommandQueue.flush failed. (" #buffer ")"))\
    {\
      return SDK_FAILURE;\
    }\
    eventStatus = CL_QUEUED;\
    while(eventStatus != CL_COMPLETE)\
    {\
      status = writeEvt.getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS,&eventStatus);\
      if(!sampleCommon->checkVal(status, CL_SUCCESS,\
        "enqueueWriteBuffer, cl:Event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS) failed. (" #buffer ")"))\
      {\
        return SDK_FAILURE;\
      }\
    }\
  }
  
/*Enqueue read buffer*/
#define ENQUEUE_READ_BUFFER(block, buffer, size, data)\
  {\
    cl_int status;\
    cl_int eventStatus = CL_QUEUED;\
    cl::Event readEvt;\
    status = commandQueue.enqueueReadBuffer(buffer, block, 0, size, data, NULL, &readEvt);\
    if(!sampleCommon->checkVal(status, CL_SUCCESS,\
      "CommandQueue::enqueueReadBuffer failed. (" #buffer ")"))\
    {\
      return SDK_FAILURE;\
    }\
    status = commandQueue.flush();\
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "cl::CommandQueue.flush failed."))\
    {\
      return SDK_FAILURE;\
    }\
    eventStatus = CL_QUEUED;\
    while(eventStatus != CL_COMPLETE)\
    {\
      status = readEvt.getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS, &eventStatus);\
      if(!sampleCommon->checkVal(status, CL_SUCCESS, \
        "enqueueReadBuffer, cl:Event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS) failed. (" #buffer ")"))\
      {\
        return SDK_FAILURE;\
      }\
    }\
  }
  
/*Set kernel argument*/
#define SET_KERNEL_ARG(kernel, arg, argNum)\
  {\
    cl_int status;\
    status = kernel.setArg(argNum, arg);\
    if(!sampleCommon->checkVal(status, CL_SUCCESS, "Kernel::setArg() failed. (" #arg ")"))\
    {\
      return SDK_FAILURE;\
    }\
  }
  
/*Target Platform Vendor*/
#define TARGET_PLATFORM_VENDOR                                "Advanced Micro Devices, Inc."

/*Target Device*/
#define TARGET_DEVICE_NAME                                    "Cayman"
//#define TARGET_DEVICE_NAME                                    "AMD Engineering Sample"

/*String used to tag statistics relevant to all kernels*/
#define KERNEL_ALL                                            "All Kernels"

/*Enable error try and catch*/
#define ERROR_TRY_CATCH_ENABLE                                0
/*Functional verification of each kernel*/
#define KERNEL_VERIFY_ENABLE                                  1
/*Enable high-level verification of sort results*/
#define SORT_VERIFY_ENABLE                                    0
#define EVENT_TIME_SLOTS                                      16
#define TOTAL_NEURON_BITS                                     16
#define WF_SIZE_WI                                            64
#define DATA_TYPE                                             float

/*Simulation parameters*/
#define SIMULATION_STEP_SIZE                                  1
#define SIMULATION_TIME_STEPS                                 5*EVENT_TIME_SLOTS

/*Reference simulation parameters*/
/*Synaptic event queue size limit per nrn: */
#define REFERENCE_EVENT_QUEUE_SIZE                            1000
#define FLIXIBLE_DELAYS_ENABLE                                1
  
/*
  Enable simulation stages
  Options: 
            - any single bit is set: unit test for that stage
            - some contiguios set of bits starting from MSB: either unit tests or integration tests 
              for enabled stages depending on the stage combination
            - all bits are set: complete integration test
*/
#define ENABLE_MASK                                           BIN_16(1111,1111,1000,0000)
//#define ENABLE_MASK                                           BIN_16(0000,0001,0000,0000)
//#define ENABLE_MASK                                           BIN_16(0000,0001,1000,0000)
#define EXPAND_EVENTS_ENABLE                                  (ENABLE_MASK&32768)
#define SCAN_ENABLE_V00                                       (ENABLE_MASK&16384)
#define GROUP_EVENTS_ENABLE_V00                               (ENABLE_MASK&8192)
#define SCAN_ENABLE_V01                                       (ENABLE_MASK&4096)
#define GROUP_EVENTS_ENABLE_V01                               (ENABLE_MASK&2048)
#define GROUP_EVENTS_ENABLE_V02                               (ENABLE_MASK&1024)
#define GROUP_EVENTS_ENABLE_V03                               (ENABLE_MASK&512)
#define MAKE_EVENT_PTRS_ENABLE                                (ENABLE_MASK&264)
#define UPDATE_NEURONS_ENABLE_V00                             (ENABLE_MASK&128)


/*
                                      General Restrictions
*/




/*
  Configuration interface for each kernel exists for tuning kernels to desired functionality and
  performance. Preprocessor allows to define parameters as literals, which reduces memory use, 
  logic and potentially increases performance.
  Configration interfaces may have several variants.
  Each variant has specific to it configuration interface parameters in addition to generic ones.
  On the host side the specific paremeters always defined by PARAMETER_NAME_VXX directives.
  On the device side the specific paremeters are enabled by VARIANT_NAME_DEVICE_VXX directives, 
  which can be passed with kernel compile flags.
  Host and device directive names are different because host code may be using several variants
  whereas the device code (kernel) is always compiled as a single variant.
  Both types of directives for the same variant have to have the same values.
*/
/***************************************************************************************************

  Configuration interface with propagation kernel for synaptic events
  
***************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
  /*Verification, Testing*/
  #define EXPAND_EVENTS_VERIFY_ENABLE                        KERNEL_VERIFY_ENABLE
  #define EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT              17
  /*Debugging*/
  #define EXPAND_EVENTS_DEBUG_ENABLE                         0
  #define EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS              (EXPAND_EVENTS_TIME_SLOTS*EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*EXPAND_EVENTS_GRID_SIZE_WG)//(1024*1024)
  /*Error tracking and codes*/
  #define EXPAND_EVENTS_ERROR_TRY_CATCH_ENABLE               ERROR_TRY_CATCH_ENABLE /*Options: 0, 1*/
  #define EXPAND_EVENTS_ERROR_TRACK_ENABLE                   1
  #define EXPAND_EVENTS_ERROR_BUFFER_SIZE_WORDS              1
  #define EXPAND_EVENTS_ERROR_CODE_1                         0x1 /*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE overflow*/
  /*Kernel parameters*/
  #define EXPAND_EVENTS_KERNEL_FILE_NAME                     "Kernel_ExpandEvents.cl"
  #define EXPAND_EVENTS_KERNEL_NAME                          "expand_events"
  #define EXPAND_EVENTS_WF_SIZE_WI                           WF_SIZE_WI /*Options: 64*/
  #define EXPAND_EVENTS_WG_SIZE_WF                           4
  #define EXPAND_EVENTS_WG_SIZE_WI                           (EXPAND_EVENTS_WG_SIZE_WF*EXPAND_EVENTS_WF_SIZE_WI)
  #define EXPAND_EVENTS_GRID_SIZE_WG                         (EXPAND_EVENTS_SPIKE_PACKETS/EXPAND_EVENTS_SPIKE_PACKETS_PER_WG)
  /*Spike data structure parameters (data in)*/
  #define EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE               63
  #define EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE             2
  #define EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS           2
  #define EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS              (EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + \
                                                            EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE * EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS)
  #define EXPAND_EVENTS_SPIKE_PACKETS                        128
  #define EXPAND_EVENTS_SPIKE_PACKETS_PER_WG                 1 /*For cases > 1 "unrolled loop" version can be implemented*/
  /*Synaptic data structure parameters*/
  #define EXPAND_EVENTS_MAX_SYNAPTIC_DATA_SIZE               512 /*1000, max number of synapses per neuron*/
  #define EXPAND_EVENTS_SYNAPTIC_DATA_UNIT_SIZE_WORDS        3 /*Size of a single synapse, 32-bit words*/
  /*Histogram of target neurons*/
  #define EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM              1
  #define EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT                  0 /*Options, LSB: 0; MSB:(EXPAND_EVENTS_TOTAL_NEURON_BITS-EXPAND_EVENTS_HISTOGRAM_BIN_BITS)*/
  #define EXPAND_EVENTS_HISTOGRAM_BIN_BITS                   4 /*Options: 4-8. Must be aligned with EXPAND_EVENTS_HISTOGRAM_BIN_MASK*/
  #define EXPAND_EVENTS_HISTOGRAM_BIN_MASK                   0xF /*Must be aligned with EXPAND_EVENTS_HISTOGRAM_BIN_BITS*/
  #define EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS                 (1<<EXPAND_EVENTS_HISTOGRAM_BIN_BITS)
  /*Time slot (bin) buffer parameters (data out)*/
  #define EXPAND_EVENTS_TIME_SLOTS                           EVENT_TIME_SLOTS /*Options: 16*/
  #define EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS               EXPAND_EVENTS_GRID_SIZE_WG
  #define EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS           3
  #define EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE  (21*1024)
  /*NN parameters*/
  #define EXPAND_EVENTS_TOTAL_NEURON_BITS                    TOTAL_NEURON_BITS /*16*/
  #define EXPAND_EVENTS_TOTAL_NEURONS                        (1<<EXPAND_EVENTS_TOTAL_NEURON_BITS)
  #define EXPAND_EVENTS_MAX_DELAY                           (EXPAND_EVENTS_TIME_SLOTS-SIMULATION_STEP_SIZE)
  #define EXPAND_EVENTS_MIN_DELAY                           1.0f
#endif



/***************************************************************************************************

  Configuration interface with scan kernel for histogram of synaptic events
  
***************************************************************************************************/
#if SCAN_ENABLE_V00 || SCAN_ENABLE_V01
/*
                                               Generic parameters
*/
  /*Verification*/
  #define SCAN_VERIFY_ENABLE                                  KERNEL_VERIFY_ENABLE
  /*Debugging*/
  #define SCAN_DEBUG_ENABLE                                   0
  #define SCAN_DEBUG_BUFFER_SIZE_WORDS                        (SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE)
  /*Error tracking and codes*/
  #define SCAN_ERROR_TRACK_ENABLE                             0
  #define SCAN_ERROR_BUFFER_SIZE_WORDS                        1
  #define SCAN_ERROR_CODE_1                                   0x1
  /*Scan kernel size parameters*/                               
  #define SCAN_KERNEL_FILE_NAME                               "Kernel_ScanHistogram.cl"
  #define SCAN_KERNEL_NAME                                    "scan_histogram"
  #define SCAN_WF_SIZE_WI                                     WF_SIZE_WI /*Options: 64*/
  #define SCAN_WG_SIZE_WF                                     2 /*Valid: 2, 4*/
  #define SCAN_WG_SIZE_WI                                     (SCAN_WG_SIZE_WF*SCAN_WF_SIZE_WI)
  #define SCAN_GRID_SIZE_WG                                   1
  /*Histogram of target neurons*/
  #define SCAN_HISTOGRAM_TIME_SLOTS                           EVENT_TIME_SLOTS /*Options: 16*/
  #define SCAN_HISTOGRAM_ELEMENTS_PER_WI                      16 /*Must div by 4. Cypress: 16, Cayman: 20*/
  #define SCAN_USE_2LEVEL_REDUCE                              1
  /*Maximum count of events per histogram item*/
  #define SCAN_HISTOGRAM_MAX_COUNT                            (0xFFFF)
/*
                                               Variants
*/
/*VARIANT_00*/
#if SCAN_ENABLE_V00
  /*Device*/
#ifdef SCAN_DEVICE_V00
  #define SCAN_HISTOGRAM_TOTAL_BINS                           (1<<4)
  #define SCAN_HISTOGRAM_BIN_SIZE                             128
  #define SCAN_HISTOGRAM_IN_TYPE_V00                          0
#endif
  /*Host and Device*/
  #define SCAN_HISTOGRAM_TOTAL_BINS_V00                       (1<<4)
  #define SCAN_HISTOGRAM_BIN_SIZE_V00                         128
  #define SCAN_HISTOGRAM_IN_TYPE_V00                          0
#endif  
/*END VARIANT_00*/
/*VARIANT_01*/
#if SCAN_ENABLE_V01
  /*Device*/
#ifdef SCAN_DEVICE_V01
  #define SCAN_HISTOGRAM_TOTAL_BINS                           (1<<4)
  #define SCAN_HISTOGRAM_BIN_SIZE                             128
  #define SCAN_HISTOGRAM_IN_TYPE                              1
#endif
  /*Host and Device*/
  #define SCAN_HISTOGRAM_TOTAL_BINS_V01                       (1<<4)
  #define SCAN_HISTOGRAM_BIN_SIZE_V01                         128
  #define SCAN_HISTOGRAM_IN_TYPE_V01                          1
  #define SCAN_HISTOGRAM_BIN_BACKETS                          128
#endif  
/*END VARIANT_01*/
/*
                                               Restrictions
*/
#if	(SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE > SCAN_WG_SIZE_WI*SCAN_HISTOGRAM_ELEMENTS_PER_WI)
  #error (Total number of histogram elements is more than kernel can handle)
#endif
#if	(SCAN_HISTOGRAM_ELEMENTS_PER_WI%4 != 0)
  #error (Parameter SCAN_HISTOGRAM_ELEMENTS_PER_WI must be divisible by 4)
#endif
#if(SCAN_WG_SIZE_WF != 2 && SCAN_WG_SIZE_WF != 4)
  #error (Parameter SCAN_WG_SIZE_WF must be either 2 or 4)
#endif
#if(!SCAN_USE_2LEVEL_REDUCE)
  #error (Set option SCAN_USE_2LEVEL_REDUCE to 1)
#endif
#endif



/***************************************************************************************************

  Configuration interface with kernel for grouping synaptic events
  
***************************************************************************************************/
#if GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03
/*
                                               Generic parameters
*/
  /*Verification*/
  #define GROUP_EVENTS_VERIFY_ENABLE                            KERNEL_VERIFY_ENABLE
  /*Debugging*/
  #define GROUP_EVENTS_DEBUG_ENABLE                             0
  #define GROUP_EVENTS_DEBUG_BUFFER_SIZE_WORDS                  (GROUP_EVENTS_TIME_SLOTS*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT*GROUP_EVENTS_HISTOGRAM_BIN_SIZE_OUT)//(1024*128)
  /*Error tracking and codes*/
  #define GROUP_EVENTS_ERROR_TRY_CATCH_ENABLE                   0
  #define GROUP_EVENTS_ERROR_TRACK_ENABLE                       0
  #define GROUP_EVENTS_ERROR_BUFFER_SIZE_WORDS                  1
  #define GROUP_EVENTS_ERROR_CODE_1                             0x1
  /*Kernel size parameters*/
  #define GROUP_EVENTS_KERNEL_FILE_NAME                         "Kernel_GroupEvents.cl"
  #define GROUP_EVENTS_KERNEL_NAME                              "group_events"
  #define GROUP_EVENTS_WF_SIZE_WI                               WF_SIZE_WI /*Options: 64*/
  #define GROUP_EVENTS_WG_SIZE_WF                               1
  #define GROUP_EVENTS_WG_SIZE_WI                               (GROUP_EVENTS_WG_SIZE_WF*GROUP_EVENTS_WF_SIZE_WI)
  #define GROUP_EVENTS_GRID_SIZE_WG                             (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS/GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG)
  /*Incoming histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_BITS                       4   /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_MASK*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_MASK                       0xF /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_BITS*/
  #define GROUP_EVENTS_HISTOGRAM_TOTAL_BINS                     (1<<GROUP_EVENTS_HISTOGRAM_BIN_BITS)
  #define GROUP_EVENTS_HISTOGRAM_BIN_SIZE                       (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS)
  #define GROUP_EVENTS_HISTOGRAM_MAX_COUNT                      (0xFFFF)/*(GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE)*/
  /*Outgoing histogram of target neurons*/
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT              1
  #define GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT                   4   /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT                   0xF /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT*/
  #define GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT                 (1<<GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT)
  #define GROUP_EVENTS_HISTOGRAM_BIN_SIZE_OUT                   (GROUP_EVENTS_GRID_SIZE_WG)
  #define GROUP_EVENTS_HISTOGRAM_OUT_ELEMENTS_PER_WG            (GROUP_EVENTS_ELEMENTS_PER_WI*GROUP_EVENTS_WG_SIZE_WI)
  #define GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE                  (GROUP_EVENTS_GRID_SIZE_WG)
  /*Synaptic event buffer parameters*/
  #define GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS                   128
  #define GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS               3
  #define GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE           (21*1024)
  #define GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE           (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS*GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE*10/10) /*Max: GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE x GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS*/
  #define GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG            1
  /*NN parameters*/
  #define GROUP_EVENTS_TOTAL_NEURON_BITS                        TOTAL_NEURON_BITS /*16*/
  #define GROUP_EVENTS_TOTAL_NEURONS                            (1<<GROUP_EVENTS_TOTAL_NEURON_BITS)
  #define GROUP_EVENTS_MAX_DELAY                                (GROUP_EVENTS_TIME_SLOTS-SIMULATION_STEP_SIZE)
  #define GROUP_EVENTS_MIN_DELAY                                1.0f
  /*Time slots for storing events*/
  #define GROUP_EVENTS_TIME_SLOTS                               EVENT_TIME_SLOTS /*Options: 16*/
  /*Each WI takes this number of data elements*/
  #define GROUP_EVENTS_ELEMENTS_PER_WI                          4
  /*Check data boundary condition before fetching*/
  #define GROUP_EVENTS_CHECK_BOUNDARY                           1
  /*Reduction level option*/
  #define GROUP_EVENTS_USE_2LEVEL_REDUCE                        1
  /*Enable local sort*/
  #define GROUP_EVENTS_LOCAL_SORT_ENABLE                        1
  /*How many threads share a counter for computing local histogram*/
  #define GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER               16
/*
                                               Variants
*/
/*
            VARIANT_00
            Loads key data from original data structure. Computes values as pointers to the original
            data. Sorts by keys the key-value pairs based on provided scanned histogram. 
            Computes new histogram for next stage.
*/
#if GROUP_EVENTS_ENABLE_V00
  /*Device*/
#ifdef GROUP_EVENTS_DEVICE_V00
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      0
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  4
  /*Source events may be delivered in different data structures*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        0 /*Block-partitioned by exapnd kernel*/
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        0
  /*Each key gets a value of its address in GM*/
  #define GROUP_EVENTS_VALUES_MODE                              1
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        1
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          0
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              0
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
  /*Host and Device*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00                  0
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V00              4
  #define GROUP_EVENTS_VALUES_MODE_V00                          1
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET_V00                    1
#endif  /*END VARIANT_00*/
/*
            VARIANT_01
            Sorts by keys the key-value pairs based on provided scanned histogram. 
            Computes new histogram for next stage.
*/
#if GROUP_EVENTS_ENABLE_V01
  /*Device*/
#ifdef GROUP_EVENTS_DEVICE_V01
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      4
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  4
  /*Source events may be delivered in different data structures*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        1 /*Contiguous*/
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        1
  /*Each key gets a value carried out from before*/
  #define GROUP_EVENTS_VALUES_MODE                              2
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        0
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          0
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              0
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
  /*Host-visible copies of device parameters*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01                  4
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01              4
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V01                    1
  #define GROUP_EVENTS_VALUES_MODE_V01                          2
#endif/*END VARIANT_01*/
/*
            VARIANT_02
            This variant is the same as VARIANT_01 except it replaces data with new key at the end.
            This allows to continue sorting with new key in the next stages.
*/
#if GROUP_EVENTS_ENABLE_V02
  /*Device*/
#ifdef GROUP_EVENTS_DEVICE_V02
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      4
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  0
  /*Source events may be delivered in different data structures*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        1 /*Contiguous*/
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        1
  /*Each key gets a value carried out from before*/
  #define GROUP_EVENTS_VALUES_MODE                              2
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        0
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          0
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              1
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
  /*Host-visible copies of device parameters*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02                  4
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V02              0
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V02                    1
  #define GROUP_EVENTS_VALUES_MODE_V02                          2
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET_V02               0
#endif/*END VARIANT_02*/
/*
            VARIANT_03
            This variant is the same as VARIANT_01 except it relocate original values in place
            of pointers. No histogram for the next stage is computed.
*/
/*TODO: need to disable histogram out*/
#if GROUP_EVENTS_ENABLE_V03
  /*Device*/
#ifdef GROUP_EVENTS_DEVICE_V03
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      4
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  0
  /*Source events may be delivered in different data structures*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        1 /*Contiguous*/
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        1
  /*Each key gets a value carried out from before*/
  #define GROUP_EVENTS_VALUES_MODE                              2
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        0
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          1
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              0
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
  /*Host-visible copies of device parameters*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03                  4
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V03              0
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V03                    1
  #define GROUP_EVENTS_VALUES_MODE_V03                          2
#endif/*END VARIANT_03*/
#endif /*GROUP_EVENTS_ENABLE*/



/***************************************************************************************************

  Kernel configuration interface for making synaptic pointer structure
  
***************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
/*
                                               Generic parameters
*/
  /*Verification*/
  /*CONTROL: enable verification*/
  #define MAKE_EVENT_PTRS_VERIFY_ENABLE                            KERNEL_VERIFY_ENABLE
  
  /*Debugging*/
  /*CONTROL: enable debug buffers*/
  #define MAKE_EVENT_PTRS_DEBUG_ENABLE                             0
  #define MAKE_EVENT_PTRS_DEBUG_BUFFER_SIZE_WORDS                  (128*MAKE_EVENT_PTRS_WG_SIZE_WI*MAKE_EVENT_PTRS_GRID_SIZE_WG)
  
  /*Error tracking*/
  /*CONTROL: enable error logging at kernel level*/
  #define MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE                       1
  #define MAKE_EVENT_PTRS_ERROR_BUFFER_SIZE_WORDS                  1
  #define MAKE_EVENT_PTRS_ERROR_CODE_1                             0x1
  
  /*Kernel parameters*/
  #define MAKE_EVENT_PTRS_KERNEL_FILE_NAME                          "Kernel_MakeEventPointers.cl"
  #define MAKE_EVENT_PTRS_KERNEL_NAME                               "make_event_pointers"
  #define MAKE_EVENT_PTRS_WF_SIZE_WI                                WF_SIZE_WI
  /*CONTROL: number of WFs in a WG*/
  #define MAKE_EVENT_PTRS_WG_SIZE_WF                                2
  #define MAKE_EVENT_PTRS_WG_SIZE_WI                                (MAKE_EVENT_PTRS_WG_SIZE_WF*MAKE_EVENT_PTRS_WF_SIZE_WI)
  /*CONTROL: each WI processes this number of data elements*/
  #define MAKE_EVENT_PTRS_ELEMENTS_PER_WI                          4
  
  /*NN parameters*/
  /*CONTROL: total neurons in the network*/
  #define MAKE_EVENT_PTRS_TOTAL_NEURON_BITS                         TOTAL_NEURON_BITS /*16*/
  #define MAKE_EVENT_PTRS_TOTAL_NEURONS                             (1<<(MAKE_EVENT_PTRS_TOTAL_NEURON_BITS))
  
  /*Buffer and memory sizes*/
  #define MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE                      2
  /*Synaptic event input buffer parameters*/
  #define MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS                   3
  /*LM allocation for scan per WF. 
  TODO: this can be reduced by 18 at least*/
  #define MAKE_EVENT_PTRS_SCAN_WF_LM_SHARE                         192
  /*Offset in GM for location of total events*/
  #define MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET                      (128*16)
  
  /*Event delivery mode:
    0 - array size of neuron count
    1 - structs*/
  #define MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE                       0
  
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  /*CONTROL: WG count*/
  #define MAKE_EVENT_PTRS_GRID_SIZE_WG                              64
  
  /*Output buffer size*/
  #define MAKE_EVENT_PTRS_STRUCTS                                   1
  #define MAKE_EVENT_PTRS_STRUCT_SIZE                               MAKE_EVENT_PTRS_TOTAL_NEURONS*\
                                                                    MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE+ \
                                                                    (2*MAKE_EVENT_PTRS_GRID_SIZE_WG*\
                                                                    MAKE_EVENT_PTRS_WG_SIZE_WF*\
                                                                    MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE)
  /*Unit test related parameters*/
  /*CONTROL: Source data size for unit test. Can be viewed as elements per neuron*/
  #define MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE                 (8*MAKE_EVENT_PTRS_TOTAL_NEURONS)
  /*CONTROL: Max neuron ID for unit test.*/
  #define MAKE_EVENT_PTRS_TEST_MAX_NEURON_ID                       (MAKE_EVENT_PTRS_TOTAL_NEURONS-1)
  
  /*CONTROL: LM allocation for computing event counts per WF*/
  #define MAKE_EVENT_PTRS_EVENT_COUNT_WF_LM_SHARE                   512
  
  /*Small post-process kernel to glue WF results together*/
  #define GLUE_EVENT_PTRS_KERNEL_FILE_NAME                          "Kernel_MakeEventPointers.cl"
  #define GLUE_EVENT_PTRS_KERNEL_NAME                               "glue_event_pointers"
  #define GLUE_EVENT_PTRS_WF_SIZE_WI                                WF_SIZE_WI
  #define GLUE_EVENT_PTRS_WG_SIZE_WF                                1
  #define GLUE_EVENT_PTRS_WG_SIZE_WI                                (GLUE_EVENT_PTRS_WG_SIZE_WF*\
                                                                    GLUE_EVENT_PTRS_WF_SIZE_WI)                            
  #define GLUE_EVENT_PTRS_GRID_SIZE_WG                              1
  #define GLUE_EVENT_WF_LM_SHARE_SIZE                               (MAKE_EVENT_PTRS_GRID_SIZE_WG*\
                                                                    MAKE_EVENT_PTRS_WG_SIZE_WF)
  #define MAKE_EVENT_PTRS_GLUE_DEFAULT_VALUE                        0xFFFFFFFF

#elif MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1

  /*CONTROL: number of pointer structs*/
  #define MAKE_EVENT_PTRS_STRUCTS                                  (64*MAKE_EVENT_PTRS_WG_SIZE_WF)
  /*CONTROL: max size of pointer struct. Depends on max number of unique neurons in the source*/
  #define MAKE_EVENT_PTRS_STRUCT_SIZE                              (128*1024*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE)
  #define MAKE_EVENT_PTRS_GRID_SIZE_WG                             (MAKE_EVENT_PTRS_STRUCTS/MAKE_EVENT_PTRS_WG_SIZE_WF)
  /*CONTROL: Source data size for unit test. Can be viewed as KB of source per output struct*/
  #define MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE                 (128*1024*MAKE_EVENT_PTRS_STRUCTS)
  /*CONTROL: Max neuron ID for unit test.*/
  #define MAKE_EVENT_PTRS_TEST_MAX_NEURON_ID                       (MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE/\
                                                                    MAKE_EVENT_PTRS_WF_SIZE_WI)
#endif/*end MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE*/

/*
                                               Restrictions
*/
#if	((MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 1) && \
     (MAKE_EVENT_PTRS_STRUCTS%MAKE_EVENT_PTRS_WG_SIZE_WF != 0))
  #error (MAKE_EVENT_PTRS_STRUCTS has to be multiple of MAKE_EVENT_PTRS_WG_SIZE_WF)
#endif
#if	((MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0) && \
     (MAKE_EVENT_PTRS_STRUCTS != 1))
  #error (MAKE_EVENT_PTRS_STRUCTS has to be 1)
#endif
/*Local memory size*/
#define MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE             (MAKE_EVENT_PTRS_WF_SIZE_WI + 1)
#if (MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE < MAKE_EVENT_PTRS_SCAN_WF_LM_SHARE)
  #undef  MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE
  #define MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE           MAKE_EVENT_PTRS_SCAN_WF_LM_SHARE
#endif
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
#if (MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE < MAKE_EVENT_PTRS_EVENT_COUNT_WF_LM_SHARE)
  #undef  MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE
  #define MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE           MAKE_EVENT_PTRS_EVENT_COUNT_WF_LM_SHARE
#endif
#endif
#endif /*MAKE_EVENT_PTRS_ENABLE*/



/***************************************************************************************************

  Kernel configuration interface for updating model variables
  
***************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00 || UPDATE_NEURONS_ENABLE_V01
/*
                                               Generic parameters
*/
  /*Verification*/
  #define UPDATE_NEURONS_VERIFY_ENABLE                            KERNEL_VERIFY_ENABLE
  /*Debugging*/
  #define UPDATE_NEURONS_DEBUG_ENABLE                             0
  #define UPDATE_NEURONS_DEBUG_BUFFER_SIZE_WORDS                  (UPDATE_NEURONS_TOTAL_NEURONS*UPDATE_NEURONS_MODEL_VARIABLES)
  /*Error tracking and codes*/
  /*CONTROL: enable error logging at kernel level*/
  #define UPDATE_NEURONS_ERROR_TRACK_ENABLE                       0
  #define UPDATE_NEURONS_ERROR_BUFFER_SIZE_WORDS                  1
  #define UPDATE_NEURONS_ERROR_CODE_1                             0x1
  /*Kernel file and name*/
  #define UPDATE_NEURONS_KERNEL_FILE_NAME                         "Kernel_UpdateNeurons.cl"
  #define UPDATE_NEURONS_KERNEL_NAME                              "update_neurons"
  /*WF size measured in WIs*/
  #define UPDATE_NEURONS_WF_SIZE_WI                               WF_SIZE_WI /*Options: 64*/
  /*Scaling parameter for the grid size*/
  #define UPDATE_NEURONS_WF_DATA_CHUNKS                           8
  /*Size of element in event pointer struct*/
  #define UPDATE_NEURONS_STRUCT_ELEMENT_SIZE                      2
  /*Event input buffer element size*/
  #define UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS                   3
  /*CONTROL: Source data size for unit test.t*/
  #define UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE                 (8*UPDATE_NEURONS_TOTAL_NEURONS)
  /*CONTROL: each WI processes this number of neurons*/
  #define UPDATE_NEURONS_ELEMENTS_PER_WI                          1
  /*Size of time delay buffer*/
  #define UPDATE_NEURONS_TIME_SLOTS                               EVENT_TIME_SLOTS /*Options: 16*/
  /*Sizes of model variables and parameters*/
  #define UPDATE_NEURONS_MODEL_VARIABLES                          4
  #define UPDATE_NEURONS_MODEL_PARAMETERS                         9
  /*Spike data structure parameters*/
  #define UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE                   63
  #define UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE                 2
  #define UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS               2
  #define UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS                  (UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + \
                                                                  UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE * UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS)
  /*Number of constant coefficient types*/
  #define UPDATE_NEURONS_CONST_COEFFICIENTS                       4
  /*Accessor to constant coefficients*/
  #define CONST_CO(mem, type, index)                              *(mem + type*(UPDATE_NEURONS_PS_ORDER_LIMIT+1)+index)
  #define CONST_TOL(mem)                                          *(mem + UPDATE_NEURONS_CONST_COEFFICIENTS*(UPDATE_NEURONS_PS_ORDER_LIMIT+1))
  #define CONST_SIZE                                              (UPDATE_NEURONS_CONST_COEFFICIENTS*(UPDATE_NEURONS_PS_ORDER_LIMIT+1) + 1)
  
  /*NN parameters*/
  #define UPDATE_NEURONS_TOTAL_NEURON_BITS                        TOTAL_NEURON_BITS /*16*/
  #define UPDATE_NEURONS_TOTAL_NEURONS                            (1<<UPDATE_NEURONS_TOTAL_NEURON_BITS)
  #define UPDATE_NEURONS_MAX_DELAY                                (UPDATE_NEURONS_TIME_SLOTS-SIMULATION_STEP_SIZE)
  #define UPDATE_NEURONS_MIN_DELAY                                1.0f
  
  /*Simulation parameters*/
  #define UPDATE_NEURONS_ZERO_TOLERANCE_ENABLE                    1
  #define UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD       1
  /*Possible values: const double tols[16] = 
  {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,1e-16};*/
  #define UPDATE_NEURONS_PS_TOLERANCE                             (1.0e-7)
  #define UPDATE_NEURONS_NR_TOLERANCE                             (1.0e-7) /*normally no more than UPDATE_NEURONS_PS_TOLERANCE*/
  #define UPDATE_NEURONS_PS_ORDER_LIMIT                           64
  #define UPDATE_NEURONS_NR_ORDER_LIMIT                           20
  #define UPDATE_NEURONS_INJECT_CURRENT_STEPS                     SIMULATION_TIME_STEPS
  /*Possible values:  const double dt_vals[16] = {(double)1/4,(double)1/6,(double)1/8,(double)1/10,
                                                  (double)1/20,(double)1/40,(double)1/60,
                                                  (double)1/80,(double)1/100,(double)1/200,
                                                  (double)1/400,(double)1/600,(double)1/800,
                                                  (double)1/1000,(double)1/2000};*/
  #define UPDATE_NEURONS_DT                                       0.25f
  /*Selected neuron model parameters, constant for all neurons*/
  #define UPDATE_NEURONS_TAU_AMPA                                 5.0f
  #define UPDATE_NEURONS_TAU_GABA                                 10.0f
  #define UPDATE_NEURONS_C                                        200.0f
  #define UPDATE_NEURONS_a                                        0.03f
/*
                                               Restrictions
*/
#if UPDATE_NEURONS_ELEMENTS_PER_WI != 1
  #error (UPDATE_NEURONS_ELEMENTS_PER_WI has to be 1)
#endif
/*
                                               Variants
*/
/*
            VARIANT_00

*/
#if UPDATE_NEURONS_ENABLE_V00
  /*Device*/
#ifdef UPDATE_NEURONS_DEVICE_V00
  /*Deliver events by placing address to events for each neuron in dedicated for it location.
    Each neuron has such a location*/
  #define UPDATE_NEURONS_EVENT_DELIVERY_MODE                      0
  /*CONTROL: number of WFs in a WG*/
  #define UPDATE_NEURONS_WG_SIZE_WF                               2
  #define UPDATE_NEURONS_WG_SIZE_WI                               (UPDATE_NEURONS_WG_SIZE_WF*\
                                                                  UPDATE_NEURONS_WF_SIZE_WI)
  #define UPDATE_NEURONS_GRID_SIZE_WG                             (UPDATE_NEURONS_TOTAL_NEURONS/\
                                                                  (UPDATE_NEURONS_ELEMENTS_PER_WI*\
                                                                  UPDATE_NEURONS_WG_SIZE_WI*\
                                                                  UPDATE_NEURONS_WF_DATA_CHUNKS))
  #define UPDATE_NEURONS_STRUCT_SIZE                              (UPDATE_NEURONS_TOTAL_NEURONS*\
                                                                  UPDATE_NEURONS_STRUCT_ELEMENT_SIZE)
#endif
  /*Host and Device*/
  #define UPDATE_NEURONS_WG_SIZE_WF_V00                           2
  #define UPDATE_NEURONS_WG_SIZE_WI_V00                           (UPDATE_NEURONS_WG_SIZE_WF_V00*\
                                                                  UPDATE_NEURONS_WF_SIZE_WI)
  #define UPDATE_NEURONS_GRID_SIZE_WG_V00                         (UPDATE_NEURONS_TOTAL_NEURONS/\
                                                                  (UPDATE_NEURONS_ELEMENTS_PER_WI*\
                                                                  UPDATE_NEURONS_WG_SIZE_WI_V00*\
                                                                  UPDATE_NEURONS_WF_DATA_CHUNKS))
#if !MAKE_EVENT_PTRS_ENABLE
  #define UPDATE_NEURONS_STRUCT_SIZE_V00                          (UPDATE_NEURONS_TOTAL_NEURONS*\
                                                                  UPDATE_NEURONS_STRUCT_ELEMENT_SIZE)
#else
  #define UPDATE_NEURONS_STRUCT_SIZE_V00                          (UPDATE_NEURONS_TOTAL_NEURONS*\
                                                                  UPDATE_NEURONS_STRUCT_ELEMENT_SIZE+\
                                                                  (2*MAKE_EVENT_PTRS_GRID_SIZE_WG*\
                                                                  MAKE_EVENT_PTRS_WG_SIZE_WF*\
                                                                  MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE))
#endif
  #define UPDATE_NEURONS_STRUCTS_V00                              1
  #define UPDATE_NEURONS_SPIKE_PACKETS_V00                        (UPDATE_NEURONS_WG_SIZE_WF_V00*\
                                                                  UPDATE_NEURONS_GRID_SIZE_WG_V00)
//TODO: get rid of it by implementing limit address in the kernel
#if(UPDATE_NEURONS_TOTAL_NEURONS%(UPDATE_NEURONS_ELEMENTS_PER_WI*UPDATE_NEURONS_WG_SIZE_WI_V00*\
   UPDATE_NEURONS_WF_DATA_CHUNKS))
  #error(UPDATE_NEURONS_TOTAL_NEURONS must be multiple of (UPDATE_NEURONS_ELEMENTS_PER_WI*\
    UPDATE_NEURONS_WG_SIZE_WI_V00*UPDATE_NEURONS_WF_DATA_CHUNKS))
#endif
#endif  /*END VARIANT_00*/
/*
            VARIANT_01

*/
#if UPDATE_NEURONS_ENABLE_V01
  /*Device*/
#ifdef UPDATE_NEURONS_DEVICE_V01
  #define UPDATE_NEURONS_EVENT_DELIVERY_MODE                      1
  /*CONTROL: number of WFs in a WG*/
  #define UPDATE_NEURONS_WG_SIZE_WF                               1
  #define UPDATE_NEURONS_GRID_SIZE_WG                             (UPDATE_NEURONS_STRUCTS/UPDATE_NEURONS_WG_SIZE_WF)
  #define UPDATE_NEURONS_WG_SIZE_WI                               (UPDATE_NEURONS_WG_SIZE_WF*UPDATE_NEURONS_WF_SIZE_WI)
#endif
  /*Host and Device*/
  #define UPDATE_NEURONS_WG_SIZE_WF_V01                           1
  #define UPDATE_NEURONS_GRID_SIZE_WG_V01                         (UPDATE_NEURONS_STRUCTS/UPDATE_NEURONS_WG_SIZE_WF_V01)
  #define UPDATE_NEURONS_WG_SIZE_WI_V01                           (UPDATE_NEURONS_WG_SIZE_WF_V01*UPDATE_NEURONS_WF_SIZE_WI)
  /*CONTROL: number of pointer structs*/
  #define UPDATE_NEURONS_STRUCTS_V01                              (64*UPDATE_NEURONS_WG_SIZE_WF_V01)
  /*CONTROL: max size of pointer struct. Depends on max number of unique neurons receiving events*/
  #define UPDATE_NEURONS_STRUCT_SIZE_V01                          (1*1024*UPDATE_NEURONS_STRUCT_ELEMENT_SIZE)
  #undef UPDATE_NEURONS_SPIKE_PACKETS_V01
  #define UPDATE_NEURONS_SPIKE_PACKETS_V01                        (UPDATE_NEURONS_WG_SIZE_WF_V01*UPDATE_NEURONS_GRID_SIZE_WG_V01)
#endif  /*END VARIANT_00*/
#endif /*UPDATE_NEURONS_ENABLE*/

#endif /*INC_DEFINITIONS_H*/
