
/*

  Definitions and macros visible in all source files including kernels

*/

#ifndef INC_DEFINITIONS_H
#define INC_DEFINITIONS_H

/** ############################################################################################# **

  I. Generic Parameters, Macros and Definitions
  
** ############################################################################################# **/



/***************************************************************************************************
  Math macros
***************************************************************************************************/
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

/*Checksum*/
/*For finding primes refer to http://www.archimedes-lab.org/primOmatic.html*/
#define CHECKSUM01(checksum, data)\
  {\
    checksum += ((1705662821u + data) % 2147483659u);\
  }
  
#define SET_RANDOM_SEED(seed)\
  {\
    time_t t = time(NULL);\
    seed = *((unsigned int *)(&t)) + srandCounter;\
    srandCounter++;\
    srand(seed);\
  }
  
#define GET_RANDOM_INT(setValue, max, minPercent, maxPercent)\
  {\
    if(minPercent > maxPercent){setValue = -1;}\
    if((minPercent > 100.0) || (maxPercent > 100.0)){setValue = -1;}\
    if((minPercent < 0.0) || (maxPercent < 0.0)){setValue = -1;}\
    if(setValue != -1)\
    {\
      if(minPercent == maxPercent){setValue = cl_uint(((double)max)*(maxPercent/100.0));}\
      else\
      {\
        setValue = cl_uint(((double)(max))*((minPercent/100.0) + \
          abs(((maxPercent-minPercent)/100.0)*((double)rand()/((double)RAND_MAX)))));\
      }\
    }\
  }
/**************************************************************************************************/



/***************************************************************************************************
  Print macros
***************************************************************************************************/
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
/**************************************************************************************************/



/***************************************************************************************************
  Memory allocation and registration
***************************************************************************************************/
/*Memory size and name registration for statistics*/
#define REGISTER_TIME(kernel_name, time, increment)\
  {\
    map<std::string, double> stats = kernelStats.execTime[#kernel_name];\
    if(stats.find("Time") ==  stats.end()){stats["Time"] = 0;}\
    if(stats.find("Count") ==  stats.end()){stats["Count"] = 0;}\
    stats["Time"] += time;\
    stats["Count"] += increment;\
    kernelStats.execTime[#kernel_name] = stats;\
    set<std::string> kernels = kernelStats.kernelNamesExecTime;\
    kernels.insert(#kernel_name);\
    kernelStats.kernelNamesExecTime = kernels;\
  }
  
/*Memory size and name registration for statistics*/
#define REGISTER_MEMORY(kernel_name, mem_type, mem_name)\
  {\
    switch (mem_type)\
    {\
      case MEM_CONSTANT:\
      {\
        map<std::string, cl_uint> memSizes = kernelStats.cmSizes[kernel_name];\
        memSizes[#mem_name] = mem_name ##SizeBytes;\
        kernelStats.cmSizes[kernel_name] = memSizes;\
      }\
        break;\
      case MEM_GLOBAL:\
      {\
        map<std::string, cl_uint> memSizes = kernelStats.gmSizes[kernel_name];\
        memSizes[#mem_name] = mem_name ##SizeBytes;\
        kernelStats.gmSizes[kernel_name] = memSizes;\
      }\
        break;\
      case MEM_LOCAL:\
      {\
        map<std::string, cl_uint> memSizes = kernelStats.lmSizes[kernel_name];\
        memSizes[#mem_name] = mem_name ##SizeBytes;\
        kernelStats.lmSizes[kernel_name] = memSizes;\
      }\
        break;\
      default:\
        std::cout << "REGISTER_MEMORY: unsupported memory type\n" << std::endl;\
        return SDK_FAILURE;\
    }\
    set<std::string> kernels = kernelStats.kernelNames;\
    kernels.insert(kernel_name);\
    kernelStats.kernelNames = kernels;\
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
/**************************************************************************************************/



/***************************************************************************************************
  OCL Macros and Constants
***************************************************************************************************/
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
        "enqueueWriteBuffer, \
        cl:Event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS) failed. (" #buffer ")"))\
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
        "enqueueReadBuffer, \
        cl:Event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS) failed. (" #buffer ")"))\
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
  
/*Define OpenCL 1.2 consants*/
#if !defined (CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT)
  #define CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT                 (1 << 7)
#endif
/**************************************************************************************************/


  
/***************************************************************************************************
  Log and Statistics Macros and Definitions
***************************************************************************************************/
/*Log report messages*/
#if !(defined(LOG_REPORT))
  #define LOG_REPORT                                          0
#endif
#if !(defined(LOG_REPORT_FILE_NAME))
  #define LOG_REPORT_FILE_NAME                                "log_report.txt"
#endif
#define LOG_MODEL_VARIABLES_NEURON_ID                         1
#define LOG_SIMULATION_FILE_NAME                              "log_simulation.txt"
#define LOG_MODEL_VARIABLES_FILE_NAME                         "log_model_variable_trace.csv"
#define LOG_MODEL_VARIABLES_FILE_HEADER                       "v,u,g_ampa,g_gaba,v_peak,I"
#define LOG_MODEL_VARIABLES_FILE_BODY(i)                      << nrn_ps[i].v << "," \
                                                              << nrn_ps[i].u << "," \
                                                              << nrn_ps[i].g_ampa << "," \
                                                              << nrn_ps[i].g_gaba << "," \
                                                              << nrn_ps[i].v_peak << "," \
                                                              << nrn_ps[i].I << std::endl
/*Compile methods for taking a snapshot of simulation state*/
#if !(defined(SIMULATION_SNAPSHOT))
  #define SIMULATION_SNAPSHOT                                 0
#endif
#if SIMULATION_SNAPSHOT
  #if !(defined(LOG_SNAPSHOT_FILE_NAME))
    #define LOG_SNAPSHOT_FILE_NAME                            "log_snapshot.txt"
  #endif
#endif
/**************************************************************************************************/



/***************************************************************************************************
  Device and platform parameters
***************************************************************************************************/
/*Wave size*/
#define WF_SIZE_WI                                            64
/*Target devices in the order of preference*/
#if !(defined(TARGET_DEVICE_NAME))
  #define TARGET_DEVICE_NAME                                  "Tahiti,Cayman"
#endif
/*Target Platform Vendor
#define TARGET_PLATFORM_VENDOR                                "Advanced Micro Devices, Inc."*/
/**************************************************************************************************/



/***************************************************************************************************
  Simulation parameters
***************************************************************************************************/
/*Network size in terms of bits (1<<TOTAL_NEURON_BITS)*/
#if !(defined(TOTAL_NEURON_BITS))
  #define TOTAL_NEURON_BITS                                   17 /*17*/
#endif
/*Max number of synapses per neuron*/
#if !(defined(MAX_SYNAPSES_PER_NEURON))
  #define MAX_SYNAPSES_PER_NEURON                             (30*64) /*(20*64)*/
#endif
/*Synapse deviation ratio from max defined by MAX_SYNAPSES_PER_NEURON*/
#define SYNAPSE_DEVIATION_RATIO                               0.5
/*Size of spike packet buffer (spikes) per WF*/
#if !(defined(SPIKE_DATA_BUFFER_SIZE))
  #define SPIKE_DATA_BUFFER_SIZE                              127 /*127*/
#endif
/*Size of event data buffer*/
#if !(defined(EVENT_DATA_BUFFER_SIZE))
  #define EVENT_DATA_BUFFER_SIZE                              (24*1024)
#endif
/*The number of spike packets*/
#if !(defined(SPIKE_PACKETS))
  #define SPIKE_PACKETS                                       128 /*128*/
#endif
/*Solver error tolerance mode (see usage below)*/
#if !(defined(TOLERANCE_MODE))
  #define TOLERANCE_MODE                                      2 /*0,1,2*/
#endif
/*Max number of time slots (define the longest delay)*/
#define EVENT_TIME_SLOTS                                      16
/*Minimum event delay latency*/
#define MINIMUM_PROPAGATION_DELAY                             1.0f
/*Populate network state with user defined parameters. Useful to start spiking right away.*/
#if !(defined(PREINITIALIZE_NETWORK_STATE))
  #define PREINITIALIZE_NETWORK_STATE                         1
#endif
#if !(defined(PREINITIALIZE_NETWORK_TIME_SLOT_DELTA))
  #define PREINITIALIZE_NETWORK_TIME_SLOT_DELTA               50
#endif
#if !(defined(PREINITIALIZE_NETWORK_TIME_SLOT_DELTA_DEVIATION))
  #define PREINITIALIZE_NETWORK_TIME_SLOT_DELTA_DEVIATION     25.0
#endif
#if !(defined(PREINITIALIZE_NETWORK_PERCENT_INHIBITORY))
  #define PREINITIALIZE_NETWORK_PERCENT_INHIBITORY            10.0
#endif
#if !(defined(PREINITIALIZE_NETWORK_MIN_SPIKE_PERCENT))
  #define PREINITIALIZE_NETWORK_MIN_SPIKE_PERCENT             0.0
#endif
#if !(defined(PREINITIALIZE_NETWORK_MAX_SPIKE_PERCENT))
  #define PREINITIALIZE_NETWORK_MAX_SPIKE_PERCENT             5.0
#endif
#if !(defined(PREINITIALIZE_NETWORK_NEURON_VARIABLES_SAMPLE_FILE))
  #define PREINITIALIZE_NETWORK_NEURON_VARIABLES_SAMPLE_FILE  "neuron_variables_sample.csv"
#endif
/*Enable overwriting spike packet data until defined simulation step. Disabled if 0.
  (useful for initiating gradually increasing spiking activity during initial steps) */
#if PREINITIALIZE_NETWORK_STATE
  #undef OVERWRITE_SPIKES_UNTILL_STEP
  #define OVERWRITE_SPIKES_UNTILL_STEP                        0
#else
  #if !(defined(OVERWRITE_SPIKES_UNTILL_STEP))
    #define OVERWRITE_SPIKES_UNTILL_STEP                      (2*EVENT_TIME_SLOTS)
  #endif
#endif
/*Spike buffer occupancy bounds during overwriting spike packet data*/
#define OVERWRITE_SPIKES_MIN_MAX_PERCENT                      0.0, 0.9
/*Enable injecting current until defined simulation step. Disabled if 0.
  (useful for initiating abruptly increasin spiking activity during initial steps) */
#if PREINITIALIZE_NETWORK_STATE
  #undef INJECT_CURRENT_UNTILL_STEP
  #define INJECT_CURRENT_UNTILL_STEP                          0
#else 
  #if !(defined(INJECT_CURRENT_UNTILL_STEP))
    #define INJECT_CURRENT_UNTILL_STEP                        0
  #endif
#endif
/*Start time profiling simulation at step*/
#if !(defined(START_PROFILING_AT_STEP))
  #if OVERWRITE_SPIKES_UNTILL_STEP < INJECT_CURRENT_UNTILL_STEP
    #define START_PROFILING_AT_STEP                           (INJECT_CURRENT_UNTILL_STEP)
  #else
    #define START_PROFILING_AT_STEP                           (OVERWRITE_SPIKES_UNTILL_STEP)
  #endif
#endif
/*Step size*/
#define SIMULATION_STEP_SIZE                                  1
/*Floating point standard*/
#define DATA_TYPE                                             float
#define CL_DATA_TYPE                                          cl_float
/*Enable error try and catch if 1 or disable if 0*/
#define ERROR_TRY_CATCH_ENABLE                                0
/*Ignore exceptions related to solver math and relay on their handling by solver routine:
  PS divergence
  PS order overflow
  Detection of more than a single spike per step
  NR divergence
  NR order overflow*/
#define IGNORE_SOLVER_EXCEPTIONS                              1
/**************************************************************************************************/



/***************************************************************************************************
  Simulation Mode Configuration
***************************************************************************************************/
/*Run simulation for this many steps*/
#if !(defined(SIMULATION_TIME_STEPS))
  #define SIMULATION_TIME_STEPS                             (5*EVENT_TIME_SLOTS)
#endif

/*Simulation mode:
  0 - Detailed verification of every kernel and unit tests. For unit tests set a single bit in 
      ENABLE_MASK. Useful in debugging and thorough verification.
  1 - Light high-level verification of the whole simulation. The verification takes much less
      time than mode 0, but still provides verification coverage.
  2 - Application-level time profiling with verification of results between host and device during
      pre-profiling steps
  3 - Application-level time profiling without verification and with relaxed math.
  4 - Kernel-level time profiling without verification and with relaxed math.
  5 - A run with non-blocking return (useful for APP Profiler)
*/

#if !(defined(SIMULATION_MODE))
  #define SIMULATION_MODE                                     3
#endif

#if SIMULATION_MODE == 0                             /*Simulation mode: verification and unit test*/
  #if !(defined(ENABLE_MASK))
    //#define ENABLE_MASK                                       BIN_16(1000,0000,0000,0000)
    #define ENABLE_MASK                                       BIN_16(1111,1111,1000,0000)
  #endif
  /*Functional verification of each kernel*/
  #define KERNEL_VERIFY_ENABLE                                1
  /*Functional high-level verification of simulation results at the end of step. It is done every
    KERNEL_ENDSTEP_VERIFY_EVERY_STEPS's step*/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   1
  #endif
  /*Enable high-level verification of sort results*/
  #define SORT_VERIFY_ENABLE                                  1
  /*Enable high-level verification of model variables, spikes, events for whole network*/
  #define NETWORK_VERIFY_ENABLE                               1
  /*Log model variables with parameters LOG_MODEL_VARIABLES_* defined above*/
  #define LOG_MODEL_VARIABLES                                 1
  /*Log simulation messages*/
  #define LOG_SIMULATION                                      1
  /*Enable recording error codes from the kernels*/
  #define ERROR_TRACK_ENABLE                                  1
  /*Access error codes from the kernels every: */
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      1
  #endif
  /*Debug mask, enables each kernel to have debug buffer r/w*/
  #define DEBUG_MASK                                          1
  /*Enable compiler flags that allow device to produce same result as host*/
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  /*Profiling mode*/
  #define PROFILING_MODE                                      0
  /*Enable gathering statistics*/
  #define STATISTICS_ENABLE                                   1

#elif SIMULATION_MODE == 1
  /*Has to be a full mask for this mode*/
  #undef  ENABLE_MASK
  #define ENABLE_MASK                                         BIN_16(1111,1111,1000,0000)
  #define KERNEL_VERIFY_ENABLE                                0
  #define SORT_VERIFY_ENABLE                                  0
  #define LOG_MODEL_VARIABLES                                 0
  #define LOG_SIMULATION                                      0
  #define DEBUG_MASK                                          0
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      0
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  #define NETWORK_VERIFY_ENABLE                               1
  #define ERROR_TRACK_ENABLE                                  1
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   1
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      1
  #endif
  
#elif SIMULATION_MODE == 2
  #if !(defined(ENABLE_MASK))
    //#define ENABLE_MASK                                       BIN_16(1000,0000,0000,0000)
    #define ENABLE_MASK                                       BIN_16(1111,1111,1000,0000)
  #endif
  #define KERNEL_VERIFY_ENABLE                                0
  #define SORT_VERIFY_ENABLE                                  0
  #define LOG_MODEL_VARIABLES                                 0
  #define LOG_SIMULATION                                      0
  #define DEBUG_MASK                                          0
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      1
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  #define NETWORK_VERIFY_ENABLE                               1
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  
#elif SIMULATION_MODE == 3
  #if !(defined(ENABLE_MASK))
    //#define ENABLE_MASK                                       BIN_16(1000,0000,0000,0000)
    #define ENABLE_MASK                                       BIN_16(1111,1111,1000,0000)
  #endif
  #define KERNEL_VERIFY_ENABLE                                0
  #define SORT_VERIFY_ENABLE                                  0
  #define LOG_MODEL_VARIABLES                                 0
  #define LOG_SIMULATION                                      0
  #define DEBUG_MASK                                          0
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      1
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      1
  #define NETWORK_VERIFY_ENABLE                               0
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  
#elif SIMULATION_MODE == 4
  #if !(defined(ENABLE_MASK))
    //#define ENABLE_MASK                                       BIN_16(1000,0000,0000,0000)
    #define ENABLE_MASK                                       BIN_16(1111,1111,1000,0000)
  #endif
  #define KERNEL_VERIFY_ENABLE                                0
  #define SORT_VERIFY_ENABLE                                  0
  #define LOG_MODEL_VARIABLES                                 0
  #define LOG_SIMULATION                                      0
  #define DEBUG_MASK                                          0
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      2
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      1
  #define NETWORK_VERIFY_ENABLE                               0
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  
#elif SIMULATION_MODE == 5
  #if !(defined(ENABLE_MASK))
    //#define ENABLE_MASK                                       BIN_16(1000,0000,0000,0000)
    #define ENABLE_MASK                                       BIN_16(1111,1111,1000,0000)
  #endif
  #define KERNEL_VERIFY_ENABLE                                0
  #define SORT_VERIFY_ENABLE                                  0
  #define LOG_MODEL_VARIABLES                                 0
  #define LOG_SIMULATION                                      0
  #define DEBUG_MASK                                          0
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      0
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  #define NETWORK_VERIFY_ENABLE                               0
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  
#else
  #error Unknown simulation mode
#endif

/*
  Enable simulation stages (use ENABLE_MASK above)
  Options: 
            - any single bit is set: unit test for that stage
            - some contiguios set of bits starting from MSB: either unit tests or integration tests 
              for enabled stages depending on the stage combination
            - all bits are set: complete integration test
*/
#define EXPAND_EVENTS_ENABLE                                  (ENABLE_MASK&32768)
#define SCAN_ENABLE_V00                                       (ENABLE_MASK&16384)
#define GROUP_EVENTS_ENABLE_V00                               (ENABLE_MASK&8192)
#define SCAN_ENABLE_V01                                       (ENABLE_MASK&4096)
#define GROUP_EVENTS_ENABLE_V01                               (ENABLE_MASK&2048)
#define GROUP_EVENTS_ENABLE_V02                               (ENABLE_MASK&1024)
#define GROUP_EVENTS_ENABLE_V03                               (ENABLE_MASK&512)
#define MAKE_EVENT_PTRS_ENABLE                                (ENABLE_MASK&264)
#define UPDATE_NEURONS_ENABLE_V00                             (ENABLE_MASK&128)
/**************************************************************************************************/



/***************************************************************************************************
  Simulation parameters for the CPU implementation used as a verification reference
***************************************************************************************************/
/*Synaptic event queue size limit per nrn: */
#define REFERENCE_EVENT_QUEUE_SIZE                            1000
#define FLIXIBLE_DELAYS_ENABLE                                1
/**************************************************************************************************/



/***************************************************************************************************
  Misceleneous definitions
***************************************************************************************************/
/*String used to tag statistics relevant to all kernels*/
#define KERNEL_ALL                                            "All Kernels"
#define OCL_COMPILER_OPTIONS_FILE_NAME                        "oclCompilerOptions.txt"
/**************************************************************************************************/



/***************************************************************************************************
  General Restrictions
***************************************************************************************************/
#if SORT_VERIFY_ENABLE && \
    !(GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 && \
     GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE)
  #undef SORT_VERIFY_ENABLE
  #define SORT_VERIFY_ENABLE  0
#endif
#if (PROFILING_MODE == 1) && KERNEL_VERIFY_ENABLE
  #error (cannot use (PROFILING_MODE = 1) if KERNEL_VERIFY_ENABLE is true)
#endif
/**************************************************************************************************/



/** ############################################################################################# **

  II. Kernel-specific Parameters, Macros and Definitions
  
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
  
** ############################################################################################# **/



/***************************************************************************************************

  Configuration interface with propagation kernel for synaptic events
  
***************************************************************************************************/
#if EXPAND_EVENTS_ENABLE
  /*Verification, Testing*/
  /*CONTROL: enable verification*/
  #define EXPAND_EVENTS_VERIFY_ENABLE                         KERNEL_VERIFY_ENABLE
  /*CONTROL: reset data each EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT time step*/
  #define EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT               17
  
  /*Debugging*/
  /*CONTROL: enable debugging*/
  #define EXPAND_EVENTS_DEBUG_ENABLE                          (0)&DEBUG_MASK
  /*CONTROL: debug buffer size*/
  #define EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS               (1024*1024)
  
  /*Error tracking and codes*/
  /*CONTROL: enable error tracking*/
  #define EXPAND_EVENTS_ERROR_TRACK_ENABLE                    ERROR_TRACK_ENABLE
  #define EXPAND_EVENTS_ERROR_BUFFER_SIZE_WORDS               1
  #define EXPAND_EVENTS_ERROR_CODE_1                          0x1 /*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE overflow*/
  
  /*CONTROL: the number of spike packets*/
  #define EXPAND_EVENTS_SPIKE_PACKETS                         (SPIKE_PACKETS)
  #define EXPAND_EVENTS_SPIKE_PACKETS_PER_WG                  1
  
  /*Kernel parameters*/
  #define EXPAND_EVENTS_KERNEL_FILE_NAME                      "Kernel_ExpandEvents.cl"
  #define EXPAND_EVENTS_KERNEL_NAME                           "expand_events"
  #define EXPAND_EVENTS_WF_SIZE_WI                            WF_SIZE_WI
  /*CONTROL: number of WFs in a WG*/
  #if !(defined(EXPAND_EVENTS_WG_SIZE_WF))
    #define EXPAND_EVENTS_WG_SIZE_WF                          4 /*Options: 1, 2, 4*/
  #endif
  #define EXPAND_EVENTS_WG_SIZE_WI                            (EXPAND_EVENTS_WG_SIZE_WF*\
                                                              EXPAND_EVENTS_WF_SIZE_WI)
  #define EXPAND_EVENTS_GRID_SIZE_WG                          (EXPAND_EVENTS_SPIKE_PACKETS/\
                                                              EXPAND_EVENTS_SPIKE_PACKETS_PER_WG)
  
  /*Spike data structure parameters (data in)*/
  #define EXPAND_EVENTS_MIN_MAX_SPIKE_PERCENT                 0.0, 10.0
  #define EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE                SPIKE_DATA_BUFFER_SIZE
  #define EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE              2
  #define EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS            2
  #define EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS              (EXPAND_EVENTS_SPIKE_TOTALS_BUFFER_SIZE + \
                                                              EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE *\
                                                              EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS)

  /*Time slot (bin) buffer parameters (data out)*/
  /*CONTROL: Max number of time slots (define the longest delay)*/
  #define EXPAND_EVENTS_TIME_SLOTS                            EVENT_TIME_SLOTS
  #define EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS                EXPAND_EVENTS_GRID_SIZE_WG
  #define EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS            3
  #define EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE   EVENT_DATA_BUFFER_SIZE
  
  /*NN parameters*/
  /*CONTROL: Network size in terms of bits (1<<EXPAND_EVENTS_TOTAL_NEURON_BITS)*/
  #define EXPAND_EVENTS_TOTAL_NEURON_BITS                     TOTAL_NEURON_BITS
  #define EXPAND_EVENTS_TOTAL_NEURONS                         (1<<EXPAND_EVENTS_TOTAL_NEURON_BITS)
  #define EXPAND_EVENTS_MAX_DELAY                             (EXPAND_EVENTS_TIME_SLOTS-\
                                                              SIMULATION_STEP_SIZE)
                                                              
  /*Synaptic data structure parameters*/
  /*CONTROL: Max number of synapses per neuron*/
  #define EXPAND_EVENTS_MAX_SYNAPTIC_DATA_SIZE                (MAX_SYNAPSES_PER_NEURON)
  #define EXPAND_EVENTS_SYNAPTIC_POINTER_SIZE                 (EXPAND_EVENTS_TOTAL_NEURONS+1)
  
  /*Histogram of target neurons*/
  #define EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM               1
  #define EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT                   0
  #define EXPAND_EVENTS_HISTOGRAM_BIN_BITS                    4 /*Must be aligned with EXPAND_EVENTS_HISTOGRAM_BIN_MASK*/
  #define EXPAND_EVENTS_HISTOGRAM_BIN_MASK                    0xF /*Must be aligned with EXPAND_EVENTS_HISTOGRAM_BIN_BITS*/
  #define EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS                  (1<<EXPAND_EVENTS_HISTOGRAM_BIN_BITS)
  
  /*CONTROL: Minimum event delay latency*/
  #define EXPAND_EVENTS_MIN_DELAY                             (MINIMUM_PROPAGATION_DELAY)
#endif
/**************************************************************************************************/



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
  #define SCAN_DEBUG_ENABLE                                   (0)&DEBUG_MASK
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
/**************************************************************************************************/



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
  #define GROUP_EVENTS_DEBUG_ENABLE                             (0)&DEBUG_MASK
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
  #define GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE           EVENT_DATA_BUFFER_SIZE
  #define GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE           (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS*\
                                                                GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE)
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
  /*Number of sort iterations for neuron ID depends on the max ID*/
#define GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS \
                              ((GROUP_EVENTS_TOTAL_NEURON_BITS/GROUP_EVENTS_HISTOGRAM_BIN_BITS) +\
                              ((GROUP_EVENTS_TOTAL_NEURON_BITS%GROUP_EVENTS_HISTOGRAM_BIN_BITS) > 0))
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
/**************************************************************************************************/



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
  #define MAKE_EVENT_PTRS_DEBUG_ENABLE                             (0)&DEBUG_MASK
  #define MAKE_EVENT_PTRS_DEBUG_BUFFER_SIZE_WORDS                  (128*MAKE_EVENT_PTRS_WG_SIZE_WI*MAKE_EVENT_PTRS_GRID_SIZE_WG)
  
  /*Error tracking*/
  /*CONTROL: enable error logging at kernel level*/
  #define MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE                       ERROR_TRACK_ENABLE
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
  #define MAKE_EVENT_PTRS_SYNAPTIC_EVENT_BUFFERS                    128
  #if !(defined(MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE))
  #define MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE                  (MAKE_EVENT_PTRS_SYNAPTIC_EVENT_BUFFERS*\
                                                                    EVENT_DATA_BUFFER_SIZE)
  #endif
  /*CONTROL: Max neuron ID for unit test.*/
  #define MAKE_EVENT_PTRS_TEST_MAX_NEURON_ID                       (MAKE_EVENT_PTRS_TOTAL_NEURONS-1)
  
  /*CONTROL: LM allocation for computing event counts per WF*/
  #define MAKE_EVENT_PTRS_EVENT_COUNT_WF_LM_SHARE                   1024
  
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
/**************************************************************************************************/



/***************************************************************************************************

  Kernel configuration interface for updating model variables
  
***************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00 || UPDATE_NEURONS_ENABLE_V01
/*
                                               Generic parameters
*/
  /*Verification*/
  #define UPDATE_NEURONS_VERIFY_ENABLE                            KERNEL_VERIFY_ENABLE
  /*Error tracking and codes*/
  /*CONTROL: enable error logging at kernel level*/
  #define UPDATE_NEURONS_ERROR_TRACK_ENABLE                       ERROR_TRACK_ENABLE
  #define UPDATE_NEURONS_ERROR_BUFFER_SIZE_WORDS                  1
  #define UPDATE_NEURONS_ERROR_NON_SOLVER_FAILURE_MASK            0x00000012
  #define UPDATE_NEURONS_ERROR_CODE_1                             1   /*UPDATE_NEURONS_PS_ORDER_LIMIT was hit*/
  #define UPDATE_NEURONS_ERROR_CODE_2                             2   /*A neuron spiked more than once*/
  #define UPDATE_NEURONS_ERROR_CODE_3                             4   /*UPDATE_NEURONS_NR_ORDER_LIMIT was hit*/
  #define UPDATE_NEURONS_ERROR_CODE_4                             8   /*NR diverged*/
  #define UPDATE_NEURONS_ERROR_CODE_5                             16  /*Spike packet size limit was hit*/
  #define UPDATE_NEURONS_ERROR_CODE_6                             32  /*PS divergence*/
  /*Kernel file and name*/
  #define UPDATE_NEURONS_KERNEL_FILE_NAME                         "Kernel_UpdateNeurons.cl"
  #define UPDATE_NEURONS_KERNEL_NAME                              "update_neurons"
  #define UPDATE_NEURONS_SPIKED_KERNEL_NAME                       "update_spiked_neurons"
  /*WF size measured in WIs*/
  #define UPDATE_NEURONS_WF_SIZE_WI                               WF_SIZE_WI /*Options: 64*/
  /*Scaling parameter for the grid size*/
  #define UPDATE_NEURONS_WF_DATA_CHUNKS                           16
  /*Size of element in event pointer struct*/
  #define UPDATE_NEURONS_STRUCT_ELEMENT_SIZE                      2
  /*Event input buffer element size*/
  #define UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS                   3
  /*CONTROL: each WI processes this number of neurons*/
  #define UPDATE_NEURONS_ELEMENTS_PER_WI                          1
  /*Size of time delay buffer*/
  #define UPDATE_NEURONS_TIME_SLOTS                               EVENT_TIME_SLOTS /*Options: 16*/
  /*Sizes of model variables and parameters*/
  #define UPDATE_NEURONS_MODEL_VARIABLES                          4
  #define UPDATE_NEURONS_MODEL_PARAMETERS                         9
  /*Spike data structure parameters*/
  #define UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE                   SPIKE_DATA_BUFFER_SIZE
  #define UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE                 2
  #define UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS               2
  #define UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS                  (UPDATE_NEURONS_SPIKE_TOTALS_BUFFER_SIZE + \
                                                                  UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE * \
                                                                  UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS)
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
  
  /*Debugging*/
  #define UPDATE_NEURONS_DEBUG_ENABLE                             (0)&DEBUG_MASK
  #define UPDATE_NEURONS_DEBUG_BUFFER_SIZE_WORDS                  (UPDATE_NEURONS_TOTAL_NEURONS*1024)
  
  /*CONTROL: Source data size for unit test*/
  #define UPDATE_NEURONS_SYNAPTIC_EVENT_BUFFERS                   128
  #if !(defined(UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE))
  #define UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE                 (UPDATE_NEURONS_SYNAPTIC_EVENT_BUFFERS*\
                                                                  EVENT_DATA_BUFFER_SIZE)
  #endif
  /*Simulation parameters*/
  /*Tolerance mode: 
  0 - zero PS tolerance (0.0), UPDATE_NEURONS_PS_TOLERANCE is ignored,
      NR tolerance is UPDATE_NEURONS_NR_ZERO_TOLERANCE, UPDATE_NEURONS_NR_TOLERANCE is ignored
  1 - tolerance is defined for all neurons by UPDATE_NEURONS_PS_TOLERANCE and 
      UPDATE_NEURONS_NR_TOLERANCE
  2 - PS and NR tolerances are the same and defined on chunk basis in a constant data structure with 
      UPDATE_NEURONS_TOLERANCE_CHUNKS except that for neurons with PS zero tolerance (0.0) the NR
      tolerance is defined by UPDATE_NEURONS_NR_ZERO_TOLERANCE*/
  #define UPDATE_NEURONS_TOLERANCE_MODE                           TOLERANCE_MODE
  #define UPDATE_NEURONS_USE_VPEAK_FOR_DIVERGENCE_THRESHOLD       1
  #define UPDATE_NEURONS_TOLERANCE_CHUNKS                         (1<<8)
  /*Possible UPDATE_NEURONS_PS_TOLERANCE values: const double tols[16] = 
  {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,1e-16};*/
  #define UPDATE_NEURONS_PS_TOLERANCE                             (1.0e-8)
  #define UPDATE_NEURONS_NR_TOLERANCE                             (1.0e-8) /*normally no more than UPDATE_NEURONS_PS_TOLERANCE*/
  #define UPDATE_NEURONS_NR_ZERO_TOLERANCE                        (1.0e-16)
  #define UPDATE_NEURONS_PS_ORDER_LIMIT                           16
  #define UPDATE_NEURONS_NR_ORDER_LIMIT                           32
  #define UPDATE_NEURONS_INJECT_CURRENT_UNTIL_STEP                INJECT_CURRENT_UNTILL_STEP
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
  
  /*Optimizations*/
  #define UPDATE_NEURONS_DT_1_0_OPTIMIZATION                      1
  
/*
                                               Restrictions
*/
#if UPDATE_NEURONS_ELEMENTS_PER_WI != 1
  #error (UPDATE_NEURONS_ELEMENTS_PER_WI has to be 1)
#endif
#if	((UPDATE_NEURONS_TOLERANCE_MODE == 2) && \
    (UPDATE_NEURONS_TOTAL_NEURONS < UPDATE_NEURONS_TOLERANCE_CHUNKS))
  #error (UPDATE_NEURONS_TOTAL_NEURONS has to be no less than UPDATE_NEURONS_TOLERANCE_CHUNKS)
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
/**************************************************************************************************/

#endif /*INC_DEFINITIONS_H*/
