
#ifndef INC_DEFINITIONS_H
#define INC_DEFINITIONS_H



/**
  @file Definitions.h

  Definitions and macros visible in all source files including kernels

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2011/09/17
*/



/** ############################################################################################# **

  Optimization guidelines
  
  SI
    Description: 
      Scalar CUs, 4 vector ALUs and 1 scalar ALU / CU. Vector instruction executs 4x longer
      than scalar. There are 512 SGPRs and 256 VGPRs / Vector ALU. SGPR allocation is in 
      increments of 8, and VGPR allocation is in increments of 4. Scalar instructions are branch,
      constant cache fetch. LDS: 64 KB / CU, max allocation is 32 KB / WG, allocation granularity 
      is 256 B. L1: R/W, 16 KB / CU.
    
    Optimization tips:
      - Have at least 1 WF/CU. Max is 10. Try to get 10.
      - Max optimal Scalar/Vector instruction ratio = 1.
      - Read coalescing does not work for 64-bit data sizes.
      - WGs with 256 WIs ensure that each CU is being used.
      - Overlap execution of different kernels where possible. Use multiple buffers.
      
    TODO:
      replace SYSTEM def with SYSTEM
  
** ############################################################################################# **/



/** ############################################################################################# **

  I. Compiler-dependent control

** ############################################################################################# **/



/***************************************************************************************************
  System encoding: COMPILER_OS_NAME = CON, where
    C - control encoding (range 1 - F)
      1 - compiler control is on
      F - compiler control is off
    O - OS encoding (range 0 - FF):
      0 - Windows
      1 - Linux
    N - compiler name encoding (range 0 - FF):
      0 - Microsoft Visual C++
      1 - MinGW64
      2 - GNU g++
***************************************************************************************************/

/*Masks*/
#define SYSTEM_CONTROL_MASK         0xF0000
#define SYSTEM_OS_MASK              0x0FF00
#define SYSTEM_COMPILER_MASK        0x000FF

/*System types*/
#define SYSTEM_WIN_VC               0x10000
#define SYSTEM_WIN_MINGW64          0x10001
#define SYSTEM_LYN_GPP              0x10102
#define SYSTEM_CONTROL_OFF          0xF0000

/*Control types*/
#define SYSTEM_COMPILER_CONTROL_ON  0x1
#define SYSTEM_COMPILER_CONTROL_OFF 0xF

/*OS types*/
#define SYSTEM_OS_WINDOWS           0x0
#define SYSTEM_OS_LINUX             0x1

/*Compiler types*/
#define SYSTEM_COMPILER_MSVC        0x0
#define SYSTEM_COMPILER_MINGW64     0x1
#define SYSTEM_COMPILER_GPP         0x2

/**************************************************************************************************/



/***************************************************************************************************
  System identification and dependencies
***************************************************************************************************/

/*Intel compiler specific stuff
#ifndef _WIN32
#if defined(__INTEL_COMPILER)
  #pragma warning(disable : 1125)
#endif
#endif*/

#if !defined (SYSTEM)
  #error (SYSTEM flag is not defined)
#endif

/*Decode system:*/
#define SYSTEM_CONTROL    ((SYSTEM_CONTROL_MASK & SYSTEM) >> (4*4))
#define SYSTEM_OS         ((SYSTEM_OS_MASK & SYSTEM) >> (4*2))
#define SYSTEM_COMPILER   (SYSTEM_COMPILER_MASK & SYSTEM)

/*Validate system:*/
#if ((SYSTEM_OS != SYSTEM_OS_WINDOWS) && (SYSTEM_OS != SYSTEM_OS_LINUX))
  #error (SYSTEM_OS is not supported)
#endif
/* *** */
#if ((SYSTEM_COMPILER != SYSTEM_COMPILER_MSVC) && (SYSTEM_COMPILER != SYSTEM_COMPILER_MINGW64) && \
     (SYSTEM_COMPILER != SYSTEM_COMPILER_GPP))
  #error (SYSTEM_COMPILER is not supported)
#endif

/*Configure system:*/
#if SYSTEM == SYSTEM_WIN_VC

  #define WARNING_CONTROL_START \
    __pragma(warning( push, 4 ))  /*\
    __pragma(warning(default : 4710 4711 4820)) */
    
  #define WARNING_CONTROL_END \
    __pragma(warning( pop ))  /*\
    __pragma(warning(disable : 4710 4711 4820)) */
    
  #define WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START
  #define WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END
  
#elif ((SYSTEM == SYSTEM_WIN_MINGW64) || (SYSTEM == SYSTEM_LYN_GPP))
  /*mingw compiler specific stuff
  #if defined(__MINGW32__) && !defined(__MINGW64_VERSION_MAJOR)
    #define _aligned_malloc __mingw_aligned_malloc 
    #define _aligned_free  __mingw_aligned_free 
  #endif*/

  #define WARNING_CONTROL_START \
    _Pragma("GCC diagnostic push")  \
    _Pragma("GCC diagnostic error \"-Wall\"") \
    _Pragma("GCC diagnostic error \"-Wextra\"") \
    _Pragma("GCC diagnostic error \"-Wtraditional\"") \
    _Pragma("GCC diagnostic error \"-Wunreachable-code\"")

  #define WARNING_CONTROL_END \
    _Pragma("GCC diagnostic pop")  
    
  #define WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START \
    _Pragma("GCC diagnostic ignored \"-Wfloat-equal\"")  
    
  #define WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END \
    _Pragma("GCC diagnostic error \"-Wfloat-equal\"")  
    
#elif (SYSTEM_CONTROL == SYSTEM_COMPILER_CONTROL_OFF)
  #define WARNING_CONTROL_START
  #define WARNING_CONTROL_END
  #define WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START
  #define WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END
#else
  #error (Unsupported SYSTEM type)
#endif

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



/** ############################################################################################# **

  II. Configuration space. Use this space to overwrite default definitions.

** ############################################################################################# **/



/** PREPEND CUE (for compilation script, do not remove)**/



/** ############################################################################################# **

  III. Independent Parameters, Macros and Definitions
  
  (Note: some of them have dependencies within this section, which are resolved by in-order 
  definitions)
  
** ############################################################################################# **/



/***************************************************************************************************
  String operations
***************************************************************************************************/

#define STRINGIFY(x)    #x
#define TOSTRING(x)     STRINGIFY(x)
    
/**************************************************************************************************/



/***************************************************************************************************
  Exceptions
***************************************************************************************************/

#define CATCH_GENERIC(stream, message, action)\
  catch(exception& e)\
  {\
    stream << #message << ": Standard Exception: " << e.what() << "\n";\
    action;\
  }\
  catch (...)\
  {\
    stream << #message << ": unknown exception occured" << "\n";\
    action;\
  }
  
#define CATCH(stream, message, action)\
  catch(SimException& e)\
  {\
    stream << #message << ": " << e.what() << "\n";\
    action;\
  }\
  CATCH_GENERIC(stream, message, action)
  
#define THROW_SIMEX(msg)\
  {\
    std::stringstream ss;\
    ss << msg << std::endl;\
    throw SimException(ss.str());\
  }
    
/**************************************************************************************************/



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
  
/**************************************************************************************************/



/***************************************************************************************************
  Random number generators
***************************************************************************************************/
#define SET_RANDOM_SEED_WITH_TIME(seed, counter)\
  {\
    time_t t = time(NULL);\
    seed = *((unsigned int *)(&t)) + counter;\
    counter++;\
    srand(seed);\
  }
  
#define SET_RANDOM_SEED_DIRECT(seed)\
  {\
    srand(seed);\
  }
  
#define SET_RANDOM_SEED_NONE(seed)
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
  Error checkers
***************************************************************************************************/

#define ASSERT_CL_SUCCESS(tested_variable, message)\
  if(tested_variable != CL_SUCCESS)\
  {\
    throw SimException(message);\
  }

/**************************************************************************************************/



/***************************************************************************************************
  Profiling
***************************************************************************************************/

/*Time registration for statistics*/
#define REGISTER_TIME(kernel_name, time, increment, stats_struct)\
  {\
    map<std::string, double> stats = stats_struct.execTime[kernel_name];\
    if(stats.find("Time") ==  stats.end()){stats["Time"] = 0;}\
    if(stats.find("Count") ==  stats.end()){stats["Count"] = 0;}\
    stats["Time"] += time;\
    stats["Count"] += increment;\
    stats_struct.execTime[kernel_name] = stats;\
    set<std::string> kernels = stats_struct.kernelNamesExecTime;\
    kernels.insert(kernel_name);\
    stats_struct.kernelNamesExecTime = kernels;\
  }
  
/*Get time stamp in ns*/
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
#define GET_TIME_NS(time)\
  {\
    QueryPerformanceCounter((LARGE_INTEGER *)&(this->performanceCounter));\
    time = ((double)(this->performanceCounter) / (double)(this->performanceFrequency));\
  }
#elif (SYSTEM_OS == SYSTEM_OS_LINUX)
#define GET_TIME_NS(time)\
  {\
    struct timespec tp;\
    clock_gettime(CLOCK_MONOTONIC, &tp);\
    time = (double) tp.tv_sec * (1000ULL * 1000ULL * 1000ULL) + (double) tp.tv_nsec;\
  }
#endif
  
/*Start timer*/
#define START_TIMER(step, starting_step, startTime)\
  {\
    if(step >= starting_step){GET_TIME_NS(startTime);}\
  }

/*End timer*/
#define END_TIMER(step, starting_step, event, tag, startTime, stats_struct)\
  if(step >= starting_step)\
  {\
    cl_int status = event.wait();\
    double endAppTime;\
    GET_TIME_NS(endAppTime);\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("END_TIMER: Failed cl:Event.wait() for " << tag << " due to error code " \
        << status);\
    }\
    REGISTER_TIME(tag, (endAppTime-startTime), 1.0, stats_struct)\
  }
  
/**************************************************************************************************/



/***************************************************************************************************
  Memory allocation and registration
***************************************************************************************************/

/*Memory types*/
#define MEM_CONSTANT                                          0
#define MEM_GLOBAL                                            1
#define MEM_LOCAL                                             2

/*Memory size and name registration for statistics*/
#define REGISTER_MEMORY(kernel_name, mem_type, mem_name)\
  {\
    cl_uint minDataTypeAlignSize = (device).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>();\
    cl_ulong memMaxAllocactionSize = (device).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();\
    size_t size = ((mem_name ##SizeBytes)/minDataTypeAlignSize + 1)*minDataTypeAlignSize;\
    \
    switch (mem_type)\
    {\
      case MEM_CONSTANT:\
      {\
        cl_ulong maxConstantMemSize = (device).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>();\
        map<std::string, size_t> memSizes = kernelStats.cmSizes[kernel_name];\
        memSizes[#mem_name] = size;\
        memSizes["TOTAL"] += size;\
        kernelStats.cmSizes[kernel_name] = memSizes;\
        if(memSizes["TOTAL"] > maxConstantMemSize)\
        {\
          THROW_SIMEX("Total allocation for constant memory in kernel " << kernel_name\
            << ", " << ((float)memSizes["TOTAL"])/(1024.0)\
            << " KB, exceeds CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "\
            << ((float)maxConstantMemSize)/(1024.0) << " KB");\
        }\
      }\
        break;\
      case MEM_GLOBAL:\
      {\
        cl_ulong maxGlobalMemSize = (device).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();\
        map<std::string, size_t> memSizes = kernelStats.gmSizes[kernel_name];\
        memSizes[#mem_name] = size;\
        memSizes["TOTAL"] += size;\
        kernelStats.gmSizes[kernel_name] = memSizes;\
        if(memSizes["TOTAL"] > maxGlobalMemSize)\
        {\
          THROW_SIMEX("Total allocation for global memory in kernel " << kernel_name\
            << ", " << ((float)memSizes["TOTAL"])/(1024.0*1024.0)\
            << " MB, exceeds CL_DEVICE_GLOBAL_MEM_SIZE, "\
            << ((float)maxGlobalMemSize)/(1024.0*1024.0) << "\n");\
        }\
      }\
        break;\
      case MEM_LOCAL:\
      {\
        cl_ulong maxLocalMemSize = (device).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();\
        map<std::string, size_t> memSizes = kernelStats.lmSizes[kernel_name];\
        memSizes[#mem_name] = size;\
        memSizes["TOTAL"] += size;\
        kernelStats.lmSizes[kernel_name] = memSizes;\
        if(memSizes["TOTAL"] > maxLocalMemSize)\
        {\
          THROW_SIMEX("Total allocation for local memory in kernel " << kernel_name\
            << ", " << ((float)memSizes["TOTAL"])/(1024.0)\
            << " KB, exceeds CL_DEVICE_LOCAL_MEM_SIZE, "\
            << ((float)maxLocalMemSize)/(1024.0) << " KB");\
        }\
      }\
        break;\
      default:\
        std::cout << "REGISTER_MEMORY: unsupported memory type\n" << std::endl;\
        throw SimException("REGISTER_MEMORY: unsupported memory type.");\
    }\
    set<std::string> kernels = kernelStats.kernelNames;\
    kernels.insert(kernel_name);\
    kernelStats.kernelNames = kernels;\
    \
    if(size > memMaxAllocactionSize)\
    {\
      THROW_SIMEX("REGISTER_MEMORY: Memory object in kernel " << kernel_name\
        << " with size identifier " << #mem_name << " and size "\
        << ((float)size)/(1024.0*1024.0) << " MB exceeds CL_DEVICE_MAX_MEM_ALLOC_SIZE, "\
        << ((float)memMaxAllocactionSize)/(1024.0*1024.0));\
    }\
  }
  
#define REGISTER_MEMORY_O(device, kernel_name, mem_type, mem_name, kernelStats)\
  {\
  size_t size_test = mem_name ##SizeBytes;\
  \
  if( size_test > 0)\
  {\
    cl_uint minDataTypeAlignSize = (device).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>();\
    cl_ulong memMaxAllocactionSize = (device).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();\
    size_t size = ((mem_name ##SizeBytes)/minDataTypeAlignSize + 1)*minDataTypeAlignSize;\
    \
    switch (mem_type)\
    {\
      case MEM_CONSTANT:\
      {\
        cl_ulong maxConstantMemSize = (device).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>();\
        map<std::string, size_t> memSizes = kernelStats.cmSizes[kernel_name];\
        memSizes[#mem_name] = size;\
        memSizes["TOTAL"] += size;\
        kernelStats.cmSizes[kernel_name] = memSizes;\
        if(memSizes["TOTAL"] > maxConstantMemSize)\
        {\
         THROW_SIMEX("Total allocation for constant memory in kernel " << kernel_name\
            << ", " << ((float)memSizes["TOTAL"])/(1024.0)\
            << " KB, exceeds CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "\
            << ((float)maxConstantMemSize)/(1024.0) << " KB");\
        }\
      }\
        break;\
      case MEM_GLOBAL:\
      {\
        cl_ulong maxGlobalMemSize = (device).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();\
        map<std::string, size_t> memSizes = kernelStats.gmSizes[kernel_name];\
        memSizes[#mem_name] = size;\
        memSizes["TOTAL"] += size;\
        kernelStats.gmSizes[kernel_name] = memSizes;\
        if(memSizes["TOTAL"] > maxGlobalMemSize)\
        {\
          THROW_SIMEX("Total allocation for global memory in kernel " << kernel_name\
            << ", " << ((float)memSizes["TOTAL"])/(1024.0*1024.0)\
            << " MB, exceeds CL_DEVICE_GLOBAL_MEM_SIZE, "\
            << ((float)maxGlobalMemSize)/(1024.0*1024.0));\
        }\
      }\
        break;\
      case MEM_LOCAL:\
      {\
        cl_ulong maxLocalMemSize = (device).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();\
        map<std::string, size_t> memSizes = kernelStats.lmSizes[kernel_name];\
        memSizes[#mem_name] = size;\
        memSizes["TOTAL"] += size;\
        kernelStats.lmSizes[kernel_name] = memSizes;\
        if(memSizes["TOTAL"] > maxLocalMemSize)\
        {\
          THROW_SIMEX("Total allocation for local memory in kernel " << kernel_name\
            << ", " << ((float)memSizes["TOTAL"])/(1024.0)\
            << " KB, exceeds CL_DEVICE_LOCAL_MEM_SIZE, "\
            << ((float)maxLocalMemSize)/(1024.0) << " KB");\
        }\
      }\
        break;\
      default:\
        throw SimException("REGISTER_MEMORY_O: unsupported memory type.");\
    }\
    \
    set<std::string> kernels = kernelStats.kernelNames;\
    kernels.insert(kernel_name);\
    kernelStats.kernelNames = kernels;\
    \
    if(size > memMaxAllocactionSize)\
    {\
      THROW_SIMEX("REGISTER_MEMORY_O: Memory object in kernel " << kernel_name\
        << " with size identifier " << #mem_name << " and size "\
        << ((float)size)/(1024.0*1024.0) << " MB exceeds CL_DEVICE_MAX_MEM_ALLOC_SIZE, "\
        << ((float)memMaxAllocactionSize)/(1024.0*1024.0));\
    }\
  }\
  }
  
/*Calloc host memory*/
#define CALLOC(name, type, size)\
  {\
  size_t size_test = size;\
  if(size_test > 0)\
  {\
    if(this->name != NULL)\
    {\
      throw SimException("CALLOC: Attempted to point not NULL pointer to memory space (" \
        #name ")");\
    }\
    this->name ##Size = size;\
    this->name ##SizeBytes = (this->name ##Size) * sizeof(type);\
    this->name = (type *)calloc(this->name ##Size, sizeof(type));\
    if(this->name == NULL)\
    {\
      throw SimException("CALLOC: Failed to allocate memory for (" #name ")");\
    }\
  }\
  }

/*Swap 2 memory references*/
#define SWAP_MEM(type, a, b)\
{\
  type tmp = a; a = b; b = tmp;\
}

/**************************************************************************************************/



/***************************************************************************************************
  OCL Macros and Constants
***************************************************************************************************/

#define FIND_TARGET_DEVICE(platformDevices, targetDevices, deviceIterator, found)\
  {\
    const char *tD = targetDevices;\
    found = false;\
    char *str = (char *) calloc(0xFFFF, sizeof(char));\
    strcpy(str, tD);\
    char *pch;\
    pch = strtok (str, ";");\
    while (pch != NULL)\
    {\
      for(deviceIterator = platformDevices.begin(); deviceIterator != platformDevices.end();\
      ++deviceIterator)\
      {\
        std::string deviceName = (*deviceIterator).getInfo<CL_DEVICE_NAME>();\
        if(strcmp(deviceName.c_str(), pch) == 0)\
        {\
          found = true; break;\
        }\
      }\
      if(found){break;}\
      pch = strtok (NULL, ";");\
    }\
    free(str);\
  }
    
/*Allocate buffer*/
#define CREATE_BUFFER(flags, buffer, size)\
  {\
    cl_int err = CL_SUCCESS;\
    buffer = cl::Buffer(context, flags, size, NULL, &err);\
    if(err != CL_SUCCESS)\
    {\
      THROW_SIMEX("CREATE_BUFFER: Failed to allocate " << #buffer << " due to error code " \
        << err);\
    }\
  }
  
#define CREATE_BUFFER_O(context, flags, buffer, size)\
  if(size > 0)\
  {\
    cl_int err = CL_SUCCESS;\
    buffer = cl::Buffer(context, flags, size, NULL, &err);\
    if(err != CL_SUCCESS)\
    {\
      THROW_SIMEX("CREATE_BUFFER_O: Failed to allocate " << #buffer << " due to error code " \
        << err);\
    }\
  }
  
/*Enqueu write buffer*/
#define ENQUEUE_WRITE_BUFFER(block, queue, buffer, size, data)\
  if(data != NULL)\
  {\
    cl_int status;\
    cl::Event writeEvt = NULL;\
    status = queue.enqueueWriteBuffer(buffer, block, 0, size, data, NULL, &writeEvt);\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("ENQUEUE_WRITE_BUFFER: Failed to enqueue buffer " << #buffer \
        << " for write due to error code " << status);\
    }\
    \
    status = queue.flush();\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("ENQUEUE_WRITE_BUFFER: Failed to flush buffer" << #buffer \
        << " for write due to error code " << status);\
    }\
    \
    cl_int eventStatus = CL_QUEUED;\
    while(eventStatus != CL_COMPLETE)\
    {\
      status = writeEvt.getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS,&eventStatus);\
      if(status != CL_SUCCESS)\
      {\
        THROW_SIMEX("ENQUEUE_WRITE_BUFFER: Failed to get successefull write event status "\
          << "for buffer " << #buffer << " due to error code " << status);\
      }\
    }\
  }
  
/*Enqueue read buffer*/
#define ENQUEUE_READ_BUFFER(block, queue, buffer, size, data)\
  if(data != NULL)\
  {\
    cl_int status;\
    cl_int eventStatus = CL_QUEUED;\
    cl::Event readEvt;\
    status = queue.enqueueReadBuffer(buffer, block, 0, size, data, NULL, &readEvt);\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("ENQUEUE_READ_BUFFER: Failed to enqueue buffer " << #buffer \
        << " for read due to error code " << status);\
    }\
    \
    status = queue.flush();\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("ENQUEUE_READ_BUFFER: Failed to flush buffer" << #buffer \
        << " for read due to error code " << status);\
    }\
    \
    eventStatus = CL_QUEUED;\
    while(eventStatus != CL_COMPLETE)\
    {\
      status = readEvt.getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS, &eventStatus);\
      if(status != CL_SUCCESS)\
      {\
        THROW_SIMEX("ENQUEUE_READ_BUFFER: Failed to get successefull read event status "\
          << "for buffer " << #buffer << " due to error code " << status);\
      }\
    }\
  }
  
/*Enqueue read buffer*/
#define IF_HIT_READ(select, selectHit, validMask)\
  if((select & selectHit) && !((validMask) & selectHit))
  
/*Set kernel argument*/
#define SET_KERNEL_ARG(kernel, arg, argNum)\
  {\
    cl_int status;\
    status = kernel.setArg(argNum, arg);\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("SET_KERNEL_ARG: Failed to set kernel argument #" << argNum << ", "\
      << #arg << " for kernel " << #kernel << " due to error code " << status);\
    }\
  }

#define ENQUEUE_KERNEL(queue, kernel, gThreads, lThreads, events, ndrEvt)\
  {\
    cl_int status;\
    status = queue.enqueueNDRangeKernel\
      (kernel, cl::NullRange, gThreads, lThreads, events, &ndrEvt);\
    if(status != CL_SUCCESS)\
    {\
      THROW_SIMEX("ENQUEUE_KERNEL: Failed to enqueue kernel " << #kernel << " due to error code " \
        << status)\
    }\
  }
  
/*Define OpenCL 1.2 consants*/
#if !defined (CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT)
  #define CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT                 (1 << 7)
#endif
/**************************************************************************************************/



/** ############################################################################################# **

  IV. Program Execution Definitions and Parameters
  
  (Note: they may utilize definitions from section III)
  
** ############################################################################################# **/



/***************************************************************************************************
  Device and platform parameters
***************************************************************************************************/
/*Wave size*/
#define WF_SIZE_WI                                            64
/*Target devices in the order of preference*/
#if !(defined(TARGET_DEVICE_NAME))
  #define TARGET_DEVICE_NAME                                  "Tahiti;Cayman"
#endif
/*Target Platform Vendor
#define TARGET_PLATFORM_VENDOR                                "Advanced Micro Devices, Inc."*/
/**************************************************************************************************/



/***************************************************************************************************
  Simulation parameters
***************************************************************************************************/
/*Run simulation for this many steps*/
#if !(defined(SIMULATION_TIME_STEPS))
  #define SIMULATION_TIME_STEPS                               1200
  #if !(defined(START_PROFILING_AT_STEP))
    #define START_PROFILING_AT_STEP                           200
  #endif
#endif
/*Network size in terms of bits (1<<TOTAL_NEURON_BITS)*/
#if !(defined(TOTAL_NEURON_BITS))
  #define TOTAL_NEURON_BITS                                   17
#endif
/*Max number of synapses per neuron*/
#if !(defined(MAX_SYNAPSES_PER_NEURON))
  #define MAX_SYNAPSES_PER_NEURON                             1536
#endif
/*Synapse deviation ratio from max defined by MAX_SYNAPSES_PER_NEURON*/
#if !(defined(SYNAPSE_DEVIATION_RATIO))
  #define SYNAPSE_DEVIATION_RATIO                             0.1
#endif
/*Gabaergic synapse ratio in respect to total synapses*/
#if !(defined(SYNAPSE_GABA_PERCENT))
  #define SYNAPSE_GABA_PERCENT                                0.5
#endif
/*Size of spike packet buffer (spikes) per WF*/
#if !(defined(SPIKE_PACKET_SIZE))
  #define SPIKE_PACKET_SIZE                                   16
#endif
/*Size of event data buffer*/
#if !(defined(EVENT_DATA_BUFFER_SIZE))
  #define EVENT_DATA_BUFFER_SIZE                              (64*1024)
#endif
/*The number of spike packets*/
#if !(defined(SPIKE_PACKETS))
  #define SPIKE_PACKETS                                       512
#endif
/*The size of bin for event count histogram*/
#if !(defined(SYNAPTIC_EVENT_BUFFERS))
  #define SYNAPTIC_EVENT_BUFFERS                              32
#endif
/*Solver error tolerance mode (see usage below)*/
#if !(defined(TOLERANCE_MODE))
  #define TOLERANCE_MODE                                      2
#endif
/*Max number of time slots (define the longest delay)*/
#define EVENT_TIME_SLOTS                                      16
/*Minimum event delay latency*/
#define MINIMUM_PROPAGATION_DELAY                             1.0
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
/*Enable overwriting spike packet data until defined simulation step. Disabled if -1.
  (useful for initiating gradually increasing spiking activity during initial steps) */
#if PREINITIALIZE_NETWORK_STATE
  #undef OVERWRITE_SPIKES_UNTILL_STEP
  #define OVERWRITE_SPIKES_UNTILL_STEP                        -1
#else
  #if !(defined(OVERWRITE_SPIKES_UNTILL_STEP))
    #define OVERWRITE_SPIKES_UNTILL_STEP                      (2*EVENT_TIME_SLOTS)
  #endif
#endif
/*Spike buffer occupancy bounds during overwriting spike packet data*/
#if !(defined(INITIALIZE_SPIKES_MIN_MAX_PERCENT))
  #define INITIALIZE_SPIKES_MIN_MAX_PERCENT                    10.0, 20.0, -1.0
#endif
/*Enable injecting current until defined simulation step. Disabled if -1.
  (useful for initiating abruptly increasin spiking activity during initial steps) */
#if PREINITIALIZE_NETWORK_STATE
  #undef INJECT_CURRENT_UNTILL_STEP
  #define INJECT_CURRENT_UNTILL_STEP                          -1
#else 
  #if !(defined(INJECT_CURRENT_UNTILL_STEP))
    #define INJECT_CURRENT_UNTILL_STEP                        -1
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
  
  0 - Detailed verification of every kernel and unit tests. For unit tests set a single bit in 
      ENABLE_MASK. Useful in debugging and thorough verification.
  1 - Light high-level verification of the whole simulation. The verification takes much less
      time than mode 0, but still provides verification coverage.
  2 - Application-level time profiling with verification of results between host and device during
      pre-profiling steps
  3 - Application-level time profiling without verification and with relaxed math.
  4 - Kernel-level time profiling without verification and with relaxed math.
  5 - A run with non-blocking return (useful for APP Profiler)
***************************************************************************************************/

#if !(defined(SIMULATION_MODE))
  #define SIMULATION_MODE                                     2
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
  /*Mode for generating random sequences
    0 - srand(1) is called once at the beginning of program execution. The seed is never changed.
    1 - srand(seed) is called with time() -based seed before data is generated with rand().
  */
  #if !(defined(RANDOM_GEN_MODE))
  #define RANDOM_GEN_MODE                                     1
  #endif
  /*Enable compiler flags that allow device to produce same result as host by disabling math 
    optimizations*/
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  /*Profiling mode
    0 - Off
    1 - Application-level profiling
    2 - Kernel-level profiling
  */
  #define PROFILING_MODE                                      0
  /*Enable gathering statistics*/
  #define STATISTICS_ENABLE                                   1
  /*Define Device-Host data sync mode:
    0 - sync is always off
    1 - sync is always on
    2 - sync is on only when required by specific mode
  */
  #if !(defined(DEVICE_HOST_DATA_COHERENCE))
  #define DEVICE_HOST_DATA_COHERENCE                          1
  #endif
  /*Enable validation of parameters, data ranges, etc in classes*/
  #define CLASS_VALIDATION_ENABLE                             1
  /*Define Unit Test mode:
    0 - unit testing: a unit test for operation is forced if relevant operator/kernel is enabled 
        to support it
    1 - integration testing: a unit test is enabled automatically if at least one of required 
        preceeding operators/kernels is not enabled
    2 - unit testing is restricted
  */
  #if !(defined(UNIT_TEST_MODE))
  #define UNIT_TEST_MODE                                      1
  #endif

#elif SIMULATION_MODE == 1
  /*Has to be a full mask for this mode*/
  #undef  ENABLE_MASK
  #define ENABLE_MASK                                         BIN_16(1111,1111,1000,0000)
  #define KERNEL_VERIFY_ENABLE                                0
  #define SORT_VERIFY_ENABLE                                  0
  #define LOG_MODEL_VARIABLES                                 0
  #define LOG_SIMULATION                                      1
  #define DEBUG_MASK                                          0
  /**/
  #if !(defined(RANDOM_GEN_MODE))
  #define RANDOM_GEN_MODE                                     1
  #endif
  /**/
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      0
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  #define ERROR_TRACK_ENABLE                                  1
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   1
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      1
  #endif
  /**/
  #if !(defined(DEVICE_HOST_DATA_COHERENCE))
  #define DEVICE_HOST_DATA_COHERENCE                          1
  #endif
  /**/
  #define CLASS_VALIDATION_ENABLE                             1
  /**/
  #if !(defined(UNIT_TEST_MODE))
  #define UNIT_TEST_MODE                                      1
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
  /**/
  #if !(defined(RANDOM_GEN_MODE))
  #define RANDOM_GEN_MODE                                     1
  #endif
  /**/
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      1
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  /**/
  #if !(defined(DEVICE_HOST_DATA_COHERENCE))
  #define DEVICE_HOST_DATA_COHERENCE                          2
  #endif
  /**/
  #define CLASS_VALIDATION_ENABLE                             0
  /**/
  #if !(defined(UNIT_TEST_MODE))
  #define UNIT_TEST_MODE                                      2
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
  /**/
  #if !(defined(RANDOM_GEN_MODE))
  #define RANDOM_GEN_MODE                                     1
  #endif
  /**/
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      1
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      1
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  /**/
  #if !(defined(DEVICE_HOST_DATA_COHERENCE))
  #define DEVICE_HOST_DATA_COHERENCE                          0
  #endif
  /**/
  #define CLASS_VALIDATION_ENABLE                             0
  /**/
  #if !(defined(UNIT_TEST_MODE))
  #define UNIT_TEST_MODE                                      2
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
  /**/
  #if !(defined(RANDOM_GEN_MODE))
  #define RANDOM_GEN_MODE                                     1
  #endif
  /**/
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      2
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      1
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  /**/
  #if !(defined(DEVICE_HOST_DATA_COHERENCE))
  #define DEVICE_HOST_DATA_COHERENCE                          0
  #endif
  /**/
  #define CLASS_VALIDATION_ENABLE                             0
  /**/
  #if !(defined(UNIT_TEST_MODE))
  #define UNIT_TEST_MODE                                      1
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
  /**/
  #if !(defined(RANDOM_GEN_MODE))
  #define RANDOM_GEN_MODE                                     1
  #endif
  /**/
  #define STATISTICS_ENABLE                                   0
  #define PROFILING_MODE                                      0
  #define COMPILER_FLAGS_OPTIMIZE_ENABLE                      0
  #define ERROR_TRACK_ENABLE                                  0
  /**/
  #if !(defined(KERNEL_ENDSTEP_VERIFY_EVERY_STEPS))
  #define KERNEL_ENDSTEP_VERIFY_EVERY_STEPS                   0
  #endif
  /**/
  #if !(defined(ERROR_TRACK_ACCESS_EVERY_STEPS))
  #define ERROR_TRACK_ACCESS_EVERY_STEPS                      0
  #endif
  /**/
  #if !(defined(DEVICE_HOST_DATA_COHERENCE))
  #define DEVICE_HOST_DATA_COHERENCE                          0
  #endif
  /**/
  #define CLASS_VALIDATION_ENABLE                             0
  /**/
  #if !(defined(UNIT_TEST_MODE))
  #define UNIT_TEST_MODE                                      1
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
#if !(defined(REFERENCE_EVENT_QUEUE_SIZE))
#define REFERENCE_EVENT_QUEUE_SIZE                            750
#endif
#define FLIXIBLE_DELAYS_ENABLE                                1
/**************************************************************************************************/



/***************************************************************************************************
  Misceleneous definitions
***************************************************************************************************/
/*String used to tag statistics relevant to all kernels*/
#define KERNEL_ALL                                            "All Kernels"

/*File name to put defines for OpenCL preprocessor*/
#define OCL_COMPILER_OPTIONS_FILE_NAME                        "oclCompilerOptions.txt"

/**************************************************************************************************/



/** ############################################################################################# **

  V. Dependent Parameters, Macros and Definitions 
  
  (Note: some of them depend on definitions in section IV and may use definitions from section III)
  
** ############################################################################################# **/



/***************************************************************************************************
  Exceptions
***************************************************************************************************/

#define THROW_ERROR_TRACK(specifier, step, getErrorDataCode, errorData)\
  if(!((step+1)%ERROR_TRACK_ACCESS_EVERY_STEPS))\
  {\
    getErrorDataCode\
    \
    if(errorData[0] != 0)\
    {\
      THROW_SIMEX(specifier << " received error code from the device: " << errorData[0]);\
    }\
  }
  
/**************************************************************************************************/



/***************************************************************************************************
  Data accessors
***************************************************************************************************/

#if SIMULATION_MODE > 1
  #define SAFE_GET(variable) variable
#else
  #define SAFE_GET(variable) (1 ? variable : 0)
#endif

/**************************************************************************************************/



/***************************************************************************************************
  Random number generation
***************************************************************************************************/

#if RANDOM_GEN_MODE == 0
  #define SET_RANDOM_SEED(seed, counter)      
#elif RANDOM_GEN_MODE == 1
  #define SET_RANDOM_SEED(seed, counter)      SET_RANDOM_SEED_WITH_TIME(seed, counter)
#else
  
#endif

#if CLASS_VALIDATION_ENABLE

#define GET_RANDOM_INT(setValue, max, minPercent, maxPercent)\
  {\
    if(minPercent > maxPercent)\
    {\
      throw SimException("GET_RANDOM_INT: (minPercent > maxPercent)");\
    }\
    if((minPercent > 100.0) || (maxPercent > 100.0))\
    {\
      throw SimException("GET_RANDOM_INT: ((minPercent > 100.0) || (maxPercent > 100.0))");\
    }\
    if((minPercent < 0.0) || (maxPercent < 0.0))\
    {\
      throw SimException("GET_RANDOM_INT: ((minPercent < 0.0) || (maxPercent < 0.0))");\
    }\
    if(minPercent == maxPercent)\
    {\
      setValue = cl_uint(((double)max)*(maxPercent/100.0));\
    }\
    else\
    {\
      setValue = cl_uint(((double)(max))*((minPercent/100.0) + \
        abs(((maxPercent-minPercent)/100.0)*((double)rand()/((double)RAND_MAX)))));\
    }\
  }
  
#else

#define GET_RANDOM_INT(setValue, max, minPercent, maxPercent)\
  {\
    if(minPercent == maxPercent)\
    {\
      setValue = cl_uint(((double)max)*(maxPercent/100.0));\
    }\
    else\
    {\
      setValue = cl_uint(((double)(max))*((minPercent/100.0) + \
        abs(((maxPercent-minPercent)/100.0)*((double)rand()/((double)RAND_MAX)))));\
    }\
  }
  
#endif

/**************************************************************************************************/



/***************************************************************************************************
  Profiling
***************************************************************************************************/

/*Kernel-level profiling*/
#if PROFILING_MODE == 2 && START_PROFILING_AT_STEP > -1
  #define KERNEL_LEVEL_PROFILING                              1
#else
  #define KERNEL_LEVEL_PROFILING                              0
#endif

/**************************************************************************************************/



/***************************************************************************************************
  Kernel handling
***************************************************************************************************/

/*Kernel enqueue*/
#if KERNEL_LEVEL_PROFILING
  
  #define ENQUEUE_KERNEL_V0(step, statsStruct, queue, kernel, globalThreads, localThreads, events, \
    event, tag)\
    {\
      double startAppTime = 0;\
      START_TIMER(step, START_PROFILING_AT_STEP, startAppTime);\
      \
      ENQUEUE_KERNEL(queue, kernel, globalThreads, localThreads, events, event);\
      \
      END_TIMER(step, START_PROFILING_AT_STEP, event, tag, startAppTime, statsStruct);\
    }
#else
  #define ENQUEUE_KERNEL_V0(step, statsStruct, queue, kernel, globalThreads, localThreads, events, \
    event, tag)\
    ENQUEUE_KERNEL(queue, kernel, globalThreads, localThreads, events, event);
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

/*Logging macro*/
#if LOG_REPORT
  #define LOG_REP(message)\
    {\
      *dataToReportLogFile << message << std::endl;\
    }
#else
  #define LOG_REP(message)
#endif

/*Logging macro*/
#if LOG_SIMULATION
  #define LOG_SIM(message)\
    {\
      time_t t; time (&t); char *s=ctime(&t); s[strlen(s)-1]=0;\
      *dataToSimulationLogFile << ctime(&t) << " " << message << std::endl;\
    }
#else
  #define LOG_SIM(message)
#endif
/**************************************************************************************************/



/** ############################################################################################# **

  VI. Kernel-specific Parameters, Macros and Definitions 
  
  (Note: they depend on definitions in section IV and may utilize definitions from section III)
  
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
  #if !(defined(EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT))
  #define EXPAND_EVENTS_RESET_DATA_AT_TIME_SLOT               17
  #endif
  /*CONTROL: test mode, 
    (-1) - all sizes of spike packets, 
    between 0 and 100 - up to target percent range (0 - 100)%
  */
  #if !(defined(EXPAND_EVENTS_TEST_MODE))
  #define EXPAND_EVENTS_TEST_MODE                             10
  #endif
  
  /*Debugging*/
  /*CONTROL: enable debugging*/
  #define EXPAND_EVENTS_DEBUG_ENABLE                          (0)&DEBUG_MASK
  /*CONTROL: debug buffer size*/
  #define EXPAND_EVENTS_DEBUG_BUFFER_SIZE_WORDS               (1024)
  
  /*Error tracking and codes*/
  /*CONTROL: enable error tracking*/
  #define EXPAND_EVENTS_ERROR_TRACK_ENABLE                    ERROR_TRACK_ENABLE
  #define EXPAND_EVENTS_ERROR_BUFFER_SIZE_WORDS               1
  #define EXPAND_EVENTS_ERROR_CODE_1                          0x1 /*EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE overflow*/
  #define EXPAND_EVENTS_ERROR_CODE_2                          0x2 /*EXPAND_EVENTS_SPIKE_BUFFER_SIZE overflow*/
  
  /*CONTROL: the number of spike packets*/
  #define EXPAND_EVENTS_SPIKE_PACKETS                         SPIKE_PACKETS

  /*Kernel parameters*/
  #define EXPAND_EVENTS_KERNEL_FILE_NAME                      "Kernel_ExpandEvents.cl"
  #define EXPAND_EVENTS_KERNEL_NAME                           "expand_events"
  #define EXPAND_EVENTS_WF_SIZE_WI                            WF_SIZE_WI
  /*CONTROL: number of WFs in a WG*/
  #if !(defined(EXPAND_EVENTS_WG_SIZE_WF))
    #define EXPAND_EVENTS_WG_SIZE_WF                          4 /*Options: 1, 2, 4*/
  #endif
  /*0 - WFs work independently in a WG, 1 - WFs work cooperatively*/
  #if !(defined(EXPAND_EVENTS_INTER_WF_COOPERATION))
    #if (EXPAND_EVENTS_WG_SIZE_WF > 1)
    #define EXPAND_EVENTS_INTER_WF_COOPERATION                1
    #else
    #define EXPAND_EVENTS_INTER_WF_COOPERATION                0
    #endif
  #endif
  #define EXPAND_EVENTS_WG_SIZE_WI                            (EXPAND_EVENTS_WG_SIZE_WF*\
                                                              EXPAND_EVENTS_WF_SIZE_WI)
  #if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
    #define EXPAND_EVENTS_GRID_SIZE_WG                        (SYNAPTIC_EVENT_BUFFERS/\
                                                              EXPAND_EVENTS_WG_SIZE_WF)
  #else
    #define EXPAND_EVENTS_GRID_SIZE_WG                        SYNAPTIC_EVENT_BUFFERS
  #endif
  #define EXPAND_EVENTS_SPIKE_PACKETS_PER_WF                  (EXPAND_EVENTS_SPIKE_PACKETS/\
                                                              (EXPAND_EVENTS_GRID_SIZE_WG*\
                                                              EXPAND_EVENTS_WG_SIZE_WF))
  /*Spike data structure parameters (data in)*/
  #define EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE                SPIKE_PACKET_SIZE
  #define EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS            2
  #define EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS              (EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE *\
                                                              EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS)
  /*Spike buffer size per WF in local memory*/
  #if !(defined(EXPAND_EVENTS_SPIKE_BUFFER_SIZE))
    #define EXPAND_EVENTS_SPIKE_BUFFER_SIZE                   128
  #endif
  
  /*Time slot (bin) buffer parameters (data out)*/
  /*CONTROL: Max number of time slots (define the longest delay)*/
  #define EXPAND_EVENTS_TIME_SLOTS                            EVENT_TIME_SLOTS
  #define EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS                SYNAPTIC_EVENT_BUFFERS
  #define EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS            3
  #define EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE   EVENT_DATA_BUFFER_SIZE
  
  /*NN parameters*/
  /*CONTROL: Network size in terms of bits (1<<EXPAND_EVENTS_TOTAL_NEURON_BITS)*/
  #define EXPAND_EVENTS_TOTAL_NEURON_BITS                     TOTAL_NEURON_BITS
  #define EXPAND_EVENTS_TOTAL_NEURONS                         (1<<EXPAND_EVENTS_TOTAL_NEURON_BITS)
  #define EXPAND_EVENTS_MAX_DELAY                             (EXPAND_EVENTS_TIME_SLOTS-\
                                                              SIMULATION_STEP_SIZE)
  /*CONTROL: Minimum event delay latency*/
  #define EXPAND_EVENTS_MIN_DELAY                             (MINIMUM_PROPAGATION_DELAY)
  
  /*Histogram of target neurons*/
  #define EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM               1
  #define EXPAND_EVENTS_HISTOGRAM_BIT_SHIFT                   0
  #define EXPAND_EVENTS_HISTOGRAM_BIN_BITS                    4 /*Must be aligned with EXPAND_EVENTS_HISTOGRAM_BIN_MASK*/
  #define EXPAND_EVENTS_HISTOGRAM_BIN_MASK                    0xF /*Must be aligned with EXPAND_EVENTS_HISTOGRAM_BIN_BITS*/
  #define EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS                  (1<<EXPAND_EVENTS_HISTOGRAM_BIN_BITS)
  
/*
                                              Local memory layout
*/
  #if EXPAND_EVENTS_INTER_WF_COOPERATION == 0

  /*Local memory offsets*/
  #if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
    #define EXPAND_EVENTS_CACHE_OFFSET_1                      (EXPAND_EVENTS_TIME_SLOTS*\
                                                              EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS*\
                                                              EXPAND_EVENTS_WG_SIZE_WF)
  #else
    #define EXPAND_EVENTS_CACHE_OFFSET_1    0
  #endif
  #define EXPAND_EVENTS_CACHE_OFFSET_2                        (EXPAND_EVENTS_CACHE_OFFSET_1 + \
                                                              EXPAND_EVENTS_TIME_SLOTS*\
                                                              EXPAND_EVENTS_WG_SIZE_WF)
  #if EXPAND_EVENTS_SPIKE_BUFFER_SIZE > EXPAND_EVENTS_SPIKE_PACKETS_PER_WF
    #define EXPAND_EVENTS_SPIKE_COUNT_PITCH                   EXPAND_EVENTS_SPIKE_BUFFER_SIZE
  #else
    #define EXPAND_EVENTS_SPIKE_COUNT_PITCH                   EXPAND_EVENTS_SPIKE_PACKETS_PER_WF
  #endif
  #define EXPAND_EVENTS_CACHE_OFFSET_3                        (EXPAND_EVENTS_CACHE_OFFSET_2 + \
                                                              (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*\
                                                              EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + \
                                                              EXPAND_EVENTS_SPIKE_COUNT_PITCH)*\
                                                              EXPAND_EVENTS_WG_SIZE_WF)
  #define EXPAND_EVENTS_CACHE_SIZE_WORDS                      EXPAND_EVENTS_CACHE_OFFSET_3

  /*Local memory accessors*/
  #define HISTOGRAM(i)                                        cache[i]
  #define TIME_SLOT_COUNTERS(i)                               cache[EXPAND_EVENTS_CACHE_OFFSET_1 + i]
  #define SPIKE_DATA(i)                                       cache[EXPAND_EVENTS_CACHE_OFFSET_2 + i]
  
  #else /*EXPAND_EVENTS_INTER_WF_COOPERATION > 0*/

  /*Local memory offsets*/
  #if (EXPAND_EVENTS_ENABLE_TARGET_HISTOGRAM)
    #define EXPAND_EVENTS_CACHE_OFFSET_1                      (EXPAND_EVENTS_TIME_SLOTS*\
                                                              EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)
  #else
    #define EXPAND_EVENTS_CACHE_OFFSET_1    0
  #endif
  #define EXPAND_EVENTS_CACHE_OFFSET_2                        (EXPAND_EVENTS_CACHE_OFFSET_1 + \
                                                              EXPAND_EVENTS_TIME_SLOTS)
  #define EXPAND_EVENTS_CACHE_OFFSET_3                        (EXPAND_EVENTS_CACHE_OFFSET_2 + 1)
  #if EXPAND_EVENTS_SPIKE_BUFFER_SIZE > (EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*EXPAND_EVENTS_WG_SIZE_WF \
    + EXPAND_EVENTS_WG_SIZE_WF)
    #define EXPAND_EVENTS_SPIKE_COUNT_PITCH                   EXPAND_EVENTS_SPIKE_BUFFER_SIZE
  #else
    #define EXPAND_EVENTS_SPIKE_COUNT_PITCH                   (EXPAND_EVENTS_SPIKE_PACKETS_PER_WF*\
                                                              EXPAND_EVENTS_WG_SIZE_WF +\
                                                              EXPAND_EVENTS_WG_SIZE_WF)
  #endif
  #define EXPAND_EVENTS_CACHE_OFFSET_4                        (EXPAND_EVENTS_CACHE_OFFSET_3 + \
                                                              (EXPAND_EVENTS_SPIKE_BUFFER_SIZE*\
                                                              EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS + \
                                                              EXPAND_EVENTS_SPIKE_COUNT_PITCH))
  #define EXPAND_EVENTS_CACHE_SIZE_WORDS                      EXPAND_EVENTS_CACHE_OFFSET_4

  /*Local memory accessors*/
  #define HISTOGRAM(i)                                        cache[i]
  #define TIME_SLOT_COUNTERS(i)                               cache[EXPAND_EVENTS_CACHE_OFFSET_1 + i]
  #define TOTAL_SPIKES                                        cache[EXPAND_EVENTS_CACHE_OFFSET_2]
  #define SPIKE_DATA(i)                                       cache[EXPAND_EVENTS_CACHE_OFFSET_3 + i]
  #endif
/*
                                              Optimizations
*/
  /*CONTROL: Reduce loops to if statements where possible*/
  #define EXPAND_EVENTS_OPTIMIZATION_REDUCE_LOOPS             1
/*
                                              Restrictions
*/
#if EXPAND_EVENTS_INTER_WF_COOPERATION == 0
#if	((SYNAPTIC_EVENT_BUFFERS%EXPAND_EVENTS_WG_SIZE_WF) != 0)
  #error Parameter SYNAPTIC_EVENT_BUFFERS must be divisible by EXPAND_EVENTS_WG_SIZE_WF if\
  EXPAND_EVENTS_INTER_WF_COOPERATION == 0
#endif
#endif
#if	(EXPAND_EVENTS_SPIKE_PACKETS%(EXPAND_EVENTS_GRID_SIZE_WG*EXPAND_EVENTS_WG_SIZE_WF) != 0)
  #error Parameter EXPAND_EVENTS_SPIKE_PACKETS must be divisible by \
    (EXPAND_EVENTS_GRID_SIZE_WG*EXPAND_EVENTS_WG_SIZE_WF)
#endif
#endif /*EXPAND_EVENTS_ENABLE*/
/*
                                               Default Definitions
*/
#if !defined (EXPAND_EVENTS_DEBUG_ENABLE)
  #define EXPAND_EVENTS_DEBUG_ENABLE                          0
#endif
/* *** */
#if !defined (EXPAND_EVENTS_ERROR_TRACK_ENABLE)
  #define EXPAND_EVENTS_ERROR_TRACK_ENABLE                    0
#endif
/* *** */
#if !defined (EXPAND_EVENTS_VERIFY_ENABLE)
  #define EXPAND_EVENTS_VERIFY_ENABLE                         0
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
  /*CONTROL: enable verification*/
  #define SCAN_VERIFY_ENABLE                                  KERNEL_VERIFY_ENABLE
  
  /*Debugging*/
  /*CONTROL: enable debugging*/
  #define SCAN_DEBUG_ENABLE                                   (0)&DEBUG_MASK
  /*CONTROL: debug buffer size*/
  #define SCAN_DEBUG_BUFFER_SIZE_WORDS                        (1024)
  
  /*Error tracking and codes*/
  /*CONTROL: enable error tracking*/
  #define SCAN_ERROR_TRACK_ENABLE                             ERROR_TRACK_ENABLE
  #define SCAN_ERROR_BUFFER_SIZE_WORDS                        1
  #define SCAN_ERROR_CODE_1                                   0x1
  
  /*Scan kernel size parameters*/                               
  #define SCAN_KERNEL_FILE_NAME                               "Kernel_ScanHistogram.cl"
  #define SCAN_KERNEL_NAME                                    "scan_histogram"
  #define SCAN_WF_SIZE_WI                                     WF_SIZE_WI
  /*CONTROL: number of WFs in a WG*/
  #if !(defined(SCAN_WG_SIZE_WF))
    #define SCAN_WG_SIZE_WF                                   2 /*Options: 1, 2, 4*/
  #endif
  #define SCAN_WG_SIZE_WI                                     (SCAN_WG_SIZE_WF*SCAN_WF_SIZE_WI)
  #define SCAN_GRID_SIZE_WG                                   1 /*Options: 1*/
  
  /*Histogram of target neurons*/
  /*CONTROL: Max number of time slots (define the longest delay)*/
  #define SCAN_HISTOGRAM_TIME_SLOTS                           EVENT_TIME_SLOTS
  
  /*Synaptic event buffer parameters*/
  /*CONTROL: number of synaptic event buffers, which store synaptic events*/
  #define SCAN_SYNAPTIC_EVENT_BUFFERS                         SYNAPTIC_EVENT_BUFFERS
  /*CONTROL: size of event data buffer*/
  #define SCAN_EVENT_DATA_MAX_SRC_BUFFER_SIZE                 EVENT_DATA_BUFFER_SIZE
/*
                                              Optimizations
*/
  /*CONTROL: Reduction/prefix sum algorithm options*/
  #if !(defined(SCAN_OPTIMIZATION_2LEVEL_REDUCE))
    #define SCAN_OPTIMIZATION_2LEVEL_REDUCE                   1
  #endif
  /*CONTROL: Local memory offset. 0-17*/
  #if !(defined(SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET))
    #define SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET         0
  #endif
/*
                                              Local memory layout
*/
  #if !(defined(SCAN_OPTIMIZATION_INTER_WF_CACHE_OFFSET))
    #define SCAN_OPTIMIZATION_INTER_WF_CACHE_OFFSET           64
  #endif
  #define SCAN_WF_CACHE_SIZE_WORDS                            (192 - \
                                                              SCAN_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET + \
                                                              SCAN_OPTIMIZATION_INTER_WF_CACHE_OFFSET)
  #define SCAN_CACHE_SIZE_WORDS                               SCAN_WF_CACHE_SIZE_WORDS*\
                                                              SCAN_WG_SIZE_WF
/*
                                               Variants
*/
/*VARIANT_00*/
#if !(defined(SCAN_DEVICE_V00))
  #define SCAN_DEVICE_V00                                     0
#endif
#if SCAN_ENABLE_V00
  /*Host and Device*/
  #define SCAN_HISTOGRAM_TOTAL_BINS_V00                       (1<<4)
  #define SCAN_HISTOGRAM_BIN_SIZE_V00                         SYNAPTIC_EVENT_BUFFERS
  #define SCAN_HISTOGRAM_IN_TYPE_V00                          0
  #define SCAN_HISTOGRAM_MAX_COUNT_FOR_TEST_V00               (0xFFFFFFFF/\
                                                              (SCAN_HISTOGRAM_TOTAL_BINS_V00*\
                                                              SCAN_HISTOGRAM_BIN_SIZE_V00))
  /*Device*/
#if SCAN_DEVICE_V00
  #define SCAN_HISTOGRAM_TOTAL_BINS                           SCAN_HISTOGRAM_TOTAL_BINS_V00
  #define SCAN_HISTOGRAM_BIN_SIZE                             SCAN_HISTOGRAM_BIN_SIZE_V00
  #define SCAN_HISTOGRAM_IN_TYPE                              SCAN_HISTOGRAM_IN_TYPE_V00
  #define SCAN_HISTOGRAM_ELEMENTS_PER_WI                      ((SCAN_HISTOGRAM_TOTAL_BINS*\
                                                              SCAN_HISTOGRAM_BIN_SIZE)/\
                                                              SCAN_WG_SIZE_WI)
/*
                                               Restrictions
*/
#if	((SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE)%SCAN_WG_SIZE_WI != 0)
  #error (Parameter (SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE) must be divisible by \
  SCAN_WG_SIZE_WI)
#endif
/* *** */
#if	(SCAN_HISTOGRAM_ELEMENTS_PER_WI%4 != 0)
  #error (Parameter SCAN_HISTOGRAM_ELEMENTS_PER_WI must be divisible by 4)
#endif
/* *** */
#endif  /*END SCAN_DEVICE_V01*/
#endif  /*END VARIANT_00*/

/*VARIANT_01*/
#if !(defined(SCAN_DEVICE_V01))
  #define SCAN_DEVICE_V01                                     0
#endif
#if SCAN_ENABLE_V01
  /*Host and Device*/
  #define SCAN_HISTOGRAM_TOTAL_BINS_V01                       (1<<4)
  #define SCAN_HISTOGRAM_BIN_SIZE_V01                         SYNAPTIC_EVENT_BUFFERS
  #define SCAN_HISTOGRAM_IN_TYPE_V01                          1
  #define SCAN_HISTOGRAM_BIN_BACKETS                          SYNAPTIC_EVENT_BUFFERS
  #define SCAN_HISTOGRAM_MAX_COUNT_FOR_TEST_V01               (0xFFFFFFFF/\
                                                              (SCAN_HISTOGRAM_TOTAL_BINS_V01*\
                                                              SCAN_HISTOGRAM_BIN_SIZE_V01*\
                                                              SCAN_HISTOGRAM_BIN_BACKETS))
  /*Device*/
#if SCAN_DEVICE_V01
  #define SCAN_HISTOGRAM_TOTAL_BINS                           SCAN_HISTOGRAM_TOTAL_BINS_V01
  #define SCAN_HISTOGRAM_BIN_SIZE                             SCAN_HISTOGRAM_BIN_SIZE_V01
  #define SCAN_HISTOGRAM_IN_TYPE                              SCAN_HISTOGRAM_IN_TYPE_V01
  #define SCAN_HISTOGRAM_ELEMENTS_PER_WI                      ((SCAN_HISTOGRAM_TOTAL_BINS*\
                                                              SCAN_HISTOGRAM_BIN_SIZE)/\
                                                              SCAN_WG_SIZE_WI)
/*
                                               Restrictions
*/
#if	((SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE)%SCAN_WG_SIZE_WI != 0)
  #error (Parameter (SCAN_HISTOGRAM_TOTAL_BINS*SCAN_HISTOGRAM_BIN_SIZE) must be divisible by \
  SCAN_WG_SIZE_WI)
#endif
/* *** */
#if	(SCAN_HISTOGRAM_ELEMENTS_PER_WI%4 != 0)
  #error (Parameter SCAN_HISTOGRAM_ELEMENTS_PER_WI must be divisible by 4)
#endif
/* *** */
#endif  /*END SCAN_DEVICE_V01*/
#endif  /*END VARIANT_01*/
/*
                                               Restrictions
*/
/*Device sections of variant definitions can be enabled only one at a time*/
#if	(SCAN_DEVICE_V00 && SCAN_DEVICE_V01)
  #error (Parameters SCAN_DEVICE_V00 && SCAN_DEVICE_V01 must be enabled one at a time)
#endif
/* *** */
#endif /*SCAN_ENABLE*/
/*
                                               Default Definitions
*/
#if !defined (SCAN_DEBUG_ENABLE)
  #define SCAN_DEBUG_ENABLE                                   0
#endif
/* *** */
#if !defined (SCAN_ERROR_TRACK_ENABLE)
  #define SCAN_ERROR_TRACK_ENABLE                             0
#endif
/* *** */
#if !defined (SCAN_VERIFY_ENABLE)
  #define SCAN_VERIFY_ENABLE                                  0
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
  /*CONTROL: enable verification*/
  #define GROUP_EVENTS_VERIFY_ENABLE                            KERNEL_VERIFY_ENABLE
  
  /*Debugging*/
  /*CONTROL: enable debugging*/
  #define GROUP_EVENTS_DEBUG_ENABLE                             (0)&DEBUG_MASK
  /*CONTROL: debug buffer size*/
  #define GROUP_EVENTS_DEBUG_BUFFER_SIZE_WORDS                  (1024)
  
  /*Error tracking and codes*/
  /*CONTROL: enable error tracking*/
  #define GROUP_EVENTS_ERROR_TRACK_ENABLE                       ERROR_TRACK_ENABLE
  #define GROUP_EVENTS_ERROR_BUFFER_SIZE_WORDS                  1
  #define GROUP_EVENTS_ERROR_CODE_1                             0x1
  
  /*Testing*/
  /*CONTROL: test mode*/
  #if !(defined(GROUP_EVENTS_TEST_MODE))
    #define GROUP_EVENTS_TEST_MODE                             10
  #endif
  
  /*Synaptic event buffer parameters*/
  /*CONTROL: number of synaptic event buffers, which store synaptic events*/
  #define GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS                   SYNAPTIC_EVENT_BUFFERS
  /*CONTROL: size of event data buffer*/
  #define GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE           EVENT_DATA_BUFFER_SIZE
  #define GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE           (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS*\
                                                                GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE)
  /*CONTROL: how many buffers a single WG is assigned to process*/
  #if !(defined(GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG))
    #define GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG          1
  #endif
  #define GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS               3
  
  /*Kernel size parameters*/
  #define GROUP_EVENTS_KERNEL_FILE_NAME                         "Kernel_GroupEvents.cl"
  #define GROUP_EVENTS_KERNEL_NAME                              "group_events"
  #define GROUP_EVENTS_WF_SIZE_WI                               WF_SIZE_WI
  /*CONTROL: Each WI takes this number of data elements*/
  #if !(defined(GROUP_EVENTS_ELEMENTS_PER_WI))
    #define GROUP_EVENTS_ELEMENTS_PER_WI                        4
  #endif
  /*CONTROL: number of WFs in a WG*/
  #if !(defined(GROUP_EVENTS_WG_SIZE_WF))
    #define GROUP_EVENTS_WG_SIZE_WF                             4
  #endif
  #define GROUP_EVENTS_WG_SIZE_WI                               (GROUP_EVENTS_WG_SIZE_WF*\
                                                                GROUP_EVENTS_WF_SIZE_WI)
  #define GROUP_EVENTS_GRID_SIZE_WG                             (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS/\
                                                                GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS_PER_WG)
  
  /*Incoming histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_BITS                       4   /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_MASK*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_MASK                       0xF /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_BITS*/
  #define GROUP_EVENTS_HISTOGRAM_TOTAL_BINS                     (1<<GROUP_EVENTS_HISTOGRAM_BIN_BITS)
  #define GROUP_EVENTS_HISTOGRAM_BIN_SIZE                       (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS)

  /*Outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT                   4   /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT*/
  #define GROUP_EVENTS_HISTOGRAM_BIN_MASK_OUT                   0xF /*Must be aligned with GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT*/
  #define GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT                 (1<<GROUP_EVENTS_HISTOGRAM_BIN_BITS_OUT)
  #define GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE                  (GROUP_EVENTS_GRID_SIZE_WG)

  /*NN parameters*/
  #define GROUP_EVENTS_TOTAL_NEURON_BITS                        TOTAL_NEURON_BITS /*16*/
  #define GROUP_EVENTS_TOTAL_NEURONS                            (1<<GROUP_EVENTS_TOTAL_NEURON_BITS)
  #define GROUP_EVENTS_MAX_DELAY                                (GROUP_EVENTS_TIME_SLOTS-SIMULATION_STEP_SIZE)
  #define GROUP_EVENTS_MIN_DELAY                                1.0f
  
  /*CONTROL: Time slots for storing events*/
  #define GROUP_EVENTS_TIME_SLOTS                               EVENT_TIME_SLOTS /*Options: 16*/
  
  /*Check data boundary condition before fetching*/
  #define GROUP_EVENTS_CHECK_BOUNDARY                           1
  
  /*CONTROL: Number of threads sharing a counter for computing local histogram*/
  #if !(defined(GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER))
    #define GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER             16
  #endif
  
  #define GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS            ((GROUP_EVENTS_WF_SIZE_WI/\
                                                                GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER)*\
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)
  /*Number of sort iterations for neuron ID depends on the max ID*/
#define GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS \
                              ((GROUP_EVENTS_TOTAL_NEURON_BITS/GROUP_EVENTS_HISTOGRAM_BIN_BITS) +\
                              ((GROUP_EVENTS_TOTAL_NEURON_BITS%GROUP_EVENTS_HISTOGRAM_BIN_BITS) > 0))
/*
                                              Optimizations
*/
  /*CONTROL: Reduce loops to if statements where possible*/
  #if !(defined(GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS))
    #define GROUP_EVENTS_OPTIMIZATION_REDUCE_LOOPS              1
  #endif
  /*CONTROL: Reduction/prefix sum algorithm options*/
  #if !(defined(GROUP_EVENTS_OPTIMIZATION_2LEVEL_REDUCE))
    #define GROUP_EVENTS_OPTIMIZATION_2LEVEL_REDUCE             1
  #endif
  /*CONTROL: Enable sort instead of atomics if 1. 
  Warning: atomics do not preserve the sort order for subsequent passes. Useful only for the
  first pass*/
  #if !(defined(GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT))
    #define GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT                1
  #endif
  /*CONTROL: keep target histogram in GM/L1,L2 instead of LDS*/
  #if !(defined(GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT))
    #define GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT  1
  #endif
  /*CONTROL: prefix sume shared memory offset to control bank conflicts.*/
  #if !(defined(GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET))
    #define GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET   17
  #endif
/*
                                              Local memory layout
*/
  #define GROUP_EVENTS_CACHE_PREFIX_SUM_WF_SIZE           (192 - \
                                                          GROUP_EVENTS_OPTIMIZATION_CACHE_PREFIX_SUM_OFFSET)
#define GROUP_EVENTS_CACHE_OFFSET_1                       (GROUP_EVENTS_CACHE_PREFIX_SUM_WF_SIZE)
#if GROUP_EVENTS_CACHE_OFFSET_1 < (GROUP_EVENTS_WF_SIZE_WI*GROUP_EVENTS_ELEMENTS_PER_WI)
  #undef GROUP_EVENTS_CACHE_OFFSET_1
  #define GROUP_EVENTS_CACHE_OFFSET_1                     (GROUP_EVENTS_WF_SIZE_WI*\
                                                          GROUP_EVENTS_ELEMENTS_PER_WI)
#endif
#if !GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT
  #undef GROUP_EVENTS_CACHE_OFFSET_1
  #define GROUP_EVENTS_CACHE_OFFSET_1                     0
#endif
#define GROUP_EVENTS_CACHE_OFFSET_2                       (GROUP_EVENTS_LOCAL_HISTOGRAM_WF_SIZE_WORDS)
#define GROUP_EVENTS_CACHE_OFFSET_3                       (2*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS)

  /*Define minimum per-WF cache size required for kernel*/
#define GROUP_EVENTS_WF_CACHE_SIZE_WORDS                  GROUP_EVENTS_CACHE_OFFSET_1
#if GROUP_EVENTS_WF_CACHE_SIZE_WORDS < GROUP_EVENTS_CACHE_OFFSET_2
  #undef GROUP_EVENTS_WF_CACHE_SIZE_WORDS
  #define GROUP_EVENTS_WF_CACHE_SIZE_WORDS                GROUP_EVENTS_CACHE_OFFSET_2
#endif
#if GROUP_EVENTS_WF_CACHE_SIZE_WORDS < GROUP_EVENTS_CACHE_OFFSET_3
  #undef GROUP_EVENTS_WF_CACHE_SIZE_WORDS
  #define GROUP_EVENTS_WF_CACHE_SIZE_WORDS                GROUP_EVENTS_CACHE_OFFSET_3
#endif

  /*Histogram offset has to be on the top of others*/
#define GROUP_EVENTS_CACHE_OFFSET_4                       ((GROUP_EVENTS_WF_CACHE_SIZE_WORDS*\
                                                          GROUP_EVENTS_WG_SIZE_WF) + \
                                                          GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*\
                                                          (GROUP_EVENTS_WG_SIZE_WF-1))
                                                          
  /*Define minimum total and per-WF cache size required for kernel*/
#define GROUP_EVENTS_CACHE_SIZE_WORDS                     GROUP_EVENTS_WF_CACHE_SIZE_WORDS*\
                                                          GROUP_EVENTS_WG_SIZE_WF
#if (GROUP_EVENTS_WG_SIZE_WF > 1)
#if GROUP_EVENTS_CACHE_SIZE_WORDS < GROUP_EVENTS_CACHE_OFFSET_4
  #undef GROUP_EVENTS_CACHE_SIZE_WORDS
  #define GROUP_EVENTS_CACHE_SIZE_WORDS                   GROUP_EVENTS_CACHE_OFFSET_4
#endif
#endif
/*
                                               Variants
*/
/*
            VARIANT_00
            Loads key data from original data structure. Computes values as pointers to the original
            data. Sorts by keys the key-value pairs based on provided scanned histogram. 
            Computes new histogram for next stage.
*/
#if !(defined(GROUP_EVENTS_DEVICE_V00))
  #define GROUP_EVENTS_DEVICE_V00                               0
#endif
#if GROUP_EVENTS_ENABLE_V00
  /*Host and Device*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00                  0
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V00              4
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V00                    0
  #define GROUP_EVENTS_VALUES_MODE_V00                          1
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET_V00                    1
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V00    0
  #define GROUP_EVENTS_RELOCATE_VALUES_V00                      0
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V00          1
  #define GROUP_EVENTS_CACHE_SIZE_WORDS_V00                     (GROUP_EVENTS_CACHE_SIZE_WORDS + \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + \
                                                                GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V00 * \
                                                                GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT * \
                                                                (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE * \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))

  /*Device*/
#if GROUP_EVENTS_DEVICE_V00
  /*Enable outgoing histogram of target neurons*/
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT              GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V00
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V00
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V00
  /*Source events may be delivered in different data structures*/
  /*Block-partitioned by exapnd kernel*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V00 
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        GROUP_EVENTS_ENABLE_STEP_SHIFT_V00
  /*Each key gets a value of its address in GM*/
  #define GROUP_EVENTS_VALUES_MODE                              GROUP_EVENTS_VALUES_MODE_V00
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        GROUP_EVENTS_SOURCE_KEY_OFFSET_V00
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          GROUP_EVENTS_RELOCATE_VALUES_V00
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              0
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
#endif  /*END VARIANT_00*/
/*
            VARIANT_01
            Sorts by keys the key-value pairs based on provided scanned histogram. 
            Computes new histogram for next stage.
*/
#if !(defined(GROUP_EVENTS_DEVICE_V01))
  #define GROUP_EVENTS_DEVICE_V01                               0
#endif
#if GROUP_EVENTS_ENABLE_V01
  /*Host-visible copies of device parameters*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01                  4
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01              4
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V01                    1
  #define GROUP_EVENTS_VALUES_MODE_V01                          2
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V01    1
  #define GROUP_EVENTS_RELOCATE_VALUES_V01                      0
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V01          1
  #define GROUP_EVENTS_CACHE_SIZE_WORDS_V01                     (GROUP_EVENTS_CACHE_SIZE_WORDS + \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + \
                                                                GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V01 * \
                                                                GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT * \
                                                                (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE * \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))

  /*Device*/
#if GROUP_EVENTS_DEVICE_V01
  /*Enable outgoing histogram of target neurons*/
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT              GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V01
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V01
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V01
  /*Source events may be delivered in different data structures*/
  /*Contiguous*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V01
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        GROUP_EVENTS_ENABLE_STEP_SHIFT_V01
  /*Each key gets a value carried out from before*/
  #define GROUP_EVENTS_VALUES_MODE                              GROUP_EVENTS_VALUES_MODE_V01
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        0
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          GROUP_EVENTS_RELOCATE_VALUES_V01
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              0
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
#endif/*END VARIANT_01*/
/*
            VARIANT_02
            This variant is the same as VARIANT_01 except it replaces data with new key at the end.
            This allows to continue sorting with new key in the next stages.
*/
#if !(defined(GROUP_EVENTS_DEVICE_V02))
  #define GROUP_EVENTS_DEVICE_V02                               0
#endif
#if GROUP_EVENTS_ENABLE_V02
  /*Host-visible copies of device parameters*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02                  4
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V02              0
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V02                    1
  #define GROUP_EVENTS_VALUES_MODE_V02                          2
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET_V02               0
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V02    1
  #define GROUP_EVENTS_RELOCATE_VALUES_V02                      0
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V02          1
  #define GROUP_EVENTS_CACHE_SIZE_WORDS_V02                     (GROUP_EVENTS_CACHE_SIZE_WORDS + \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + \
                                                                GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V02 * \
                                                                GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT * \
                                                                (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE * \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))

  /*Device*/
#if GROUP_EVENTS_DEVICE_V02
  /*Enable outgoing histogram of target neurons*/
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT              GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V02
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V02
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V02
  /*Source events may be delivered in different data structures*/
  /*Contiguous*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V02
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        GROUP_EVENTS_ENABLE_STEP_SHIFT_V02
  /*Each key gets a value carried out from before*/
  #define GROUP_EVENTS_VALUES_MODE                              GROUP_EVENTS_VALUES_MODE_V02
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        0
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          GROUP_EVENTS_RELOCATE_VALUES_V02
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              1
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   GROUP_EVENTS_REPLACEMENT_KEY_OFFSET_V02
#endif
#endif/*END VARIANT_02*/
/*
            VARIANT_03
            This variant is the same as VARIANT_01 except it relocate original values in place
            of pointers. No histogram for the next stage is computed.
*/
#if !(defined(GROUP_EVENTS_DEVICE_V03))
  #define GROUP_EVENTS_DEVICE_V03                               0
#endif
/*TODO: need to disable histogram out*/
#if GROUP_EVENTS_ENABLE_V03
  /*Host-visible copies of device parameters*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03                  4
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V03              0
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT_V03                    1
  #define GROUP_EVENTS_VALUES_MODE_V03                          2
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V03    1
  #define GROUP_EVENTS_RELOCATE_VALUES_V03                      1
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V03          0
  #define GROUP_EVENTS_CACHE_SIZE_WORDS_V03                     (GROUP_EVENTS_CACHE_SIZE_WORDS + \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS + \
                                                                GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V03 * \
                                                                GROUP_EVENTS_OPTIMIZATION_LDS_TARGET_HISTOGRAM_OUT * \
                                                                (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE * \
                                                                GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT))

  /*Device*/
#if GROUP_EVENTS_DEVICE_V03
  /*Enable outgoing histogram of target neurons*/
  #define GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT              GROUP_EVENTS_ENABLE_TARGET_HISTOGRAM_OUT_V03
  /*Bit shift used for computing the incoming histogram of target neurons passed as a parameter*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT                      GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_V03
  /*Bit shift used for computing the outgoing histogram of target neurons*/
  #define GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT                  GROUP_EVENTS_HISTOGRAM_BIT_SHIFT_OUT_V03
  /*Source events may be delivered in different data structures*/
  /*Contiguous*/
  #define GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE        GROUP_EVENTS_SOURCE_EVENTS_DATA_STRUCTURE_TYPE_V03
  /*Linking step parameter to shift parameter used in masking bits for sort*/
  #define GROUP_EVENTS_ENABLE_STEP_SHIFT                        GROUP_EVENTS_ENABLE_STEP_SHIFT_V03
  /*Each key gets a value carried out from before*/
  #define GROUP_EVENTS_VALUES_MODE                              GROUP_EVENTS_VALUES_MODE_V03
  /*Key offset: which element to use as a key from GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS*/
  #define GROUP_EVENTS_SOURCE_KEY_OFFSET                        0
  /*Relocate original values based on pointers stored in their place (useful for multiple values/key)*/
  #define GROUP_EVENTS_RELOCATE_VALUES                          GROUP_EVENTS_RELOCATE_VALUES_V03
  /*Replace key. Useful if need to sort with a new key.*/
  #define GROUP_EVENTS_REPLACE_KEY                              0
  #define GROUP_EVENTS_REPLACEMENT_KEY_OFFSET                   0
#endif
#endif/*END VARIANT_03*/
/*
                                               Restrictions
*/
/*Device sections of variant definitions can be enabled only one at a time*/
#if	(GROUP_EVENTS_DEVICE_V00 && GROUP_EVENTS_DEVICE_V01)
  #error \
    (Parameters GROUP_EVENTS_DEVICE_V00 && GROUP_EVENTS_DEVICE_V01 must be enabled one at a time)
#endif
#if	(GROUP_EVENTS_DEVICE_V00 && GROUP_EVENTS_DEVICE_V02)
  #error \
    (Parameters GROUP_EVENTS_DEVICE_V00 && GROUP_EVENTS_DEVICE_V02 must be enabled one at a time)
#endif
#if	(GROUP_EVENTS_DEVICE_V00 && GROUP_EVENTS_DEVICE_V03)
  #error \
    (Parameters GROUP_EVENTS_DEVICE_V00 && GROUP_EVENTS_DEVICE_V03 must be enabled one at a time)
#endif
#if	(GROUP_EVENTS_DEVICE_V01 && GROUP_EVENTS_DEVICE_V02)
  #error \
    (Parameters GROUP_EVENTS_DEVICE_V01 && GROUP_EVENTS_DEVICE_V02 must be enabled one at a time)
#endif
#if	(GROUP_EVENTS_DEVICE_V01 && GROUP_EVENTS_DEVICE_V03)
  #error \
    (Parameters GROUP_EVENTS_DEVICE_V01 && GROUP_EVENTS_DEVICE_V03 must be enabled one at a time)
#endif
#if	(GROUP_EVENTS_DEVICE_V02 && GROUP_EVENTS_DEVICE_V03)
  #error \
    (Parameters GROUP_EVENTS_DEVICE_V02 && GROUP_EVENTS_DEVICE_V03 must be enabled one at a time)
#endif
/* *** */
#if GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT && (GROUP_EVENTS_ELEMENTS_PER_WI != 4)
  #error GROUP_EVENTS_ELEMENTS_PER_WI has to be 4 if GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT is enabled
#endif
/* *** */
#if GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*2 > GROUP_EVENTS_WF_SIZE_WI
  #error GROUP_EVENTS_HISTOGRAM_TOTAL_BINS*2 cannot be more than GROUP_EVENTS_WF_SIZE_WI
#endif
/* *** */
#if GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT && (GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER != 16)
  #error BUG: GROUP_EVENTS_WIs_PER_BIN_COUNTER_BUFFER has to be 16 if \
    GROUP_EVENTS_OPTIMIZATION_LOCAL_SORT is set
#endif
/* *** */
#endif /*GROUP_EVENTS_ENABLE*/
/*
                                               Default Definitions
*/
#if !defined (GROUP_EVENTS_DEBUG_ENABLE)
  #define GROUP_EVENTS_DEBUG_ENABLE                                   0
#endif
/* *** */
#if !defined (GROUP_EVENTS_ERROR_TRACK_ENABLE)
  #define GROUP_EVENTS_ERROR_TRACK_ENABLE                             0
#endif
/* *** */
#if !defined (GROUP_EVENTS_VERIFY_ENABLE)
  #define GROUP_EVENTS_VERIFY_ENABLE                                  0
#endif
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
  #define MAKE_EVENT_PTRS_DEBUG_BUFFER_SIZE_WORDS                  (1024)
  
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
  #if !(defined(MAKE_EVENT_PTRS_WG_SIZE_WF))
    #define MAKE_EVENT_PTRS_WG_SIZE_WF                             4
  #endif
  #define MAKE_EVENT_PTRS_WG_SIZE_WI                                (MAKE_EVENT_PTRS_WG_SIZE_WF*\
                                                                    MAKE_EVENT_PTRS_WF_SIZE_WI)
  /*CONTROL: each WI processes this number of data elements*/
  #if !(defined(MAKE_EVENT_PTRS_ELEMENTS_PER_WI))
    #define MAKE_EVENT_PTRS_ELEMENTS_PER_WI                        8
  #endif
  
  /*NN parameters*/
  /*CONTROL: total neurons in the network*/
  #define MAKE_EVENT_PTRS_TOTAL_NEURON_BITS                         TOTAL_NEURON_BITS
  #define MAKE_EVENT_PTRS_TOTAL_NEURONS                             (1<<(MAKE_EVENT_PTRS_TOTAL_NEURON_BITS))
  
  /*Buffer and memory sizes*/
  #define MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE                      2
  /*Synaptic event input buffer parameters*/
  #define MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS                   3
  /*LM allocation for scan per WF. 
  TODO: this can be reduced by 17 at least*/
  #define MAKE_EVENT_PTRS_SCAN_WF_LM_SHARE                         192
  /*Offset in GM for location of total events*/
  #define MAKE_EVENT_PTRS_HISTOGRAM_BIN_BITS                       4   
  #define MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET                      (SYNAPTIC_EVENT_BUFFERS*(1<<4))
  
  /*Event delivery mode:
    0 - array size of neuron count
    1 - structs*/
  #define MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE                       0
  
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  /*CONTROL: WG count*/
  #if !(defined(MAKE_EVENT_PTRS_GRID_SIZE_WG))
    #define MAKE_EVENT_PTRS_GRID_SIZE_WG                            128
  #endif
  
  /*Output buffer size*/
  #define MAKE_EVENT_PTRS_STRUCTS                                   1
  #define MAKE_EVENT_PTRS_STRUCT_SIZE                               MAKE_EVENT_PTRS_TOTAL_NEURONS*\
                                                                    MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE+ \
                                                                    (2*MAKE_EVENT_PTRS_GRID_SIZE_WG*\
                                                                    MAKE_EVENT_PTRS_WG_SIZE_WF*\
                                                                    MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE)
  /*Unit test related parameters*/
  /*CONTROL: Source data size for unit test. Can be viewed as elements per neuron*/
  #if !(defined(MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE))
  #define MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE                  (SYNAPTIC_EVENT_BUFFERS*\
                                                                    EVENT_DATA_BUFFER_SIZE)
  #endif
  /*CONTROL: Max neuron ID for unit test.*/
  #if !(defined(MAKE_EVENT_PTRS_TEST_MAX_NEURON_ID))
    #define MAKE_EVENT_PTRS_TEST_MAX_NEURON_ID                     (MAKE_EVENT_PTRS_TOTAL_NEURONS-1)
  #endif
  
  /*CONTROL: LM allocation for computing event counts per WF*/
  #if !(defined(MAKE_EVENT_PTRS_EVENT_COUNT_WF_LM_SHARE))
    #define MAKE_EVENT_PTRS_EVENT_COUNT_WF_LM_SHARE                 1024
  #endif
  
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
/*
                                               Default Definitions
*/
#if !defined (MAKE_EVENT_PTRS_DEBUG_ENABLE)
  #define MAKE_EVENT_PTRS_DEBUG_ENABLE                                   0
#endif
/* *** */
#if !defined (MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE)
  #define MAKE_EVENT_PTRS_ERROR_TRACK_ENABLE                             0
#endif
/* *** */
#if !defined (MAKE_EVENT_PTRS_VERIFY_ENABLE)
  #define MAKE_EVENT_PTRS_VERIFY_ENABLE                                  0
#endif
/**************************************************************************************************/



/***************************************************************************************************

  Kernel configuration interface for updating model variables
  
***************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
/*
                                               Generic parameters
*/
  /*Verification*/
  /*CONTROL: enable verification*/
  #define UPDATE_NEURONS_VERIFY_ENABLE                            KERNEL_VERIFY_ENABLE
  
  /*Error tracking and codes*/
  /*CONTROL: enable error logging at kernel level*/
  #define UPDATE_NEURONS_ERROR_TRACK_ENABLE                       ERROR_TRACK_ENABLE
  #define UPDATE_NEURONS_ERROR_BUFFER_SIZE_WORDS                  1
  #define UPDATE_NEURONS_ERROR_NON_SOLVER_FAILURE_MASK            0x000000D2
  #define UPDATE_NEURONS_ERROR_CODE_1                             1   /*UPDATE_NEURONS_PS_ORDER_LIMIT was hit*/
  #define UPDATE_NEURONS_ERROR_CODE_2                             2   /*A neuron spiked more than once*/
  #define UPDATE_NEURONS_ERROR_CODE_3                             4   /*UPDATE_NEURONS_NR_ORDER_LIMIT was hit*/
  #define UPDATE_NEURONS_ERROR_CODE_4                             8   /*NR diverged*/
  #define UPDATE_NEURONS_ERROR_CODE_5                             16  /*Spike packet size limit was hit*/
  #define UPDATE_NEURONS_ERROR_CODE_6                             32  /*PS divergence*/
  #define UPDATE_NEURONS_ERROR_CODE_7                             64  /*Spike count in NN overflow*/
  #define UPDATE_NEURONS_ERROR_CODE_8                             128 /*Event count overflow*/
  
  /*Kernel file and name*/
  #define UPDATE_NEURONS_KERNEL_FILE_NAME                         "Kernel_UpdateNeurons.cl"
  #define UPDATE_NEURONS_KERNEL_NAME                              "update_neurons"
  #define UPDATE_NEURONS_SPIKED_KERNEL_NAME                       "update_spiked_neurons"
  /*WF size measured in WIs*/
  #define UPDATE_NEURONS_WF_SIZE_WI                               WF_SIZE_WI /*Options: 64*/
  
  /*Size of element in event pointer struct*/
  #define UPDATE_NEURONS_STRUCT_ELEMENT_SIZE                      2
  /*Event input buffer element size*/
  #define UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS                   3
  /*CONTROL: each WI processes this number of neurons*/
  #if !(defined(UPDATE_NEURONS_ELEMENTS_PER_WI))
    #define UPDATE_NEURONS_ELEMENTS_PER_WI                        1
  #endif
  /*Size of time delay buffer*/
  #define UPDATE_NEURONS_TIME_SLOTS                               EVENT_TIME_SLOTS /*Options: 16*/
  
  /*Sizes of model variables and parameters*/
  #define UPDATE_NEURONS_MODEL_VARIABLES                          4
  #define UPDATE_NEURONS_MODEL_PARAMETERS                         9
  
  /*Spike data structure parameters*/
  #define UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE                   SPIKE_PACKET_SIZE
  #define UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS               2
  #define UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS                  (UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE * \
                                                                  UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS)
  /*NN parameters*/
  #define UPDATE_NEURONS_TOTAL_NEURON_BITS                        TOTAL_NEURON_BITS /*16*/
  #define UPDATE_NEURONS_TOTAL_NEURONS                            (1<<UPDATE_NEURONS_TOTAL_NEURON_BITS)
  #define UPDATE_NEURONS_MAX_DELAY                                (UPDATE_NEURONS_TIME_SLOTS-SIMULATION_STEP_SIZE)
  #define UPDATE_NEURONS_MIN_DELAY                                1.0f
  
  /*Debugging*/
  /*CONTROL: enable debugging*/
  #define UPDATE_NEURONS_DEBUG_ENABLE                             (0)&DEBUG_MASK
  /*CONTROL: debug buffer size*/
  #define UPDATE_NEURONS_DEBUG_BUFFER_SIZE_WORDS                  (1024)
  
  /*CONTROL: Source data size for unit test*/
  #if !(defined(UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE))
  #define UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE                 (SYNAPTIC_EVENT_BUFFERS*\
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
  
  /*Number of constant coefficient types*/
  #define UPDATE_NEURONS_CONST_COEFFICIENTS                       4
  /*Accessor to constant coefficients*/
  #define CONST_CO(mem, type, index)                              *(mem + type*\
                                                                  (UPDATE_NEURONS_PS_ORDER_LIMIT+1)+\
                                                                  index)
  #define CONST_TOL(mem)                                          *(mem + UPDATE_NEURONS_CONST_COEFFICIENTS*\
                                                                  (UPDATE_NEURONS_PS_ORDER_LIMIT+1))
  #define CONST_SIZE                                              (UPDATE_NEURONS_CONST_COEFFICIENTS*\
                                                                  (UPDATE_NEURONS_PS_ORDER_LIMIT+1) + 1)
  
  /*Optimizations*/
  #define UPDATE_NEURONS_DT_1_0_OPTIMIZATION                      1
/*
                                               Variants
*/
/*
            VARIANT_00

*/
#if UPDATE_NEURONS_ENABLE_V00
  /*Wisible to Host and Device*/
  /*CONTROL: number of WFs per WG*/
  #if !(defined(UPDATE_NEURONS_WG_SIZE_WF_V00))
  #define UPDATE_NEURONS_WG_SIZE_WF_V00                           1
  #endif
  /*CONTROL: number of spike packets generated by kernel*/
  #define UPDATE_NEURONS_SPIKE_PACKETS_V00                        SPIKE_PACKETS
  #define UPDATE_NEURONS_WG_SIZE_WI_V00                           (UPDATE_NEURONS_WG_SIZE_WF_V00*\
                                                                  UPDATE_NEURONS_WF_SIZE_WI)
  #define UPDATE_NEURONS_GRID_SIZE_WG_V00                         (UPDATE_NEURONS_SPIKE_PACKETS_V00/\
                                                                  UPDATE_NEURONS_WG_SIZE_WF_V00)
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

  /*Visible to Device*/
#ifdef UPDATE_NEURONS_DEVICE_V00
  #define UPDATE_NEURONS_WG_SIZE_WF                               UPDATE_NEURONS_WG_SIZE_WF_V00
  #define UPDATE_NEURONS_WG_SIZE_WI                               UPDATE_NEURONS_WG_SIZE_WI_V00
  #define UPDATE_NEURONS_GRID_SIZE_WG                             UPDATE_NEURONS_GRID_SIZE_WG_V00
#endif
/*
                                      Restrictions, VARIANT_00
*/
#if(UPDATE_NEURONS_TOTAL_NEURONS%(UPDATE_NEURONS_ELEMENTS_PER_WI*UPDATE_NEURONS_WG_SIZE_WI_V00*\
   UPDATE_NEURONS_GRID_SIZE_WG_V00))
  #error(UPDATE_NEURONS_TOTAL_NEURONS must be multiple of (UPDATE_NEURONS_ELEMENTS_PER_WI*\
    UPDATE_NEURONS_WG_SIZE_WI_V00*UPDATE_NEURONS_GRID_SIZE_WG_V00))
#endif
/* *** */
#if(UPDATE_NEURONS_SPIKE_PACKETS_V00%UPDATE_NEURONS_WG_SIZE_WF_V00)
  #error(UPDATE_NEURONS_SPIKE_PACKETS_V00 must be multiple of UPDATE_NEURONS_WG_SIZE_WF_V00)
#endif
/* *** */
#endif  /*END VARIANT_00*/
/*
                                               Restrictions
*/
#if UPDATE_NEURONS_PS_ORDER_LIMIT < 2
  #error (UPDATE_NEURONS_PS_ORDER_LIMIT must be >= 2)
#endif
/* *** */
#if UPDATE_NEURONS_ELEMENTS_PER_WI != 1
  #error (UPDATE_NEURONS_ELEMENTS_PER_WI has to be 1)
#endif
/* *** */
#if	((UPDATE_NEURONS_TOLERANCE_MODE == 2) && \
    (UPDATE_NEURONS_TOTAL_NEURONS < UPDATE_NEURONS_TOLERANCE_CHUNKS))
  #error (UPDATE_NEURONS_TOTAL_NEURONS has to be no less than UPDATE_NEURONS_TOLERANCE_CHUNKS)
#endif
/* *** */
#endif /*UPDATE_NEURONS_ENABLE*/
/*
                                               Default Definitions
*/
#if !defined (UPDATE_NEURONS_DEBUG_ENABLE)
  #define UPDATE_NEURONS_DEBUG_ENABLE                                   0
#endif
/* *** */
#if !defined (UPDATE_NEURONS_ERROR_TRACK_ENABLE)
  #define UPDATE_NEURONS_ERROR_TRACK_ENABLE                             0
#endif
/* *** */
#if !defined (UPDATE_NEURONS_VERIFY_ENABLE)
  #define UPDATE_NEURONS_VERIFY_ENABLE                                  0
#endif
/**************************************************************************************************/



/** ############################################################################################# **

  VII. Integration Restrictions
  
  (Note: This section depends on definitions in section VI)
  
  The purpose of Integration Restrictions is to verify consistency between related definitions
  between kernels when these kernels are integrated in a single program
  
** ############################################################################################# **/



/*If both EXPAND_EVENTS_ENABLE and UPDATE_NEURONS_ENABLE_V00 are enabled*/
#if EXPAND_EVENTS_ENABLE && (UPDATE_NEURONS_ENABLE_V00)
/* *** */
#if EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS != UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS
  #error (EXPAND_EVENTS_SPIKE_DATA_UNIT_SIZE_WORDS != UPDATE_NEURONS_SPIKE_DATA_UNIT_SIZE_WORDS)
#endif
/* *** */
#if EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE != UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE
  #error (EXPAND_EVENTS_SPIKE_DATA_BUFFER_SIZE != UPDATE_NEURONS_SPIKE_DATA_BUFFER_SIZE)
#endif
/* *** */
#if EXPAND_EVENTS_TOTAL_NEURONS != UPDATE_NEURONS_TOTAL_NEURONS
  #error (EXPAND_EVENTS_TOTAL_NEURONS != UPDATE_NEURONS_TOTAL_NEURONS)
#endif
/* *** */
#if EXPAND_EVENTS_SPIKE_PACKETS != UPDATE_NEURONS_SPIKE_PACKETS_V00
  #error (EXPAND_EVENTS_SPIKE_PACKETS != UPDATE_NEURONS_SPIKE_PACKETS_V00)
#endif
/* *** */
#if EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS != UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS
  #error (EXPAND_EVENTS_SPIKE_PACKET_SIZE_WORDS != UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS)
#endif
/* *** */
#endif

/*If both EXPAND_EVENTS_ENABLE and GROUP_EVENTS_ENABLE_VXX are enabled*/
#if (EXPAND_EVENTS_ENABLE && (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 ||\
    GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03))
/* *** */
#if (GROUP_EVENTS_TIME_SLOTS != EXPAND_EVENTS_TIME_SLOTS)
  #error (GROUP_EVENTS_TIME_SLOTS != EXPAND_EVENTS_TIME_SLOTS)
#endif
/* *** */
#if (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
  #error (GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
#endif
/* *** */
#if (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE != \
     EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
  #error (GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE != \
          EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
#endif
/* *** */
#if (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS != EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)
  #error (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS != EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS)
#endif
/* *** */
#if (GROUP_EVENTS_HISTOGRAM_BIN_SIZE != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
  #error (GROUP_EVENTS_HISTOGRAM_BIN_SIZE != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
#endif
/* *** */
#if (GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS != EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
  #error (GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS != EXPAND_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
/* *** */
#if (GROUP_EVENTS_TOTAL_NEURON_BITS && EXPAND_EVENTS_TOTAL_NEURON_BITS) &&\
    (GROUP_EVENTS_TOTAL_NEURON_BITS != EXPAND_EVENTS_TOTAL_NEURON_BITS)
  #error (GROUP_EVENTS_TOTAL_NEURON_BITS != EXPAND_EVENTS_TOTAL_NEURON_BITS)
#endif
/* *** */
#endif

/*If both SCAN_ENABLE_V00 and GROUP_EVENTS_ENABLE_VXX are enabled*/
#if (SCAN_ENABLE_V00 && (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 ||\
    GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03))
/* *** */
#if (GROUP_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
  #error (GROUP_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
#endif
/* *** */
#if (SCAN_SYNAPTIC_EVENT_BUFFERS != GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS)
  #error (SCAN_SYNAPTIC_EVENT_BUFFERS != GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS)
#endif
/* *** */
#if (SCAN_EVENT_DATA_MAX_SRC_BUFFER_SIZE != GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE)
  #error (SCAN_EVENT_DATA_MAX_SRC_BUFFER_SIZE != GROUP_EVENTS_EVENT_DATA_MAX_SRC_BUFFER_SIZE)
#endif
/* *** */
#if (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
  #error (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
#endif
/* *** */
#if (GROUP_EVENTS_HISTOGRAM_BIN_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V00)
  #error (GROUP_EVENTS_HISTOGRAM_BIN_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V00)
#endif
/* *** */
#endif

/*If both SCAN_ENABLE_V00 and EXPAND_EVENTS_ENABLE are enabled*/
#if (SCAN_ENABLE_V00 && EXPAND_EVENTS_ENABLE)
/* *** */
#if (EXPAND_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
  #error (EXPAND_EVENTS_TIME_SLOTS != SCAN_HISTOGRAM_TIME_SLOTS)
#endif
/* *** */
#if (SCAN_SYNAPTIC_EVENT_BUFFERS != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
  #error (SCAN_SYNAPTIC_EVENT_BUFFERS != EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS)
#endif
/* *** */
#if (SCAN_EVENT_DATA_MAX_SRC_BUFFER_SIZE != EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
  #error (SCAN_EVENT_DATA_MAX_SRC_BUFFER_SIZE != EXPAND_EVENTS_SYNAPTIC_EVENT_DATA_MAX_BUFFER_SIZE)
#endif
/* *** */
#if (EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
  #error (EXPAND_EVENTS_HISTOGRAM_TOTAL_BINS != SCAN_HISTOGRAM_TOTAL_BINS_V00)
#endif
/* *** */
#if (EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS != SCAN_HISTOGRAM_BIN_SIZE_V00)
  #error (EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS != SCAN_HISTOGRAM_BIN_SIZE_V00)
#endif
/* *** */
#endif

/*If both MAKE_EVENT_PTRS_ENABLE and GROUP_EVENTS_ENABLE_VXX are enabled*/
#if (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03) && MAKE_EVENT_PTRS_ENABLE
/* *** */
#if MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS
  #error (MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
/* *** */
#if GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE
  #error (GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE)
#endif
/* *** */
#if (MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET != \
    (GROUP_EVENTS_HISTOGRAM_BIN_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS))
  #error (MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET != \
         (GROUP_EVENTS_HISTOGRAM_BIN_SIZE*GROUP_EVENTS_HISTOGRAM_TOTAL_BINS))
#endif
/* *** */
#if (MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
  #error (MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
/* *** */
#endif

/*If both UPDATE_NEURONS_ENABLE_V00 and GROUP_EVENTS_ENABLE_VXX are enabled*/
#if (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03) && UPDATE_NEURONS_ENABLE_V00
/* *** */
#if UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS
  #error (UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS)
#endif
/* *** */
#if UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE
  #error (UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE)
#endif
/* *** */
#endif

/*If both UPDATE_NEURONS_ENABLE_V00 and MAKE_EVENT_PTRS_ENABLE are enabled*/
#if MAKE_EVENT_PTRS_ENABLE && UPDATE_NEURONS_ENABLE_V00
/* *** */
#if UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS
  #error (UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS != MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS)
#endif
/* *** */
#if UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE
  #error (UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE != MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE)
#endif
/* *** */
#endif

/*If both SCAN_ENABLE_V01 and GROUP_EVENTS_ENABLE_VXX are enabled*/
#if (((GROUP_EVENTS_ENABLE_V00) && SCAN_ENABLE_V01) ||\
     ((GROUP_EVENTS_ENABLE_V01) && SCAN_ENABLE_V01) ||\
     ((GROUP_EVENTS_ENABLE_V02) && SCAN_ENABLE_V01) ||\
     ((GROUP_EVENTS_ENABLE_V03) && SCAN_ENABLE_V01))
/* *** */
#if (GROUP_EVENTS_GRID_SIZE_WG != SCAN_HISTOGRAM_BIN_BACKETS)
  #error (GROUP_EVENTS_GRID_SIZE_WG != SCAN_HISTOGRAM_BIN_BACKETS)
#endif
/* *** */
#if (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V01)
  #error (GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE != SCAN_HISTOGRAM_BIN_SIZE_V01)
#endif
/* *** */
#if (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT != SCAN_HISTOGRAM_TOTAL_BINS_V01)
  #error (GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT != SCAN_HISTOGRAM_TOTAL_BINS_V01)
#endif
/* *** */
#endif

/*Verification of sort*/
#if SORT_VERIFY_ENABLE && \
    !(GROUP_EVENTS_ENABLE_V03 && GROUP_EVENTS_ENABLE_V02 && GROUP_EVENTS_ENABLE_V01 && \
     GROUP_EVENTS_ENABLE_V00 && SCAN_ENABLE_V01 && MAKE_EVENT_PTRS_ENABLE)
  #undef SORT_VERIFY_ENABLE
  #define SORT_VERIFY_ENABLE  0
#endif

/*Profiling mode*/
#if (PROFILING_MODE == 1) && KERNEL_VERIFY_ENABLE
  #error (cannot use (PROFILING_MODE = 1) if KERNEL_VERIFY_ENABLE is true)
#endif



/** ############################################################################################# **

  VIII. Enabling/Disabling Operators
  
  (Note: This section depends on definitions in section VI)
  
** ############################################################################################# **/



/*Enable Operator_Expand*/
#define ENABLE_OPERATOR_EXPAND            (EXPAND_EVENTS_ENABLE)

/*Enable Operator_Scan*/
#define ENABLE_OPERATOR_SCAN              (SCAN_ENABLE_V00 || SCAN_ENABLE_V01)

/*Enable Operator_Scan*/
#define ENABLE_OPERATOR_GROUP             (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || \
                                           GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03)

/*Enable Operator_Sort*/
#define ENABLE_OPERATOR_SORT              (ENABLE_OPERATOR_SCAN || ENABLE_OPERATOR_GROUP)



/** ############################################################################################# **

  IX. Enabling/Disabling Unit Tests
  
  (Note: This section depends on definitions in section VI)
  
** ############################################################################################# **/



/*Validate UNIT_TEST_MODE*/
#if !((UNIT_TEST_MODE == 0) || (UNIT_TEST_MODE == 1) || (UNIT_TEST_MODE == 2))
  #error (UNIT_TEST_MODE: provided mode is not supported)
#endif

/*Enable unit test for Operator_Expand*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_EXPAND_EVENTS                    (EXPAND_EVENTS_ENABLE)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_EXPAND_EVENTS                    (EXPAND_EVENTS_ENABLE && \
                                                            !UPDATE_NEURONS_ENABLE_V00)
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_EXPAND_EVENTS                    (EXPAND_EVENTS_ENABLE && \
                                                            !UPDATE_NEURONS_ENABLE_V00)
#if ENABLE_UNIT_TEST_EXPAND_EVENTS
  #error (ENABLE_UNIT_TEST_EXPAND_EVENTS: cannot enable integration tests in this mode)
#endif
#endif

/*Enable unit test for Operator_Group V00*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V00                 (GROUP_EVENTS_ENABLE_V00)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V00                 (GROUP_EVENTS_ENABLE_V00 && \
                                                            (!EXPAND_EVENTS_ENABLE || \
                                                            !SCAN_ENABLE_V00))
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V00                 (GROUP_EVENTS_ENABLE_V00 && \
                                                            (!EXPAND_EVENTS_ENABLE || \
                                                            !SCAN_ENABLE_V00))
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V00
  #error (ENABLE_UNIT_TEST_GROUP_EVENTS_V00: cannot enable integration tests in this mode)
#endif
#endif

/*Enable unit test for Operator_Group V01*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V01                 (GROUP_EVENTS_ENABLE_V01)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V01                 (GROUP_EVENTS_ENABLE_V01 && \
                                                            (!GROUP_EVENTS_ENABLE_V00 || \
                                                            !SCAN_ENABLE_V01 || \
                                                            !GROUP_EVENTS_ENABLE_V02))
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V01                 (GROUP_EVENTS_ENABLE_V01 && \
                                                            (!GROUP_EVENTS_ENABLE_V00 || \
                                                            !SCAN_ENABLE_V01 || \
                                                            !GROUP_EVENTS_ENABLE_V02))
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V01
  #error (ENABLE_UNIT_TEST_GROUP_EVENTS_V01: cannot enable integration tests in this mode)
#endif
#endif

/*Enable unit test for Operator_Group V02*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V02                 (GROUP_EVENTS_ENABLE_V02)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V02                 (GROUP_EVENTS_ENABLE_V02 && \
                                                            (!GROUP_EVENTS_ENABLE_V00 || \
                                                            !SCAN_ENABLE_V01 || \
                                                            !GROUP_EVENTS_ENABLE_V01))
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V02                 (GROUP_EVENTS_ENABLE_V02 && \
                                                            (!GROUP_EVENTS_ENABLE_V00 || \
                                                            !SCAN_ENABLE_V01 || \
                                                            !GROUP_EVENTS_ENABLE_V01))
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V02
  #error (ENABLE_UNIT_TEST_GROUP_EVENTS_V02: cannot enable integration tests in this mode)
#endif
#endif

/*Enable unit test for Operator_Group V03*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V03                 (GROUP_EVENTS_ENABLE_V03)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V03                 (GROUP_EVENTS_ENABLE_V03 && \
                                                            (!GROUP_EVENTS_ENABLE_V00 || \
                                                            !GROUP_EVENTS_ENABLE_V01 || \
                                                            !GROUP_EVENTS_ENABLE_V02 || \
                                                            !SCAN_ENABLE_V01))
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_GROUP_EVENTS_V03                 (GROUP_EVENTS_ENABLE_V03 && \
                                                            (!GROUP_EVENTS_ENABLE_V00 || \
                                                            !GROUP_EVENTS_ENABLE_V01 || \
                                                            !GROUP_EVENTS_ENABLE_V02 || \
                                                            !SCAN_ENABLE_V01))
#if ENABLE_UNIT_TEST_GROUP_EVENTS_V03
  #error (ENABLE_UNIT_TEST_GROUP_EVENTS_V03: cannot enable integration tests in this mode)
#endif
#endif

/*Enable unit test for Operator_Group*/
#define ENABLE_UNIT_TEST_GROUP_EVENTS                       (ENABLE_UNIT_TEST_GROUP_EVENTS_V00 || \
                                                            ENABLE_UNIT_TEST_GROUP_EVENTS_V01 || \
                                                            ENABLE_UNIT_TEST_GROUP_EVENTS_V02 || \
                                                            ENABLE_UNIT_TEST_GROUP_EVENTS_V03)

/*Enable unit test for Operator_Scan V00*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_SCAN_V00                         (SCAN_ENABLE_V00)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_SCAN_V00                         (SCAN_ENABLE_V00 && \
                                                            !EXPAND_EVENTS_ENABLE)
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_SCAN_V00                         (SCAN_ENABLE_V00 && \
                                                            !EXPAND_EVENTS_ENABLE)
#if ENABLE_UNIT_TEST_SCAN_V00
  #error (ENABLE_UNIT_TEST_SCAN_V00: cannot enable integration tests in this mode)
#endif
#endif

/*Enable unit test for Operator_Scan V01*/
#if (UNIT_TEST_MODE == 0)
  #define ENABLE_UNIT_TEST_SCAN_V01                         (SCAN_ENABLE_V01)
#elif (UNIT_TEST_MODE == 1)
  #define ENABLE_UNIT_TEST_SCAN_V01                         (SCAN_ENABLE_V01 && (\
                                                            !(GROUP_EVENTS_ENABLE_V00) || \
                                                            !(GROUP_EVENTS_ENABLE_V01 && \
                                                            GROUP_EVENTS_ENABLE_V02)))
#elif (UNIT_TEST_MODE == 2)
  #define ENABLE_UNIT_TEST_SCAN_V01                         (SCAN_ENABLE_V01 && (\
                                                            !(GROUP_EVENTS_ENABLE_V00) || \
                                                            !(GROUP_EVENTS_ENABLE_V01 && \
                                                            GROUP_EVENTS_ENABLE_V02)))
#if ENABLE_UNIT_TEST_SCAN_V01
  #error (ENABLE_UNIT_TEST_SCAN_V01: cannot enable integration tests in this mode)
#endif
#endif



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif /*INC_DEFINITIONS_H*/
