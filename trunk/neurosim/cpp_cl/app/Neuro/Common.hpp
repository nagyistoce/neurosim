
#ifndef COMMON_HPP_
#define COMMON_HPP_



/**
  @file Common.hpp

  Defines methods used in any class

  @author Dmitri Yudanov, dxy7370@gmail.com

  @date 2012/05/17
*/



/***************************************************************************************************
  Includes
***************************************************************************************************/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>

#ifdef _WIN32
  #include <windows.h>
  #include <winbase.h>
#else
  #include <sys/time.h>
  #include <linux/limits.h>
#endif

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <ctime>
#include <exception>

#include "Definitions.h"
#include "integ_util.h"
#include "iz_util.h"
#include <CL/cl.hpp>
#include <CL/opencl.h>
#include <SDKCommon.hpp>
#include <SDKApplication.hpp>
#include <SDKFile.hpp>
#include <SDKBitMap.hpp>
#include "OperatorScan.hpp"
#include "SpikeEvents.hpp"
#include "Connectome.hpp"
#include "SynapticEvents.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Namespaces
***************************************************************************************************/

using namespace std;

/**************************************************************************************************/



/***************************************************************************************************
  Structs
***************************************************************************************************/

struct kernelStatistics 
{
  /*map{kernel_name -> map{alloc_name -> alloc_size}}*/
  map<std::string, map<std::string, size_t> > gmSizes;
  map<std::string, map<std::string, size_t> > cmSizes;
  map<std::string, map<std::string, size_t> > lmSizes;
  set<std::string> kernelNames;
  map<std::string, map<std::string, double> > execTime;
  set<std::string> kernelNamesExecTime;
};



struct SimException : public std::exception
{
  std::string s;
  SimException(std::string ss) : s(ss) {}
  const char* what() const throw() { return s.c_str(); }
};

/**************************************************************************************************/



/***************************************************************************************************
  Templates
***************************************************************************************************/

template<typename T>
void swap2(T& a, T& b)
{
  T tmp = a;
  a = b;
  b = tmp;
}

/**************************************************************************************************/



/***************************************************************************************************
  Methods
***************************************************************************************************/

/**************************************************************************************************/
  /**
    Returns the path of executable that calls this method.
  */
  std::string
  getPath
  ();
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Returns the contents of source code file by reference.
  */
  bool
  getSource
  (
    const char*,
    std::string&
  );
/**************************************************************************************************/



/**************************************************************************************************/
  /**
    Replaces new line with spaces.
  */
  void 
  replaceNewlineWithSpaces
  (
    std::string &string
  );
/**************************************************************************************************/

  

/**************************************************************************************************/
  /**
    Creates an OpenCL kernel and returns it by reference
  */
  void
  createKernel
  (
    cl::Context&,
    cl::Device&,
    cl::Kernel&,
    std::string,
    const char*,
    std::string,
    size_t,
    size_t
  );
/**************************************************************************************************/

#endif // COMMON_HPP_
