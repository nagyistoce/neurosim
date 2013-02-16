
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

#include "Definitions.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
  #include <windows.h>
  #include <winbase.h>
#elif (SYSTEM_OS == SYSTEM_OS_LINUX)
  #include <sys/time.h>
  #include <linux/limits.h>
  #include <unistd.h>
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
#include <CL/cl.hpp>
#include <CL/opencl.h>

#include "integ_util.h"
#include "iz_util.h"
#include "Data_Connectome.hpp"
#include "Data_SpikeEvents.hpp"
#include "Data_SynapticEvents.hpp"
#include "Operator_Expand.hpp"
#include "Operator_Group.hpp"
#include "Operator_Scan.hpp"
#include "Operator_Sort.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Namespaces
***************************************************************************************************/

using namespace std;

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



/***************************************************************************************************
  Structs
  
  TODO: may need to redefine these with constructor/destructor. 
  For example (from http://stackoverflow.com/questions/4203010/how-to-initialize-member-struct-
  in-initializer-list-of-c-class):
  
  struct Foo {

    Foo(int const a, std::initializer_list<char> const b, short* c)
      : x(a), y(c) {
      assert(b.size() >= 24, "err");
      std::copy(b.begin(), b.begin() + 24, array);
    }

    ~Foo() { delete y; }

    int x;
    char array[24];
    short* y;
  };

  // Class, which has Foo as its field:
  class Bar {

    Bar() : x(5), foo(5, {'a', 'b', ..., 'y', 'z'},
      new short(5)) { }

    private:

    int x;
    Foo foo;
  };
  
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
  
  /*Constructor*/
  SimException
  (
    std::string ss
  ) : 
    s(ss) 
  {};
  
  /*Destructor*/
  ~SimException()  throw ()
  {};
  
  const char* what() const throw()
  { 
    return s.c_str(); 
  }
};

/**************************************************************************************************/



/***************************************************************************************************
  Variables
***************************************************************************************************/



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
#if LOG_SIMULATION
    std::stringstream*,
#endif
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



/**************************************************************************************************/
  /**
    Find target device from the list of target devices.
  */
  bool 
  findTargetDevice
  (
    vector<cl::Device>                  platformDevices,
    const char                          *targetDevices,
    std::vector<cl::Device>::iterator   *d
  );
/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/



#endif // COMMON_HPP_
