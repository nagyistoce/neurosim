/* ===============================================================================================

  code snippets:
  
  for(cl_uint i = 0; i < EXPAND_EVENTS_GRID_SIZE_WG*EXPAND_EVENTS_WG_SIZE_WF; i++)
  {
    std::stringstream ss;
    for(cl_uint j = 0; j < EXPAND_EVENTS_WF_SIZE_WI; j++)
    {
      ss << dataExpandEventsDebugDevice[i*EXPAND_EVENTS_WF_SIZE_WI + j] << ",";
    }
    LOG(ss.str(), 0)
  }
  
  =============================================================================================== */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <map>
#include <set>
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <windows.h>
#include <winbase.h>
#include <exception>

#include "Definitions.h"
#include "integ_util.h"
#include "iz_util.h"
#include <CL/cl.hpp>
#include <SDKCommon.hpp>
#include <SDKApplication.hpp>
#include <SDKFile.hpp>
#include <SDKBitMap.hpp>
#include "SpikeEvents.hpp"
#include "Connectome.hpp"



using namespace std;



/** ############################################################################################# **

  Structs
  
** ############################################################################################# **/

struct kernelStatistics 
{
  /*map{kernel_name -> map{alloc_name -> alloc_size}}*/
  map<std::string, map<std::string, cl_uint> > gmSizes;
  map<std::string, map<std::string, cl_uint> > cmSizes;
  map<std::string, map<std::string, cl_uint> > lmSizes;
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



/** ############################################################################################# **

  Templates
  
** ############################################################################################# **/

template<typename T>
void swap2(T& a, T& b)
{
  T tmp = a;
  a = b;
  b = tmp;
}

#endif // DEFINITIONS_H_
