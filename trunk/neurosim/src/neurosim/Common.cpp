
/* ===============================================================================================

  code snippets:
  
  for(cl_uint i = 0; i < EXPAND_EVENTS_GRID_SIZE_WG*EXPAND_EVENTS_WG_SIZE_WF; i++)
  {
    std::stringstream ss;
    for(cl_uint j = 0; j < EXPAND_EVENTS_WF_SIZE_WI; j++)
    {
      ss << dataExpandEventsDebugDevice[i*EXPAND_EVENTS_WF_SIZE_WI + j] << ",";
    }
    LOG_SIM(ss.str(), 0)
  }
  
  =============================================================================================== */
  

  
/***************************************************************************************************
  Includes
***************************************************************************************************/

#include "Common.hpp"

/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_START

/**************************************************************************************************/



std::string
getPath
()
/**************************************************************************************************/
{
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
  char buffer[MAX_PATH];
#ifdef UNICODE
  if(!GetModuleFileName(NULL, (LPWCH)buffer, sizeof(buffer)))
      throw std::string("GetModuleFileName() failed!");
#else
  if(!GetModuleFileName(NULL, buffer, sizeof(buffer)))
      throw std::string("GetModuleFileName() failed!");
#endif
  std::string str(buffer);
  /* '\' == 92 */
  size_t last = str.find_last_of((char)92);
#elif (SYSTEM_OS == SYSTEM_OS_LINUX)
  char buffer[PATH_MAX + 1];
  ssize_t len;
  if((len = readlink("/proc/self/exe",buffer, sizeof(buffer) - 1)) == -1)
      throw std::string("readlink() failed!");
  buffer[len] = '\0';
  std::string str(buffer);
  /* '/' == 47 */
  size_t last = str.find_last_of((char)47);
#endif
  return str.substr(0, last + 1);
}
/**************************************************************************************************/



bool
getSource
(
  const char  *fileName,
  std::string &source
)
/**************************************************************************************************/
{
  std::streampos size;
  char*         str;

  /* Open file stream*/
  std::fstream f(fileName, (std::fstream::in | std::fstream::binary));

  /* Check if we have opened file stream*/
  if (f.is_open()) 
  {
    std::streampos sizeFile;
    /* Find the stream size*/
    f.seekg(0, std::fstream::end);
    size = sizeFile = f.tellg();
    f.seekg(0, std::fstream::beg);

    if(size < 0) 
    {
      f.close();
      return  false;
    }
    
    str = new char[(size_t)size + 1];
    
    if(!str) 
    {
      f.close();
      return  false;
    }

    /* Read file*/
    f.read(str, sizeFile);
    f.close();
    str[size] = '\0';

    source = str;

    delete[] str;

    return true;
  }

  return false;
}
/**************************************************************************************************/



void 
replaceNewlineWithSpaces
(
  std::string &string
)
/**************************************************************************************************/
{
  size_t pos = string.find_first_of('\n', 0);
  
  while(pos != string::npos)
  {
      string.replace(pos, 1, " ");
      pos = string.find_first_of('\n', pos + 1);
  }
  
  pos = string.find_first_of('\r', 0);
  
  while(pos != string::npos)
  {
      string.replace(pos, 1, " ");
      pos = string.find_first_of('\r', pos + 1);
  }
}
/**************************************************************************************************/



void 
createKernel
(
#if LOG_SIMULATION
  std::stringstream   *dataToSimulationLogFile,
#endif
  cl::Context         &context,
  cl::Device          &device,
  cl::Kernel          &kernel,
  std::string         kernelFileName,
  const char          *kernelName,
  std::string         flagsStr,
  size_t              blockSizeX,
  size_t              blockSizeY
)
/**************************************************************************************************/
{
  cl_int err = CL_SUCCESS;
  cl::Program program;

  /* create a CL program using the kernel source */

  std::string kernelPath = getPath();
  std::string kernelSrc("");
  
  std::string defPath = getPath();
  std::string defSrc("");

  std::string flagsPath = OCL_COMPILER_OPTIONS_FILE_NAME;
  std::string flagsSrc("");

  kernelPath.append(kernelFileName);
  
  if(!getSource(kernelPath.c_str(), kernelSrc))
  {
    THROW_SIMEX("Neurosim::createKernel: Failed to load kernel file: " << kernelPath);
  }
  
  defPath.append("Definitions.h");
  
  if(!getSource(defPath.c_str(), defSrc))
  {
    THROW_SIMEX("Neurosim::createKernel: Failed to load kernel file: " << defPath);
  }

  /*To nearest 4B chunk*/
  size_t sourceCodeSize = ((defSrc.size() + kernelSrc.size())/sizeof(int)+1)*sizeof(int);
  /*size_t sourceCodeSize = defSrc.size() + kernelSrc.size();*/
  LOG_SIM("Source code size (" << kernelName << "): " << ((float)sourceCodeSize)/1024.0 << " KB");
    
  char *sourceCode = (char *) calloc(sourceCodeSize, sizeof(char));
  sourceCode[0] = '\0';
  strcat ( sourceCode, defSrc.data() );
  strcat ( sourceCode, kernelSrc.data() );
  /*std::cout << "Source code:\n" << sourceCode << std::endl;*/
  cl::Program::Sources programSource(1, std::make_pair(sourceCode, sourceCodeSize));

  program = cl::Program(context, programSource, &err);
  
  free(sourceCode);
  
  if(err != CL_SUCCESS)
  {
    THROW_SIMEX("Neurosim::createKernel: Failed cl::Program(Source) due to error code " << err);
  }

  if(!getSource(flagsPath.c_str(), flagsSrc))
  {
    THROW_SIMEX("Neurosim::createKernel: Failed to load flags file: " << flagsPath);
  }
  
  replaceNewlineWithSpaces(flagsSrc);
  const char *flags = flagsSrc.c_str();
  flagsStr.append(flags);
  if(flagsStr.size() != 0)
  {
    LOG_SIM("Build Options (" << kernelName << "): " << flagsStr.c_str());
  }
  
  vector<cl::Device> deviceVector;
  deviceVector.push_back(device);
  LOG_SIM("Building for device: " << (device).getInfo<CL_DEVICE_NAME>());
  err = program.build(deviceVector, flagsStr.c_str());
  
  if(err != CL_SUCCESS)
  {
    if(err == CL_BUILD_PROGRAM_FAILURE)
    {
      std::string str = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
      
      std::cout << " \n\tBUILD LOG for kernel " << kernelName << "\n";
      std::cout << "*******************************************************************\n";
      std::cout << str << std::endl;
      std::cout << "*******************************************************************\n";
    }
  }

  if(err != CL_SUCCESS)
  {
    THROW_SIMEX("Neurosim::createKernel: Failed cl::Program:build due to error code " << err);
  }
  
  /* Create kernel */
  kernel = cl::Kernel(program, kernelName, &err);

  if(err != CL_SUCCESS)
  {
    THROW_SIMEX("Neurosim::createKernel: Failed cl::Kernel() due to error code " << err);
  }

  /* Check group size against group size returned by kernel */
  size_t kernelWorkGroupSize = 
    kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device, &err);
    
  if(err != CL_SUCCESS)
  {
    THROW_SIMEX("Neurosim::createKernel: Failed Kernel::getWorkGroupInfo() due to error code " 
      << err);
  }

  if((blockSizeX * blockSizeY) > kernelWorkGroupSize)
  {
    THROW_SIMEX("Neurosim::createKernel: Out of Resources!" << std::endl
      << "Group Size specified : " << blockSizeX * blockSizeY << std::endl
      << "Max Group Size supported on the kernel : " << kernelWorkGroupSize);
  }
}
/**************************************************************************************************/



bool 
findTargetDevice
(
  vector<cl::Device>                  platformDevices,
  const char                          *targetDevices,
  std::vector<cl::Device>::iterator   *d
)
/**************************************************************************************************/
{
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
/**************************************************************************************************/



/***************************************************************************************************
  Warning control
***************************************************************************************************/

WARNING_CONTROL_END

/**************************************************************************************************/
