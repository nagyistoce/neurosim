
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



#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
void
Neurosim::allocateHostData
(
  cl::Device  &device
)
/**************************************************************************************************/
{
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
  cl_uint size = 0;
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
  memset((*synapticEvents).dataSortedEventsHistogram, 0, (*synapticEvents).dataSortedEventsHistogramSizeBytes);
  memset((*synapticEvents).dataSortedEvents, 0, (*synapticEvents).dataSortedEventsSizeBytes);
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

  (*synapticEvents).dataSortedEventsHistogram[MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET] = totalEvents;
  
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
    (cl_uint)((*synapticEvents).dataSortedEventsSize)/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS,
    (*synapticEvents).dataSortedEvents
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
    (*synapticEvents).dataSortedEventsSize/MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS,
    (*synapticEvents).dataSortedEvents
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
    memset((*synapticEvents).dataSortedEvents, 0, (*synapticEvents).dataSortedEventsSizeBytes);
    result = initializeSortedEvents
    (
      2,
      (UPDATE_NEURONS_TOTAL_NEURONS-1),
      30.0,
      5.0,
      gabaRatio,
      (cl_uint)((*synapticEvents).dataSortedEventsSize)/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
      0,
      0,
      0,
      (cl_uint)((*synapticEvents).dataSortedEventsSize)/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
      (*synapticEvents).dataSortedEvents
    );
    if(result != 1){return result;}
    
    memset(dataMakeEventPtrsStruct, 0, dataMakeEventPtrsStructSizeBytes);
    
    result = initializeEventPointers
    (
      true,
      (cl_uint)((*synapticEvents).dataSortedEventsSize)/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
      UPDATE_NEURONS_TOTAL_NEURONS,
      UPDATE_NEURONS_STRUCT_ELEMENT_SIZE,
      (*synapticEvents).dataSortedEvents,
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



#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
void
Neurosim::registerLocalMemory
(
  cl::Device  &device
)
/**************************************************************************************************/
{
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
  {
  size_t lmGenericMakeEventPtrs = sizeof(cl_uint)*MAKE_EVENT_PTRS_WF_LM_SHARE_SIZE; 
  size_t lmGenericMakeEventPtrsSizeBytes = lmGenericMakeEventPtrs;
  REGISTER_MEMORY(MAKE_EVENT_PTRS_KERNEL_NAME, MEM_LOCAL, lmGenericMakeEventPtrs);
  
#if MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE == 0
  size_t lmLastElementMakeEventPtrs = sizeof(cl_uint)*
    (MAKE_EVENT_PTRS_WG_SIZE_WF*MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE); 
  size_t lmLastElementMakeEventPtrsSizeBytes = lmLastElementMakeEventPtrs;
  REGISTER_MEMORY(MAKE_EVENT_PTRS_KERNEL_NAME, MEM_LOCAL, lmLastElementMakeEventPtrs);
  
  size_t lmGenericGlueEventPtrs = sizeof(cl_uint)*GLUE_EVENT_WF_LM_SHARE_SIZE; 
  size_t lmGenericGlueEventPtrsSizeBytes = lmGenericGlueEventPtrs;
  REGISTER_MEMORY(GLUE_EVENT_PTRS_KERNEL_NAME, MEM_LOCAL, lmGenericGlueEventPtrs);
#endif
  }
#endif
/**************************************************************************************************/
#if UPDATE_NEURONS_ENABLE_V00
  {
  std::string kernelTag = std::string(UPDATE_NEURONS_KERNEL_NAME) + std::string("_V00");

  size_t lmSpikePackets = sizeof(cl_uint)*UPDATE_NEURONS_WG_SIZE_WF_V00*
    UPDATE_NEURONS_SPIKE_PACKET_SIZE_WORDS; 
  size_t lmSpikePacketsSizeBytes = lmSpikePackets;
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
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
  this->allocateHostData(*d);
#endif

    /*Initialize data objects*/

#if EXPAND_EVENTS_ENABLE || GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || \
  GROUP_EVENTS_ENABLE_V02 || GROUP_EVENTS_ENABLE_V03
  
  cl_uint sortedEventBufferCount = 0;
  cl_uint sortedEventBufferSize = 0;
  cl_uint sortedEventsHistogramBacketCount = 0;
  cl_uint sortedEventsHistogramBinSize = 0;
  cl_uint sortedEventsHistogramBinCount = 0;

#if UPDATE_NEURONS_ENABLE_V00
  sortedEventBufferSize = UPDATE_NEURONS_TEST_MAX_SRC_BUFFER_SIZE;
  sortedEventBufferCount = UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS;
#elif MAKE_EVENT_PTRS_ENABLE
  sortedEventBufferSize = MAKE_EVENT_PTRS_TEST_MAX_SRC_BUFFER_SIZE;
  sortedEventBufferCount = MAKE_EVENT_PTRS_EVENT_DATA_PITCH_WORDS;
#elif (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03)
  sortedEventBufferSize = GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE;
  sortedEventBufferCount = GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS;
#endif

#if MAKE_EVENT_PTRS_ENABLE
  sortedEventsHistogramBacketCount = MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET+1;
  sortedEventsHistogramBinSize = 1;
  sortedEventsHistogramBinCount = 1;
#endif

#if (GROUP_EVENTS_ENABLE_V00 || GROUP_EVENTS_ENABLE_V01 || GROUP_EVENTS_ENABLE_V02 ||\
    GROUP_EVENTS_ENABLE_V03)
  sortedEventsHistogramBacketCount = GROUP_EVENTS_GRID_SIZE_WG;
  sortedEventsHistogramBinSize = GROUP_EVENTS_HISTOGRAM_OUT_GRID_SIZE;
  sortedEventsHistogramBinCount = GROUP_EVENTS_HISTOGRAM_TOTAL_BINS_OUT;
#endif

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
    EXPAND_EVENTS_SYNAPTIC_EVENT_BUFFERS,
    sortedEventsHistogramBacketCount,
    sortedEventsHistogramBinSize,
    sortedEventsHistogramBinCount,
    sortedEventBufferCount,
    sortedEventBufferSize
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
    GROUP_EVENTS_HISTOGRAM_BIN_SIZE,
    sortedEventsHistogramBacketCount,
    sortedEventsHistogramBinSize,
    sortedEventsHistogramBinCount,
    sortedEventBufferCount,
    sortedEventBufferSize
  );
#endif
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
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF),
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
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF),
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
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -cl-unsafe-math-optimizations "
      "-cl-finite-math-only -cl-fast-relaxed-math",
#else
      /*"-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt -cl-opt-disable",*/
      /*"-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1",*/
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1",
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
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -cl-unsafe-math-optimizations "
      "-cl-finite-math-only -cl-fast-relaxed-math",
#else
      /*"-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt -cl-opt-disable",*/
      /*"-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt",*/
      "-D SYSTEM=" TOSTRING(SYSTEM_CONTROL_OFF) " -D UPDATE_NEURONS_DEVICE_V00=1 "
      "-cl-fp32-correctly-rounded-divide-sqrt",
#endif
      blockSizeX_kernelUpdateNeuronsV00,
      blockSizeY_kernelUpdateNeuronsV00
    );
  }
#endif

    /*Initialize operator objects*/
    
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

#if ENABLE_OPERATOR_SORT
  {
    cl_uint sortByTimeIterationCount = 0;
    cl_uint sortByNeuronIterationCount = 0;
    
#if ENABLE_OPERATOR_GROUP
    sortByTimeIterationCount = (32/GROUP_EVENTS_HISTOGRAM_BIN_BITS);
    sortByNeuronIterationCount = GROUP_EVENTS_NEURON_ID_SORT_ITERATIONS;
#endif
    
    this->operatorSort = new Operator_Sort
    (
#if ENABLE_UNIT_TEST_GROUP_EVENTS
      GROUP_EVENTS_TEST_MODE,
#endif
#if (SCAN_DEBUG_ENABLE || SCAN_ERROR_TRACK_ENABLE || ENABLE_UNIT_TEST_SCAN_V00 || \
  ENABLE_UNIT_TEST_SCAN_V01) || (GROUP_EVENTS_DEBUG_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE)
      CL_FALSE,
      commandQueue,
#endif
      sortByTimeIterationCount,
      sortByNeuronIterationCount,
      context,
      *d,
      kernelStats,
      synapticEvents,
      dataToSimulationLogFile,
      dataToReportLogFile
    );
  }
#endif

  /* Register local memory for stats */
#if MAKE_EVENT_PTRS_ENABLE || UPDATE_NEURONS_ENABLE_V00
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
      -1.0
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
  ENQUEUE_WRITE_BUFFER(CL_FALSE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, 
    (*synapticEvents).dataSortedEventsSizeBytes, (*synapticEvents).dataSortedEvents);
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
#if ENABLE_OPERATOR_SORT
    (*operatorSort).sortEvents
    (
#if KERNEL_LEVEL_PROFILING
      kernelStats,
#endif
#if GROUP_EVENTS_VERIFY_ENABLE
      currentTimeSlot,
#endif
#if KERNEL_LEVEL_PROFILING || SCAN_ERROR_TRACK_ENABLE || GROUP_EVENTS_ERROR_TRACK_ENABLE || \
(ENABLE_OPERATOR_SCAN && SCAN_ENABLE_V00) || (ENABLE_OPERATOR_GROUP && GROUP_EVENTS_ENABLE_V00)
      currentTimeStep,
#endif
      commandQueue,
      ndrEvt
    );
#endif
/**************************************************************************************************/    
#if SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataUnsortedEventWeightsBuffer, 
        (*synapticEvents).dataUnsortedEventWeightsSizeBytes, 
        (*synapticEvents).dataUnsortedEventWeights);
    }
#endif
/**************************************************************************************************/
#if MAKE_EVENT_PTRS_ENABLE
    {
/*Unit test initialization*/
#if !(GROUP_EVENTS_ENABLE_V00 && GROUP_EVENTS_ENABLE_V01 && GROUP_EVENTS_ENABLE_V02 &&\
      GROUP_EVENTS_ENABLE_V03 && SCAN_ENABLE_V01)
#error (MAKE_EVENT_PTRS_ENABLE unit test is temporary unavailable because *synapticEvents cannot be initialized)
    /*Use milder mode if UPDATE_NEURONS_ENABLE_V00*/
    int mode = (currentTimeStep == 0) ? 0 : (UPDATE_NEURONS_ENABLE_V00 ? 1 : 2);
    
    if(initializeDataForKernelMakeEventPtrs(mode, currentTimeStep-1) != 1)
    {
      throw SimException("Neurosim::run: Failed initializeDataForKernelMakeEventPtrs");
    }
    
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsHistogramWriteBuffer, 
      (*synapticEvents).dataSortedEventsHistogramSizeBytes, (*synapticEvents).dataSortedEventsHistogram);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, 
      (*synapticEvents).dataSortedEventsSizeBytes, (*synapticEvents).dataSortedEvents);
/*End unit test initialization*/

#elif MAKE_EVENT_PTRS_VERIFY_ENABLE
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsHistogramReadBuffer, 
      (*synapticEvents).dataSortedEventsHistogramSizeBytes, (*synapticEvents).dataSortedEventsHistogram);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, 
      (*synapticEvents).dataSortedEventsSizeBytes, (*synapticEvents).dataSortedEvents);
#endif

#if (MAKE_EVENT_PTRS_DEBUG_ENABLE)
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsDebugHostBuffer, 
      dataMakeEventPtrsDebugHostSizeBytes, dataMakeEventPtrsDebugHost);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsDebugDeviceBuffer, 
      dataMakeEventPtrsDebugDeviceSizeBytes, dataMakeEventPtrsDebugDevice);
#endif

    SET_KERNEL_ARG(kernelMakeEventPtrs, (*synapticEvents).dataSortedEventsHistogramWriteBuffer, argNumMakeEventPtrs);
    SET_KERNEL_ARG(kernelMakeEventPtrs, (*synapticEvents).dataSortedEventsReadBuffer, argNumMakeEventPtrs+1);
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
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, 
      (*synapticEvents).dataSortedEventsSizeBytes, (*synapticEvents).dataSortedEvents);
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, 
      dataMakeEventPtrsStructSizeBytes, dataMakeEventPtrsStruct);

    (*operatorSort).captureUnsortedEvents
    (
      currentTimeSlot,
      GROUP_EVENTS_SYNAPTIC_EVENT_BUFFERS,
      GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
      GROUP_EVENTS_MAX_DELAY,
      GROUP_EVENTS_MIN_DELAY,
      commandQueue,
      (*synapticEvents)
    );

    (*operatorSort).verifySortedEvents
    (
      currentTimeStep,
      2,
      GROUP_EVENTS_EVENT_DATA_UNIT_SIZE_WORDS,
      GROUP_EVENTS_EVENT_DATA_MAX_DST_BUFFER_SIZE,
      MAKE_EVENT_PTRS_STRUCTS,
      MAKE_EVENT_PTRS_STRUCT_SIZE,
      MAKE_EVENT_PTRS_STRUCT_ELEMENT_SIZE,
      MAKE_EVENT_PTRS_EVENT_DELIVERY_MODE,
      GROUP_EVENTS_MAX_DELAY,
      GROUP_EVENTS_MIN_DELAY,
      (*synapticEvents).dataSortedEvents,
      dataMakeEventPtrsStruct
    );
#elif SIMULATION_SNAPSHOT
    if(takeSimSnapshot)
    {
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, 
        (*synapticEvents).dataSortedEventsSizeBytes, (*synapticEvents).dataSortedEvents);
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
#error (UPDATE_NEURONS_ENABLE_V00 unit test is temporary unavailable because *synapticEvents cannot be initialized)
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
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
      (*synapticEvents).dataSortedEvents);
    ENQUEUE_WRITE_BUFFER(CL_TRUE, commandQueue, dataMakeEventPtrsStructBuffer, dataMakeEventPtrsStructSizeBytes, 
      dataMakeEventPtrsStruct);
# else
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
      (*synapticEvents).dataSortedEvents);
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
    SET_KERNEL_ARG(kernelUpdateNeuronsV00, (*synapticEvents).dataSortedEventsReadBuffer, argNumUpdateNeuronsV00);
    SET_KERNEL_ARG(kernelUpdateNeuronsV00, currentTimeStep, argNumUpdateNeuronsV00+1);

    ENQUEUE_KERNEL_V0(currentTimeStep, kernelStats, commandQueue, kernelUpdateNeuronsV00, 
      globalThreadsUpdateNeuronsV00, localThreadsUpdateNeuronsV00, NULL, ndrEvt, 
      "kernelUpdateNeuronsV00");
      
#if DEVICE_HOST_DATA_COHERENCE
    (*spikeEvents).invalidateEvents();
#endif
/*END enqueue update kernel*/

/*Set loop-dependent arguments for update kernel for spiked neurons and equeue it*/
    SET_KERNEL_ARG(kernelUpdateSpikedNeuronsV00, (*synapticEvents).dataSortedEventsReadBuffer, 
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
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
        (*synapticEvents).dataSortedEvents);
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
          (*synapticEvents).dataSortedEventsSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          (*synapticEvents).dataSortedEvents,
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
          (cl_uint)((*synapticEvents).dataSortedEventsSize)/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          (*synapticEvents).dataSortedEvents,
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
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
        (*synapticEvents).dataSortedEvents);

      if
      (
        verifyKernelUpdateNeurons
        (
          !COMPILER_FLAGS_OPTIMIZE_ENABLE,
          true,
          true,
          currentTimeStep,
          (*synapticEvents).dataSortedEventsSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          (*synapticEvents).dataSortedEvents,
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
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
        (*synapticEvents).dataSortedEvents);

      if
      (
        verifyKernelUpdateNeurons
        (
          !COMPILER_FLAGS_OPTIMIZE_ENABLE,
          false,
          true,
          currentTimeStep,
          (*synapticEvents).dataSortedEventsSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          (*synapticEvents).dataSortedEvents,
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
      ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
        (*synapticEvents).dataSortedEvents);
        
      if
      (
        verifyKernelUpdateNeurons
        (
          !COMPILER_FLAGS_OPTIMIZE_ENABLE,
          false,
          true,
          currentTimeStep,
          (*synapticEvents).dataSortedEventsSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
          dataMakeEventPtrsStruct,
          (*synapticEvents).dataSortedEvents,
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
        (*synapticEvents).dataSortedEventsSize/UPDATE_NEURONS_EVENT_DATA_PITCH_WORDS,
        dataMakeEventPtrsStruct,
        (*synapticEvents).dataSortedEvents,
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
    ENQUEUE_READ_BUFFER(CL_TRUE, commandQueue, (*synapticEvents).dataSortedEventsReadBuffer, (*synapticEvents).dataSortedEventsSizeBytes, 
      (*synapticEvents).dataSortedEvents);
        
    takeSimulationSnapshot
    (
      currentTimeStep,
      1000,
      dataMakeEventPtrsStruct,
      (*synapticEvents).dataSortedEvents,
      modelVariables,
      (*spikeEvents),
      (*synapticEvents),
      commandQueue
    );
  }
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



#if MAKE_EVENT_PTRS_ENABLE
int 
Neurosim::verifyKernelMakeEventPtrs()
/**************************************************************************************************/
{
  int result = 1;
  
  /*Load total event for the test*/
  cl_uint totalEvents = (*synapticEvents).dataSortedEventsHistogram[MAKE_EVENT_PTRS_TOTAL_EVENTS_OFFSET];

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
      cl_uint previousElement = (*synapticEvents).dataSortedEvents[previousElementOffset];
      
      for(cl_uint j = 0; j < totalEvents; j++)
      {
        /*Detect boundary*/
        if((*synapticEvents).dataSortedEvents[j] != previousElement)
        {
          cl_uint offsetVerify = dataMakeEventPtrsStruct[(*synapticEvents).dataSortedEvents[j]*
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
              << (*synapticEvents).dataSortedEvents[j] << ": pointer (" << offsetVerify << " vs " << j 
              << "), count (" << countVerify << " vs " << j-previousElementOffset << ")"
              << std::endl;
            result = 0;
          }
          
          previousElement = (*synapticEvents).dataSortedEvents[j];
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
  cl_uint error = ((*synapticEvents).dataSortedEvents[j] != dataMakeEventPtrsStruct[1 + offsetPtr]) && count;
  error += (j != dataMakeEventPtrsStruct[1 + offsetPtr + 1]) && count;
  
  if(!error && (result != 0))
  {
    for(j = 1; j < totalEvents; j++)
    {
      /*Detect boundary*/
      if((*synapticEvents).dataSortedEvents[j] != (*synapticEvents).dataSortedEvents[j-1])
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
        error = ((*synapticEvents).dataSortedEvents[j] != 
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
      << ", global element ID " << j << ": key " << (*synapticEvents).dataSortedEvents[j] << " vs " 
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
#if UPDATE_NEURONS_ENABLE_V00
int 
Neurosim::injectSortedEvents
(
  bool          verify,
  bool          resetEventBuffer,
  cl_uint       totalNeurons,
  cl_uint       eventQueueSize,
  cl_uint       pointersPitch,
  size_t        sortedEventsSize,
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
#endif



/*
  Inject unsorted synaptic events into event queue
*/
#if UPDATE_NEURONS_ENABLE_V00
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
#endif



/*
  Propagation of spike events to synaptic events for PS method
*/
#if UPDATE_NEURONS_ENABLE_V00
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
#endif



/*
  Verification of synaptic event between host and device.
*/
#if UPDATE_NEURONS_ENABLE_V00
int 
Neurosim::verifyEvents
(
  bool          ignoreWarnings,
  bool          correctWeightPositionMismatch,
  unsigned int  totalNeurons,
  unsigned int  structElementSize,
  size_t        sortedEventsSize,
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
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
      bool test = (tH != tD);
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
      if(test)
      {
        std::cerr << "ERROR, verifyEvents, event time mismatch for neuron " << nH << ": " 
          << tH << " != " << tD << " ("; PRINT_HEX(4, tH); std::cerr << " != "; PRINT_HEX(4, tD);
          std::cerr << ")" << std::endl;
        result = 0;
        break;
      }
      
      if(result == 0){break;}
      
      /*Verify event weight*/
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
      test = (wH != wD);
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
      if(test)
      {
        /*Count events with equal time*/
        unsigned int timeCount = 0;
        for(unsigned int i = e; i < entryCount; i++)
        {
          unsigned int nD1 = sortedEvents[entryAddress + i];
          DATA_TYPE tH1 = nrn[nH].in_t[i];
          DATA_TYPE tD1 = *((DATA_TYPE *)(&sortedEvents[entryAddress + i + sortedEventsSize]));
          
          WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
          if((nD1 == nH) && (tD1 == tD) && (tH1 == tH)){timeCount++;}
          else{break;}
          WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
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
              WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
              test = ((wH1 == wD1) && !weightCheck[j]);
              WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
              if(test)
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
#endif



/*
  Main stepper routine for PS method on Izhikevich neuron - runs full step
*/
#if UPDATE_NEURONS_ENABLE_V00
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
  WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
  if(tol != 0){nrTolerance = tol;}
  WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
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
#endif



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
  size_t            sortedEventsSize,
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
    
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
    bool test = (underTestType1 != reference);
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
    
    if(test)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable v for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
    underTestType1 = modelVariables[UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].u;
    
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
    test = (underTestType1 != reference);
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
    
    if(test)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable u for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
    underTestType1 = modelVariables[2*UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].g_ampa;
    
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
    test = (underTestType1 != reference);
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
    
    if(test)
    {
      std::cerr << "ERROR, verifyKernelUpdateNeurons, step " << step 
        << ", mismatch of variable g_ampa for neuron " << i 
        << ": " << underTestType1 << " != " << reference << " ("; PRINT_HEX(4, underTestType1); 
        std::cerr << " != "; PRINT_HEX(4, reference); std::cerr << ")" << std::endl;
      result = 0;
    }
    
    underTestType1 = modelVariables[3*UPDATE_NEURONS_TOTAL_NEURONS+i];
    reference = nrn_ps[i].g_gaba;
    
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
    test = (underTestType1 != reference);
    WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
    
    if(test)
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
        
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
      bool test = (te_ps[spiked_neuron] != spike_time);
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
      
      if(test)
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
    
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_START;
      bool test = (te_ps[n] != 0.0 && spikeCheck[n] == 0);
      WARNING_CONTROL_IGNORE_FLOAT_EQUAL_END;
      
      if(test)
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
    char env[] = "GPU_DUMP_DEVICE_KERNEL=3";
    
#if (SYSTEM_OS == SYSTEM_OS_WINDOWS)
    int putenv_res = _putenv(env);
#elif (SYSTEM_OS == SYSTEM_OS_LINUX)
    int putenv_res = putenv(env);
#endif

    if(putenv_res != 0)
    {
      THROW_SIMEX("Neurosim::main: unable to set env variable " << env << ", return code: " 
        << putenv_res);
    }

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
