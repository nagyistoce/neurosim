/*Integration routines for the Parker-Sochacki, Runge-Kutta & Bulirsch-Stoer methods*/

#include "Kernel_Primitives.h"



/**
* TODO: This source was never tested. Add description
*/
void segmentedBufferedLoad
(
  uint      segmentCount,
  uint      dataStride,
  uint      *dataSegmentCounts,
  uint      *dataLoadSpace,
  uint      *dataStoreSpace,
){
  __local uint transactionBufferUtilization[wg_size_wfs];
  uint wi_id = get_local_id(0);
  
  uint p = 0, i = 0, dataItemId[2], dataLoadIndex[2], dataStoreIndex[2],
    offsetLoad = wi_id%PRIMITIVES_WF_SIZE_WI, offsetStore = offsetLoad;
  int c = 0;
  
  /*WF: iterate over segments*/
  while(p < segmentCount)
  {
    /*WI: the offset falls into current count*/
    if(dataSegmentCounts[p] > offsetLoad)
    {
      /*WI: buffer count index and the offset*/
      dataItemId[i] = p;
      dataLoadIndex[i] = offsetLoad;
      dataStoreIndex[i] = offsetStore;
      i = !(i); c++;
      offsetLoad += PRIMITIVES_WF_SIZE_WI;
      offsetStore += PRIMITIVES_WF_SIZE_WI;
      atomic_inc(&transactionBufferUtilization[wi_id/PRIMITIVES_WF_SIZE_WI]);
    }
    /*WI: the offset is out of bounds of the current count*/
    else
    {
      /*WI: incrment the count pointer and adjust the offset*/
      offsetLoad -= dataSegmentCounts[p];
      p++;
    }
    
    /*Issue access transaction when the buffer has a full WF*/
    if(transactionBufferUtilization[wi_id/PRIMITIVES_WF_SIZE_WI] >= PRIMITIVES_WF_SIZE_WI)
    {
#if (PRIMITIVES_ERROR_TRACK_ENABLE)
      if(c == 0 || c > 2)
      {
      
      }
#endif
      /*Always store the oldest value. Use logical XOR not to calculate the pointer.*/
      c--; cl_uint a = !(!i != !c);
      
      dataStoreSpace[dataStoreIndex[a]] = 
        dataLoadSpace[dataItemId[a]*dataStride + dataLoadIndex[a]];
        
      if(wi_id%PRIMITIVES_WF_SIZE_WI == 0)
      {
        transactionBufferUtilization[wi_id/PRIMITIVES_WF_SIZE_WI] -= PRIMITIVES_WF_SIZE_WI;
      }
    }
  }
  
  /*Always store the oldest value. Use logical XOR not to calculate the pointer.*/
  c--; cl_uint a = !(!i != !c);
  
  /*Issue access transaction(s) for whatever left in the buffer*/
  while(c >= 0)
  {
    dataStoreSpace[dataStoreIndex[a]] = 
      dataLoadSpace[dataItemId[a]*dataStride + dataLoadIndex[a]];
    a = !(a);
  }
}
