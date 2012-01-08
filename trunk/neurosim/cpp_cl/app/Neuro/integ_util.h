
#ifndef INC_INTEG_UTIL_H
#define INC_INTEG_UTIL_H

#include "Definitions.h"
/*Models: */
#include "iz_util.h"
#include "integ_util.h"
#include <CL/cl.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <SDKCommon.hpp>
#include <SDKApplication.hpp>
#include <SDKFile.hpp>
#include <SDKBitMap.hpp>



extern void ps_update(DATA_TYPE **, int, int, DATA_TYPE, DATA_TYPE *);

extern int ps_step
(
#if (UPDATE_NEURONS_TOLERANCE_MODE != 0)
  DATA_TYPE *,
#endif
  DATA_TYPE **,
  DATA_TYPE **,
  DATA_TYPE *,
  DATA_TYPE *,
  DATA_TYPE *,
  void (*)(DATA_TYPE **,DATA_TYPE **,DATA_TYPE *), 
	void (*)(DATA_TYPE **,DATA_TYPE **,DATA_TYPE *,int),
  int, 
  int, 
  int,
  int &
);
  
extern int ps_step0(double **,double **,double *,double *,double *,
  void (*)(double **,double **,double *), 
	void (*)(double **,double **,double *,int), int, int, int, int);

extern void cauchy_prod(int, double *, double, double *, double, double *);

extern void cauchy_prod_basic(int, double *, double *, double *);

extern void cauchy_prod_signed(int, double *, double, int, double *, double, 
  int, double *);

extern void cauchy_prod_gain(int, double *, double, double, double *, double, 
  double, double *);

extern void series_div(int, double, double *, double, double *, double);

extern void series_pow(int, double *, double, double *, double, double);

extern void rk_step(double *, double *, double *, int, double *, double,
	void (*)(double *, double *, double *));
  
extern void mmid(double *, double *, int, double, int, double *, double *,
	void (*)(double *, double *, double *));

extern void pzextr(int,double,double *,double *,double *,int,double *,
  double **);

extern void rzextr(int,double,double *,double *,double *,int,double *,
  double **);

extern int bs_step(double *, double *, double *, int, double, double *, 
  double *, void (*)(double *, double *, double *));
  
#endif /* INC_INTEG_UTIL_H */
