/*Header file to accompany iz_util.c*/
/*Written by Dr Robert Stewart for Stewart & Bair, 2009*/

#ifndef INC_IZ_UTIL_H
#define INC_IZ_UTIL_H

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



typedef struct 
{
	double E; 											/*Electrical elastance 1/C*/
	double vr;											/*Resting mebrane potential*/
	double k;												/*Scaling constant*/
	double l;												/*-k*vt*/
	double a;												/*Recovery variable rate constant*/
	double b;												/*Scaling constant*/
	double v_peak;									/*Peak voltage during spike*/
	double v_reset;									/*Post-spike reset potential*/
	double u_step;									/*Post-spike recovery variable step*/
	double I;												/*Input current*/
	double g_ampa;									/*AMPA conductance*/
	double g_gaba;									/*GABA_A conductance*/
	double E_ampa;									/*AMPA reversal potential*/
	double E_gaba;									/*GABA reversal potential*/
	double v;												/*Membrane voltage*/
	double u;												/*Recovery variable*/
	int n_in;												/*Number of synaptic events in input buffers*/
	double in_t[REFERENCE_EVENT_QUEUE_SIZE];						/*Time input buffer*/
	float in_w[REFERENCE_EVENT_QUEUE_SIZE];							/*Weight input buffer*/
} neuron_iz;

typedef struct 
{
	DATA_TYPE E; 											/*Electrical elastance 1/C*/
	DATA_TYPE vr;											/*Resting mebrane potential*/
	DATA_TYPE k;											/*Scaling constant*/
	DATA_TYPE l;											/*-k*vt*/
	DATA_TYPE a;											/*Recovery variable rate constant*/
	DATA_TYPE b;											/*Scaling constant*/
	DATA_TYPE v_peak;									/*Peak voltage during spike*/
	DATA_TYPE v_reset;								/*Post-spike reset potential*/
	DATA_TYPE u_step;									/*Post-spike recovery variable step*/
	DATA_TYPE I;											/*Input current*/
	DATA_TYPE g_ampa;									/*AMPA conductance*/
	DATA_TYPE g_gaba;									/*GABA_A conductance*/
	DATA_TYPE E_ampa;									/*AMPA reversal potential*/
	DATA_TYPE E_gaba;									/*GABA reversal potential*/
	DATA_TYPE v;											/*Membrane voltage*/
	DATA_TYPE u;											/*Recovery variable*/
	int n_in;												/*Number of synaptic events in input buffers*/
	DATA_TYPE in_t[REFERENCE_EVENT_QUEUE_SIZE];						/*Time input buffer*/
	float in_w[REFERENCE_EVENT_QUEUE_SIZE];							/*Weight input buffer*/
} neuron_iz_ps;

typedef struct 
{
	DATA_TYPE E; 										/*Electrical elastance 1/C*/
	DATA_TYPE vr;										/*Resting mebrane potential*/
	DATA_TYPE k;											/*Scaling constant*/
	DATA_TYPE l;											/*-k*vt*/
	DATA_TYPE a;											/*Recovery variable rate constant*/
	DATA_TYPE b;											/*Scaling constant*/
	DATA_TYPE v_peak;								/*Peak voltage during spike*/
	DATA_TYPE v_reset;								/*Post-spike reset potential*/
	DATA_TYPE u_step;								/*Post-spike recovery variable step*/
	DATA_TYPE I;											/*Input current*/
	DATA_TYPE g_ampa;								/*AMPA conductance*/
	DATA_TYPE g_gaba;								/*GABA_A conductance*/
	DATA_TYPE E_ampa;								/*AMPA reversal potential*/
	DATA_TYPE E_gaba;								/*GABA reversal potential*/
	DATA_TYPE v;											/*Membrane voltage*/
	DATA_TYPE u;											/*Recovery variable*/
	int n_in;												/*Number of synaptic events in input buffers*/
	DATA_TYPE in_t[REFERENCE_EVENT_QUEUE_SIZE];					/*Time input buffer*/
	float in_w[REFERENCE_EVENT_QUEUE_SIZE];							/*Weight input buffer*/
} neuron_iz_gpu;

typedef struct 
{
	unsigned int n; 								/*Target or source neuron*/
	float w; 												/*Weight*/
	float d;						            /*Delay interval normalized to dt*/
} synapse;

void iz_derivs(double *, double *, double *);
void iz_first(DATA_TYPE **, DATA_TYPE **, DATA_TYPE *);
void iz_iter(DATA_TYPE **, DATA_TYPE **, DATA_TYPE *, int);

#endif /* INC_IZ_UTIL_H */
