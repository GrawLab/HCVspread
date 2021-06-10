// ######################################################################
// PROJECT: SPATIAL SIMULATION OF VIRAL SPREAD WITHIN A 2D LAYER OF CELLS
//
// random number generator
// ###############################################################

// ###############################################################
// MAIN HEADER FILES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

// ###############################################################
// RANDOM NUMBER FUNCTIONS

void randbin(double *output, double p){
	GetRNGstate();
	double b=rbinom(1,p);
	PutRNGstate();
	*output=b;
}

void randuniform(double *output, double min, double max){
	GetRNGstate();
	double b=runif(min,max);
	PutRNGstate();
	*output=b;
}

void randpoisson(double *output, double a){
	GetRNGstate();
	double b=rpois(a);
	PutRNGstate();
	*output=b;
}

void randgamma(double *output, double lT, double scale){
	GetRNGstate();
	double b=rgamma(lT,scale);
	PutRNGstate();
	*output=b;
}

void densgamma(double *output, double vir, double shape, double scale){
	GetRNGstate();
	double b=dgamma(vir,shape,scale,0);
	PutRNGstate();
	*output=b;
}

void randbeta(double *output, double shape1, double shape2){
	GetRNGstate();
	double b=rbeta(shape1,shape2);
	PutRNGstate();
	*output=b;
}

void randnorm(double *output, double mean, double sd){
    GetRNGstate();
    double b=rnorm(mean,sd);
    PutRNGstate();
    *output=b;
}

// ################################################################
