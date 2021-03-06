/* Wesley Ellington
 * Math 4370
 * Semester Project
 * Serial FFT class */

#ifndef FFT_DEFINED__
#define FFT_DEFINED__

/*STL inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

/* User Defined Dependancies*/
#include "vec2d_b.h"
#include "cvec1d.h"

/* Boolean values */
#define True 1
#define False 0


/* Create onesided power spectrum from FFT result */
vec2d_b * powerSpec(cvec1d * invec);

/* ============ One Dimensional FFT ============ */

/* 1D Base call helper function */
/* This function should be called to run any program using 
 * FFT, as it first runs check methods to verify input, then 
 * applies tranform */
cvec1d* fft_1d(cvec1d* input);

/* Experimental 2d fft algorithm */
cvec1d* fft_2d(cvec1d* input);

/* Experimental 2d ifft algorithm */
cvec1d* ifft_2d(cvec1d* input);

/* Function to reorder values for initial comp */
cvec1d* reverseCopy(cvec1d* inVec);

/* Reverse int by bit order */
int reverseBits(int num, int size);

/* 1D input verification function for complex inputs */
int fft_1dVerif(cvec1d* input);

/* 1D Transform function */
/* Actual FFT algorithm, should not be called directly*/
cvec1d* fft_1dFunc(cvec1d * x, int dir);

/* inverse FFT transform using posative root of unity*/
cvec1d* ifft_1d(cvec1d* input);

/* Parallelizable 1-2D FFT wrapper function*/
cvec1d* fft(cvec1d* signal);

/* 2D verification Function*/
int fftVerif(cvec1d* input);

/* 2D iterative FFT using inplace solve*/
cvec1d* fftFunc(cvec1d* input, int dir, int rowCol);

/* 2D inverse FFT with verification */
cvec1d* ifft(cvec1d* input);

/* Whole Vector reverse copy */
void  revCpy(cvec1d* input, int rowCol);
#endif
