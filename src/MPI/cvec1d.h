/* Wesley Ellington
 * Math 4370
 * Serial FFT implementation
 * Complex number array */

#ifndef __CVEC1D
#define __CVEC1D

#include "vec2d_b.h"
#include <complex.h>
#include <stdio.h>
#include <string.h>

/* ============ Complex Array ================== */
typedef struct _cvec1d{
	/* size Variable */
	long int n;

	/* Size variable used only in pseudo-2D representations */
	long int m;

	/* data array */
	complex * data;
} cvec1d;

/* Base constructor */
cvec1d * cvec1dNew (long int size);

/* Destructor */
void cvec1dDestroy(cvec1d * x);

/* Cast funciton for moving real components into usable space */
cvec1d * cplxCast(vec2d_b * invec);

/* Cast function for 2D vectors */
cvec1d * cplxCast2D(vec2d_b * invec);

/* Output to file */
int cvec1dWriteFile(cvec1d* x, const char* outfile);

/* Read input from file*/
cvec1d* cvec1dReadFile(const char * infile, int n);
#endif
