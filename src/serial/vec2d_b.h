/* Wesley Ellington
 * Math 4370 
 * Homework 2
 * 2 Dimensional Vector Struct 
 * Using linear(* type) Vector memory*/

#ifndef VEC2D_DEFINED__
#define VEC2D_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*Boolean values*/
#define True 1
#define False 0

/* Definition of member variables*/
typedef struct _vec2d_b{
	/*size variables*/
	long int m;
	long int n;
	
	/*value array*/
	double * data;
} vec2d_b;

/* Base constructor, allocates and zeros memory */
vec2d_b* vec2d_bNew(long int mlen, long int nlen);

/* Linear spaced constructor */
vec2d_b* vec2d_bLinspace(double a, double b, 
	long int m, long int n);

/* Randomly populated vector constructor */
vec2d_b* vec2d_bRandom(long int m, long int n);

/* Destructor, deallocates memory */
void vec2d_bDestroy(vec2d_b* inVec);

/* Create vector object from input array and size */
vec2d_b* vec2d_bVectorize(double * data, int size, int mlen, int nlen);

/* I/O functions */
/* Write to stdout*/
int vec2d_bWrite(vec2d_b* x);

/* Write to file */
int vec2d_bWriteFile(vec2d_b* x, const char* outfile);

/* Read in from file */
vec2d_b* vec2d_bReadFile(const char* infile, int m, int n);
 
/* Linear sum function */
/* x = a*y + b*z */
int vec2d_bLinearSum(vec2d_b* x, double a, vec2d_b* y,
		double b, vec2d_b* z );

/* Scale function */
/* x = a*x */
int vec2d_bScale(vec2d_b* x, double a);

/* Vector copy */
/* x = y */
int vec2d_bCopy(vec2d_b* x, vec2d_b* y);

/* Vector Constant */
/* sets all entries to a */
int vec2d_bConstant(vec2d_b* x, double a);

/* Scalar values */
double vec2d_bMin(vec2d_b* x);
double vec2d_bMax(vec2d_b* x);
double vec2d_bDot(vec2d_b* x, vec2d_b* y);
double vec2d_bTwoNorm(vec2d_b* x);
double vec2d_bRmsNorm(vec2d_b* x);
double vec2d_bMaxNorm(vec2d_b* x);

#endif
