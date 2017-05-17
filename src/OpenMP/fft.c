/* Wesley Ellington
 *  * Math 4370
 *   * Semester Project
 *    * Serial Fast Fourier Transform */

#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "vec2d_b.h"
#include "fft.h"
#include "cvec1d.h"
#include "get_time.h"

#define pi 3.14159265358979323846264338

/* Implementation of FFT alogrithm and helper functions
 *  * for 1 dimensional arrays. This set of functions
 *   * relies on the vec2d_b struct using indexing [i*n + j] notation 
 *    * for casting to complex number vector*/

/* NOTE: all one dimensional vectors are column formatted. 
 *  * Using the vec2d_b class, this would be of m by 1 formatting */

/* Create one sided, normalized power spectrum */
vec2d_b * powerSpec(cvec1d * invec){
	vec2d_b * x = vec2d_bNew(invec->n, 1);
	int i,j;
	double max = 0;
	complex powerVal;
	
	/* iterate over values in complex array */
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			
			/* calculate memeber wise conjugate * complex */
			powerVal = invec->data[i*x->n + j] *
				conj(invec->data[i*x->n + j]);
			
			/* assign to power specturm */
			x->data[i*x->n + j] = (double)creal(powerVal);
		}
	}

	/* Normalize to maximum value of spectrum */	
	max = vec2d_bMax(x);

	for(i = 0; i < x->m; i++){
		//x->data[i] = x->data[i]/max;
	}

	return x;
}

/* ================== FFT Implementation =============== */

/* Basic wrapper function for verification and execution */
cvec1d* fft_1d(cvec1d* timespec){
	/* Run verificaiton on vectors */
	if(fft_1dVerif(timespec) != 0){
		fprintf(stderr, "1D FFT Err, Verification failed!\n");
		return NULL;
	}
	
	/* Run transform */
	return fft_1dFunc(timespec, 0);
}

/* Reorder values in vector to bit reversed initial positions */
cvec1d* reverseCopy(cvec1d* inVec){
	int i =0;

	cvec1d* result = cvec1dNew(inVec->n);
	result->n = inVec->n;
	for(i = 0; i < inVec->n; i++){
		result->data[reverseBits(i, log2(inVec->n))] = inVec->data[i];
	}
	return result;
}

/* Reverse bits of num */
/* numSize corresponds to the number of bits in largest value in set */
int reverseBits(int num, int numSize){
	int reverse = 0;
     	int mask = 1;
	int remain = numSize - 1;
	while(remain){
		reverse |= num & 1; 
		
		remain --;
		num >>= 1;
		reverse <<= 1;       
	}
	reverse |= num & 1;
	return reverse;
}

/* Verification for Complex Vector */
int fft_1dVerif(cvec1d * input){	
	int isPower = False;
	int power = 0;
	int powerSize = 0;

	/* Check if matrix is power of two size */
	while(input->n > powerSize){
		powerSize = pow(2, power);
		if(powerSize == input->n){
			isPower = True;
			break;
		}
		power ++;
	}

	/* Return error if vector is not power of two*/
	if(isPower == False){
		fprintf(stderr,"1D FFT Err, Vectors not of size 2^n!\n");
		return 1;
	}

	return 0;
}

/* Iterative FFT implementation */
/* dir = 0 implies forward, 1 is reverse*/
/* for OpenMP implementation, barriers are needed after each loop to ensure
 * sync from one level of butterfly to the next */
cvec1d* fft_1dFunc(cvec1d * x, int dir){
	int i, j, k, m, n;
	complex Wm, W, u, t;
	double stime, etime;	
	#pragma omp parallel default(shared) private(i, j, k, u, t, W, Wm, m)
	{
	
	#pragma omp single
	{
	n = x->n;
	/* Reorder data into bit reversed order for forward tranfrom */
	x = reverseCopy(x);
	} /* end single */
	/* Divide over lg(n) subtrees */

	for(i = 0; i < log2(n) + 1; i++){
		/* set initial values */
		m = 1 << i;

		/* Set principle root of unity */
		Wm =  cexp(2*pi*-_Complex_I/(complex)m);
		W = 1;
		
		/* Set transform direction */
		if(dir){
			Wm = cexp(2*pi*_Complex_I/(complex)m) ;
		}
		
		/* For each value associated with tree span */
		for(j = 0; j < m/2; j ++){
				
			#pragma omp for schedule(static, 1)
			for(k = j; k < n -1; k += m){
				/* Calculate butterfly for each branch */
				t = W * x->data[k + m/2];
				u = x->data[k];
				x->data[k] = u + t;
				x->data[k + m/2] = u - t;
			}

			/* Raise power of Root */
			W = W * Wm;
		}

	}
	

	} /* end parallel region */	
	/* Return complex vector after full transform */
	return x;
}

/* Inverse FFT using forward FFT with posative principle root swapping */
cvec1d* ifft_1d(cvec1d* x){
	int i;

	/* Calculate scale down factor*/
	double scaleDown = 1.0/(double)(x->n);
			
	/* Run forward FFT */
	x = fft_1dFunc(x, 1);

	/* Scale down points */
	#pragma omp parallel for
	for(i = 0; i < x->n; i++){
		x->data[i] = x->data[i] * scaleDown;
	}

	return x;
}

/* 2D FFT */
cvec1d* fft_2d(cvec1d* x){
	int i, j;
	cvec1d* temp;
	
	/* verify valid inputs */
	if (x->m != x->n){
		fprintf(stderr, "2D FFT not runnable on non square image\n");
		return NULL;
	}
	
	#pragma omp parallel default(shared) private(i, j, temp)

	/* For each row of vector */
	#pragma omp for schedule(dynamic)
	for(i = 0; i < x->m; i++){

		/* Create vector for row */
		temp = cvec1dNew(x->n);
		for(j = 0; j < x->n; j++){
			temp->data[j] = x->data[i*x->n + j];
		}
	
		/* Run fft on  row  */
		temp = fft_1d(temp);
	
		/* Assign back to master signal */
		for(j = 0; j < x->n; j++){
			x->data[i*x->n + j] = temp->data[j];
		}
		cvec1dDestroy(temp);
	}
	
	#pragma omp barrier

	/* For each col of vector */
	#pragma omp for schedule(dynamic)
	for(j = 0; j < x->n; j++){

		/* Create vector for col */
		temp = cvec1dNew(x->n);
		for(i = 0; i < x->m; i++){
			temp->data[i] = x->data[i*x->n + j];
		}
	
		/* Run fft on col */
		temp = fft_1d(temp);
	
		/* Assign back to master signal */
		for(i = 0; i < x->m; i++){
			x->data[i*x->n + j] = temp->data[i];
		}
		cvec1dDestroy(temp);
	}

	#pragma omp barrier
	return x;
}

/* Inverse FFt for 2d Vector */
cvec1d* ifft_2d(cvec1d* x){

	int i, j;
	cvec1d* temp;
	/* verify valid inputs */

	if (x->m != x->n){
		fprintf(stderr, "2D FFT not runnable on non square image\n");
		return NULL;
	}

	#pragma omp parallel default(shared) private(i, j, temp)
	
	/* For each col of vector */
	#pragma omp for schedule(dynamic)
	for(j = 0; j < x->n; j++){

		/* Create vector for col */
		temp = cvec1dNew(x->n);
		for(i = 0; i < x->m; i++){
			temp->data[i] = x->data[i*x->n + j];
		}
	
		/* Run fft on col */
		temp = ifft_1d(temp);
	
		/* Assign back to master signal */
		for(i = 0; i < x->m; i++){
			x->data[i*x->n + j] = temp->data[i];
		}
		cvec1dDestroy(temp);
	}

	#pragma omp barrier

	/* For each row of vector */
	#pragma omp for schedule(dynamic)
	for(i = 0; i < x->m; i++){

		/* Create vector for row */
		temp = cvec1dNew(x->n);
		for(j = 0; j < x->n; j++){
			temp->data[j] = x->data[i*x->n + j];
		}
	
		/* Run fft on  row  */
		temp = ifft_1d(temp);
	
		/* Assign back to master signal */
		for(j = 0; j < x->n; j++){
			x->data[i*x->n + j] = temp->data[j];
		}
		/* Clean up temp variables */
		cvec1dDestroy(temp);
	}
	
	#pragma omp barrier

	return x;
}
