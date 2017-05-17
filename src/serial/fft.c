/* Wesley Ellington
 * Math 4370
 * Semester Project
 * Serial Fast Fourier Transform */

#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "vec2d_b.h"
#include "fft.h"
#include "cvec1d.h"

#define pi 3.14159265358979323846264338

/* Implementation of FFT alogrithm and helper functions
 * for 1 dimensional arrays. This set of functions
 * relies on the vec2d_b struct using indexing [i*n + j] notation 
 * for casting to complex number vector*/

/* NOTE: all one dimensional vectors are column formatted. 
 * Using the vec2d_b class, this would be of m by 1 formatting */

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
		x->data[i] = x->data[i]/max;
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
cvec1d* fft_1dFunc(cvec1d * x, int dir){
	int n = x->n;
	int i, j, k;
	/* Reorder data into bit reversed order for forward tranfrom */
	x = reverseCopy(x);

	/* Divide over lg(n) subtrees */
	for(i = 0; i < log2(n) + 1; i++){
		int m = pow(2, i);
 
		/* Set principle root of unity */
		complex Wm =  cexp(2*pi*-_Complex_I/(complex)m);
		complex W = 1;
		
		/* Set transform direction */
		if(dir){
			Wm = cexp(2*pi*_Complex_I/(complex)m) ;
		}

		/* For each value associated with tree span */
		for(j = 0; j < m/2; j ++){
			for(k = j; k < n - 1; k += m){
				/* Calculate butterfly for each branch */
				complex t = W * x->data[k + m/2];
				complex u = x->data[k];
				x->data[k] = u + t;
				x->data[k + m/2] = u - t;
			}
			/* Raise power of Root */
			W = W * Wm;
		}
	}
		
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
	for(i = 0; i < x->n; i++){
		x->data[i] = x->data[i] * scaleDown;
	}

	return x;
}

/* ======================= Parallel Implementations =========================*/

/* Parallelizable 1-2D FFT wrapper function*/
cvec1d* fft(cvec1d* signal){
	/* Verify proper size and non broken data */
	if(fftVerif(signal) != 0){
		fprintf(stderr, "FFT Verification Failed, aborting!\n");
		return NULL;
	}

	/* run 1-2D forward fft on signal */ 
		
	/* Reorganize data accross row*/
	revCpy(signal, 0);
	/* Run FFT on each row */
	fftFunc(signal, 0, 0);	
	/* Reorganize data accross col*/
	revCpy(signal, 1);
	/* Run FFT on each col */
	fftFunc(signal, 0, 1);

	return signal;
}

/* 2D inverse FFT with verification */
cvec1d* ifft(cvec1d* coeffs){
	/* Verify proper size and non broken data */
	if(fftVerif(coeffs) != 0){
		fprintf(stderr, "IFFT Verification Failed, aborting!\n");
		return NULL;
	}
	
	/* run 1-2D inverse fft on signal */ 
		
	/* Reorganize data accross row*/
	revCpy(coeffs, 0);
	/* Run FFT on each row */
	fftFunc(coeffs, 1, 0);	
	/* Reorganize data accross col*/
	revCpy(coeffs, 1);
	/* Run FFT on each col */
	fftFunc(coeffs, 1, 1);

	return coeffs;
}

/* 2D verification Function*/
int fftVerif(cvec1d* input){
	int m,n;
	/* Check for null data pointer */
	if(input->data == NULL){
		fprintf(stderr, "Null data pointer!\n");
		return 1;
	}
	
	/* Ensure m and n are powers of Two */
	if(!((input->m != 0) && !(input->m & (input->m - 1)))){
		fprintf(stderr, "M is not a power of two!\n");
		return 1;
	}
	if(!((input->n != 0) && !(input->n & (input->n - 1)))){
		fprintf(stderr, "N is not a power of two!\n");
		return 1;
	}
	
	/* No errors */
	return 0;
}

/* 2D iterative FFT using inplace solve*/
/* dir sets direction of transform from forward to inverse */
/* rowCol sets iteration by row or column*/
cvec1d* fftFunc(cvec1d* x, int dir, int rowCol){
	int row, col, i, j, k;
	int xsize, ysize;
	double scaleDown = 1.0;	
	/* if iterating over rows */
	if(rowCol == 0){
		/* row update size */
		xsize = x->n;
		for(row = 0; row < x->m; row++){	
			for(i = 0; i < log2(xsize) + 1; i++){
				int m = pow(2, i);

				/* Set principle root of unity */
				complex Wm =  cexp(2*pi*-_Complex_I/(complex)m);
				complex W = 1;
				
				/* Set transform direction */
				if(dir){
					Wm = cexp(2*pi*_Complex_I/(complex)m) ;
				}

				/* For each value associated with tree span */
				for(j = 0; j < m/2; j ++){
					for(k = j; k < xsize - 1; k += m){
						/* Calculate butterfly for each branch */
						complex t = W * x->data[(row*xsize) + (k + m/2)];
						complex u = x->data[(row*xsize) + k];
						x->data[(row*xsize) + k] = u + t;
						x->data[(row*xsize) + (k + m/2)] = u - t;
					}
					/* Raise power of Root */
					W = W * Wm;
				}
			}
			
			/* if running inverse, scale down*/
			if(dir){
				for(j = 0; j < x->n; j++){
					x->data[row*xsize + j] = x->data[row*xsize + j] * (1.0/xsize);
				}
			}
		}
	}
	
	/* for iterating over cols */
	else{
		/* row update size */
		ysize = x->m;
		for(col = 0; col < x->m; col++){	
			for(i = 0; i < log2(ysize) + 1; i++){
				int m = pow(2, i);

				/* Set principle root of unity */
				complex Wm =  cexp(2*pi*-_Complex_I/(complex)m);
				complex W = 1;
				
				/* Set transform direction */
				if(dir){
					Wm = cexp(2*pi*_Complex_I/(complex)m) ;
				}

				/* For each value associated with tree span */
				for(j = 0; j < m/2; j ++){
					for(k = j; k < ysize; k += m){
						/* Calculate butterfly for each branch */
						complex t = W * x->data[col + (k + m/2)*ysize];
						complex u = x->data[col + k*ysize];
						x->data[col + k*ysize] = u + t;
						x->data[col + (k + m/2)*ysize] = u - t;
					}
					/* Raise power of Root */
					W = W * Wm;
				}
			}
			
			/* if running inverse, scale down*/
			if(dir){
				for(j = 0; j < x->n; j++){
					x->data[col + j*ysize] = x->data[col + j*ysize] * (1.0/ysize);
				}
			}
		}


	}
	return x;
}

/* In place reverse copy */
/* rowCol = 0, reverse each row, rowCol = 1, reverse col */
void revCpy(cvec1d* input, int rowCol){
	
	int i, j, m, n, swapInt, step;
	double complex tmp;
	m = input->m;
	n = input->n;

	/* if runnng reverse on rows */
	if (rowCol == 0){
		/* For the first half of the indicies int a row*/
		for(j = 0; j < m/2; j++){
			/* Calculate swap buddies */
			swapInt = reverseBits(j,log2(m));
			if(swapInt == j){
				/* skip iteration of same location*/
				continue;
			}
			/* For each row in Vect*/
			for(i = 0; i < n; i++){
				/* Swap pairs */
				tmp = input->data[i*n + j];
				input->data[i*n + j] = input->data[i*n + swapInt];
				input->data[i*n + swapInt] = tmp;
			}
		}		
	}

	/* if running reverse on Cols */
	else {
		/* For the first half of the indicies in a col */
		for(j = 0; j < n/2; j++){
			/* Calculate swap buddies */
			swapInt = reverseBits(j,log2(n));
			if(swapInt == j){
				/* skip iteration of same location*/
				continue;
			}
			/* For each col in vect */
			for(i = 0; i < m; i++){
				/* Swap pairs */
				tmp = input->data[j*n + i];
				input->data[j*n + i] = input->data[swapInt*n + i];
				input->data[swapInt*n + i] = tmp;
			}
		}		
	}	
}

