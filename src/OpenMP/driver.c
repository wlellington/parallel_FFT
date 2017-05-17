/* Wesley Ellington
 * Math 4370
 * Semester Project
 * Serial Fast Fourier Transform */

/* This is a base driver for testing and running the FFT 
 * Algorithm based on the Cooley Tukey design. It has both 
 * a one dimensional and two dimensional implementation. */ 

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vec2d_b.h"
#include "fft.h"
#include "get_time.h"


/* main method */
int main(int argc, const char* argv){
	omp_set_nested(0);
	double stime, etime, optime;
	vec2d_b * x = vec2d_bReadFile("samples.txt", 1, 8192 * 256);
	vec2d_b * y = vec2d_bReadFile("2Dcount.txt", 1024, 1024);
	cvec1d * signal2d;
	cvec1d * signal;
	vec2d_b * pwrSpec;
	

	
	printf("Running 1D FFT of size %d\n", x->n);
	signal = cplxCast(x);

	cvec1dWriteFile(signal, "signal.txt");
	stime = get_time();
	signal = fft_1d(signal);
	etime = get_time();
	optime = etime - stime;
	printf("\t1D FFT Run in: %f sec\n", optime);
	cvec1dWriteFile(signal, "coeffs.txt");
	pwrSpec = powerSpec(signal);
	vec2d_bWriteFile(pwrSpec, "powerSpec.txt");
	stime = get_time(); 
	signal = ifft_1d(signal);
	etime = get_time();
	optime = etime - stime;
	printf("\t1D IFFT Run in: %f sec\n", optime);
	cvec1dWriteFile(signal, "rSignal.txt");


	printf("Running 2D FFT of size %d by %d\n", y->m, y->n);
	signal2d = cplxCast(y);

	cvec1dWriteFile(signal2d, "signal2d.txt");
	stime = get_time();
	signal2d = fft_2d(signal2d);
	etime = get_time();
	optime = etime - stime;
	printf("\t2D FFT Run in: %f sec\n", optime);
	cvec1dWriteFile(signal2d, "coeffs2d.txt");
	stime = get_time();
	signal2d = ifft_2d(signal2d);	
	etime = get_time();
	optime = etime - stime;
	printf("\t2D IFFT Run in: %f sec\n", optime);
	cvec1dWriteFile(signal2d, "rSignal2d.txt");

	return 0;
}
