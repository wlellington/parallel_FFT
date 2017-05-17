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
#include "mpi.h"

#include "vec2d_b.h"
#include "fft_mpi.h"


/* utility function for error checking */
void err_check(int ierr, int myid, const char* err){
	/* if no error exists, do nothing */
	if(ierr == 0){
		return;
	}
	/* if error is detected, print string and kill Comm */
	else{
		fprintf(stderr, "Err on %d : %s\n", myid, err);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/* utility function for checking in from procs */
void check_in(int myid, const char* msg){
	printf("Proc %d says: %s \n", myid, msg);
	fflush(stdout);
}


/* main method */
int main(int argc, char* argv[]){
	/* mpi proc info */
	int ierr, nprocs, myid;
	double stime, ftime;

	/* problem size information */
	int m_total, n_total, n_loc, m_loc;

	/* input file name */
	char * input_file = "2Dcount.txt";

	/* root process structs for IO */
	vec2d_b * real_signal;
	cvec1d * coeffs;
	cvec1d * signal;
	cvec1d * rsignal;

	/* initialize MPI */
	ierr = MPI_Init(&argc, &argv);
	err_check(ierr, 0, "Initialization Err");

	/* set problem total problem size */
	m_total = 4096;
	n_total = 4096;
	

	/* assign process specific vars */
	ierr += MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	ierr += MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	err_check(ierr, myid, "Process declaration Err");
	
	/* set local problem size, since problems are all of size 2^n, we 
 	*  need not deal with edge cases assuming nprocs is multiple of two*/
	n_loc = n_total / nprocs;
	m_loc = m_total / nprocs;

	/* load all problem data into process zero */
	if(myid == 0){
		vec2d_b * real_signal = vec2d_bReadFile(input_file, m_total, n_total);	
		signal = cplxCast(real_signal);
		coeffs = cvec1dNew(m_total * n_total);
		coeffs->m = m_total;
		coeffs->n = n_total;
		rsignal = cvec1dNew(m_total * n_total);
		rsignal->m = m_total;
		rsignal->n = n_total;
	}

	/* main fft control flow */
	/* run and time forward FFT */
	if(myid == 0 ){
		printf("Running FFT on %d by %d with %d processes\n", m_total, n_total, nprocs);
		stime = MPI_Wtime();	
	}

	fft_mpi(signal, coeffs, 0, myid, nprocs, n_total,
        	m_total, n_loc, m_loc);

	if(myid == 0 ){
		ftime = MPI_Wtime();
		printf("\tExecution time of fft: %.16e\n", ftime-stime);
	}

	/* run and time reverse FFT */
	if(myid == 0 ){
		stime = MPI_Wtime();	
	}

	fft_mpi(coeffs, rsignal, 1, myid, nprocs, n_total,
        	m_total, n_loc, m_loc);	
	
	if(myid == 0 ){
		ftime = MPI_Wtime();
		printf("\tExecution time of ifft: %.16e\n", ftime-stime);
	}

	/* output to file if root */
	if(myid == 0){
		cvec1dWriteFile(signal, "signal.txt");
		cvec1dWriteFile(coeffs, "coeffs.txt");
		cvec1dWriteFile(rsignal, "rsignal.txt");
		cvec1dDestroy(coeffs);
		cvec1dDestroy(signal);
		cvec1dDestroy(rsignal);
	}

	ierr += MPI_Finalize();

	return 0;
}
