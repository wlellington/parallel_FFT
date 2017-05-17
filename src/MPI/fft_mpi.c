/* Wesley Ellington
 * Math 4370
 * Semester Project
 * Serial Fast Fourier Transform */

#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "vec2d_b.h"
#include "fft_mpi.h"
#include "cvec1d.h"
#include "mpi.h"

#define pi 3.14159265358979323846264338

/* Implementation of FFT alogrithm and helper functions
 * for 1 dimensional arrays. This set of functions
 * relies on the vec2d_b struct using indexing [i*n + j] notation 
 * for casting to complex number vector*/

/* NOTE: all one dimensional vectors are column formatted. 
 * Using the vec2d_b class, this would be of m by 1 formatting */

/* ============== MPI FFT ============== */

/* Control flow for all procs, all params are valid for all procs minus invec
 *  *  this should be null for all but root */
cvec1d* fft_mpi(cvec1d* invec, cvec1d* outvec, int dir, int myid, int nprocs, int n_total,
        	int m_total, int n_loc, int m_loc){

	/* declarations */
	int i, j, k, ierr, sector, segment, allsize;

	double complex * clone;
	double complex * sbuf;
	double complex * rbuf;
	double complex * signals;
	double complex tmpswp;
	
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Request send, recv;
	MPI_Status status;

	/* allocate send and receive buffers for all procs */
	sbuf = (double complex * ) malloc(n_loc * m_total * sizeof(double complex));
	rbuf = (double complex * ) malloc(n_loc * m_total * sizeof(double complex));
	signals = (double complex *) malloc(n_loc * m_total * sizeof(double complex));

	/* set problem size information for message packaging*/
	sector = n_total / nprocs;
	segment = sector * sector;
	allsize = (n_total * n_total) / (nprocs * nprocs); 

	/* copy data into holding buffer on root */
	if(myid == 0){
		clone = (double complex * ) malloc(n_total * m_total * sizeof(double complex));
		for(i = 0; i < m_total * n_total; i++){
			clone[i] = invec->data[i];
		}
	}
	
	/* for each process, scatter it subsection of data*/
	ierr = MPI_Scatter(clone, n_loc * m_total, MPI_C_DOUBLE_COMPLEX,
			rbuf, n_loc * m_total, MPI_C_DOUBLE_COMPLEX, 0, comm);
	
	err_check(ierr, myid, "Scatter issue");
	
	/* calculate fft for each "row" given by root */
	fft_single(signals, dir, m_total, n_loc);	

	/* copy values form send into working buffer */	
	for(i = 0; i < n_loc * m_total; i++){
		signals[i] = rbuf[i];
	}

	/* create send buffer using reorder over signal */
	send_reorder(signals, sbuf, n_total, nprocs, n_loc);
	
	/* all to all data, receive from all procs */
	ierr = MPI_Alltoall(sbuf, allsize, MPI_C_DOUBLE_COMPLEX,
				rbuf, allsize, MPI_C_DOUBLE_COMPLEX, comm);

	err_check(ierr, myid, "All to All Send or Recieve Issue!");

	/* re order data from send */
	receive_reorder(signals, rbuf, n_total, sector, segment); 	

	/* run fft over each new col */
	fft_single(signals, dir, m_total, n_loc);
	
	/* reassign values into send buffer for gather */
	for(i = 0; i < n_loc * m_total; i++){
		sbuf[i] = signals[i];
	}

	/* prepare buffer for gather on root*/
	if(myid == 0){
		free(rbuf);
		rbuf = (double complex *) malloc(n_total * m_total * sizeof(double complex));
	}	

	/* send final values back to root */
	ierr = MPI_Gather(sbuf, m_total * n_loc, MPI_C_DOUBLE_COMPLEX, 
			rbuf, m_total * n_loc, MPI_C_DOUBLE_COMPLEX, 0, comm);
	
	err_check(ierr, myid, "Gather issue!");

	/* Transpose received matrix */	
	if(myid == 0){
		for(i = 0; i < n_total; i++){
			for(j = 0; j < m_total; j++){
				outvec->data[i * m_total + j] = rbuf[j*n_total + i];
			}
		}
	}

	/* free buffers */
	free(sbuf);
	free(rbuf);
	free(signals);
	
	/* free holding buffer on root */
	if(myid == 0){
		free(clone);
	}

	/* wait for all procs before return */
	MPI_Barrier(comm);

	return outvec;
}

/* process for calculating multiple ffts over several rows */
double complex * fft_single(double complex * input, int dir, int length, int n_lines){
	int g, i, j, k, n, m, linescl;
	double complex Wm, W, t, u;

	/* run reverse order copy over elements of array */
	flat_reverse_copy(input, length, n_lines);

	/* perform fft*/
	n = length;

	/* iterate over rows of signal group */
	for(g = 0; g < n_lines; g++){
		/* sets offset for operataions on each line */
		linescl = g * length;
		/* Divide over lg(n) subtrees */
		for(i = 0; i < log2(n) + 1; i++){
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
				for(k = j; k < n - 1; k += m){
					/* Calculate butterfly for each branch */
					t = W * input[linescl + k + m/2];
					u = input[linescl + k];
					input[linescl + k] = u + t;
					input[linescl + k + m/2] = u - t;
				}
				/* Raise power of Root */
				W = W * Wm;
			}
		}
	}
	
	/* if running inverse function, scale down */
	if(dir){
		for(i = 0; i < length * n_lines; i++){
			input[i] = input[i] / (double)length;	
		}
	}

	/* Return complex vector after full transform */
	return input;
}

double complex * flat_reverse_copy(double complex * input, int length, int n_lines){
	/* for each index that would exist in one row */	
	int i, j, swaploc;

	/* allocate temp destination array */
	double complex * newdata;
	newdata = (double complex * ) malloc(length * n_lines * sizeof(double complex));

	/* copy data to tmp array for reorder */
	for(i  = 0; i < length * n_lines; i++){
		newdata[i] = input[i];
	}
	
	/* for each item in one row */	
	for(i = 0; i < length; i++){
		/* calculate swap location */
		swaploc = reverseBits(i, log2(length));
		/* apply operation to all rows of subsection of matrix */
		for(j = 0; j < n_lines; j++){
			/* place data in final desination */
			input[j * length + swaploc] = newdata[j * length + i]; 
			//input[swaploc] = newdata[i]; 
		}
	}

	/* free temp reorder memory */
	free(newdata);
	
	return input;
}

/* re orders data on send for all to all */
complex * send_reorder(double complex * to_send, double complex * sbuf, 
			int length, int nprocs, int n_lines){
	int i, j, k, counter, chunksize;

	chunksize = n_lines * n_lines;
	counter = 0;
	for(k = 0; k < nprocs; k++){
		for( j = 0; j < n_lines; j++) {
			for(i = 0; i < n_lines; i++){
				/* reorder data in small squares of total stripe */
				sbuf[counter] = to_send[(k * n_lines) + (j * length) + i];
				counter ++;
			}
		}
	}

	return sbuf;
}


/* re orders data on receive from all to all */
complex * receive_reorder(double complex * from_receive, double complex * rbuf, 
			int length, int sector, int segement){
	int i, j, nprocs, counter;

	counter = 0;
	for(i = 0; i < length; i ++){
		for(j = 0; j < sector; j++){
			/* reorder all data received in all to all in column wise fashion */
			from_receive[j * length + i] = rbuf[counter];
			counter ++;		
		}
	}

	return from_receive;
}

/* Create onesided power spectrum from FFT result */
vec2d_b * powerSpec(cvec1d * invec);

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
