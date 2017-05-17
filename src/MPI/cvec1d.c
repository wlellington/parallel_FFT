/* Wesley Ellington
 * Math 4370
 * Serial FFT implementation
 * Complex Vector struct for complex space work
 */

/* Inclusions */
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "vec2d_b.h"
#include "cvec1d.h"

/* =================== Complex Vector ================ */

/* Base constructor */
cvec1d * cvec1dNew(int long size){
	cvec1d * x = (cvec1d*) malloc(sizeof(cvec1d));
	x->data = (complex *) calloc(size,sizeof(complex));
	x->n = size;
	x->m = 1;
	return x;
}

/* Destructor */
void cvec1dDestroy(cvec1d * x){
	free(x->data);
	x->n = 0;
	x->m = 0;
}

/* Cast real components to into  imaginary space */
cvec1d * cplxCast(vec2d_b * invec){
	int i, j;
	
	/* Size variable for new array */
	long int size = invec->m * invec->n;
	
	/* Allocate vector pointer */
	cvec1d * x = (cvec1d * ) malloc(sizeof(cvec1d));
	x->n = invec->n;
	x->m = invec->m;
	
	/* Allocate array space */
	x->data = (complex *) malloc(sizeof(complex) * size);	

	/* Copy data from reals into cplx vec */
	for(i = 0; i < invec->m; i++){
		for(j = 0; j < invec->n; j++){
			x->data[i*invec->n+j] = invec->data[i*invec->n+j] + 0.0*I;
		}
	}
	return x;
}

/* Cast real componenets onto 2D complex vector */
cvec1d * cplxCast2D(vec2d_b * invec){
	int i, j;
	cvec1d * x = (cvec1d*)malloc(sizeof(cvec1d));
	
	/* Set sizing of interal value arrays */
	x->n = invec->n;
	x->m = invec->m;

	/* Allocate array space */
	x->data = (complex* ) malloc(sizeof(complex) * x->n * x->m);

	/* Copy data from real to complex vector */
	for(i = 0; i < invec->m; i++){
		for(j = 0; j < invec->n; j++){
			x->data[i*invec->n+j] = invec->data[i*invec->n+j] + 0.0*I;
		}
	}
	return x;
}

/* Write complex Vector to file*/
int cvec1dWriteFile(cvec1d* x, const char* outfile){
	long int i, j;
	FILE * fptr = NULL;

	/* Quit if vector has no alloced array*/
	if (x->data == NULL){
		fprintf(stderr, "Write to file error, no array\n");
		return 1;
	}

	/* Return if no file string specified */
	if (strlen(outfile) < 1){
		fprintf(stderr, "Write for file error, no name\n");
		return 1;
	}

	/* Open file and verify valid */
	fptr = fopen(outfile, "w");
	if(fptr == NULL){
		fprintf(stderr, "Write to file error, unable to open %s \n", outfile);
		return 1;
	}

	/* Print to file, seperating real and imaginary components with space */
	for(i = 0; i < x->n; i++){
		for(j= 0; j < x->m; j++){
			fprintf(fptr, "(%f %fj) ", creal(x->data[i*x->m+j]), \
				cimag(x->data[i*x->m+j]));
		}
		fprintf(fptr, "\n");
	}

	/* Close file */
	fclose(fptr);
	return 0;
} 

/* Read in imaginary Vector */
cvec1d* cvec1dReadFile(const char * infile, int n){
	
	/* declare vars */
	FILE * input;
	int i;

	/* Create new zeroed out vector of m by n */
	cvec1d* x = cvec1dNew(n);

	/* Check for error conditions */
	if(x == NULL){
		fprintf(stderr,"Read File error, new array null\n");
		return NULL;
	}
	if(x->data == NULL){
		fprintf(stderr, "Read File error, array pointer null\n");
		return NULL;
	}
	if(strlen(infile) < 1){
		fprintf(stderr, "Read File error, no file specified\n");
		return NULL;
	}

	/* Open file and check success */
	input = fopen(infile, "r");

	if(input == NULL){
		fprintf(stderr, "Read File error, file not opened\n");
		return NULL;
	}

	double tmpreal, tmpimag;
	char nwln;

	/* Read in real and imaginary number pairs from file*/
	for(i = 0; i < x->n; i++){
		fscanf(input, "(%f %fj) ", &tmpreal, &tmpimag);
		x->data[i] = tmpreal + tmpimag*I;

		/* Check if eof reached prematurely */
		if (feof(input)){
			fprintf(stderr, "File read error, eof reached early\n");
			return NULL;
	}
		
		/* Grab newline at end of each line */
		fscanf(input, "%c\n", &nwln); 
	}	

	fclose(input);

	return x;
}

