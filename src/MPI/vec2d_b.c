/* Wesley Ellington
 * Math 4370 
 * Homework 2
 * 2 Dimensional Vector Struct 
 * Using linear (* type) vector memory*/

#include "vec2d_b.h"
#include "time.h"
#include <math.h>
#include <stdio.h>

/* Implementation of 2D vector structs
 * Based on a double * linear type 
 * Index [i][j] is reached by [i*n + j]*/


/* Base constructor, allocates and zeros memory */
vec2d_b* vec2d_bNew(long int mlen, long int nlen){
	/* ensure lengths are legal*/
	int i;
	if (mlen < 1 || nlen < 1){
		fprintf(stderr, "vec2dNew error, illegal vector length");
		return NULL;
	}

	/* begin to allocate memory */

	/* allocate struct */
	vec2d_b *x = (vec2d_b*) malloc(sizeof(vec2d_b));
	
	/* if allocation fails, return null */
	if (x == NULL){
		return NULL;
	}
	
	/* assign struct vars */
	x->m = mlen;
	x->n = nlen;

	/* allocate data memory and set to zero*/
	x->data = (double*) calloc((mlen*nlen), sizeof(double*));

	return x;
}

/* Linear spaced constructor */
vec2d_b* vec2d_bLinspace(double a, double b, 
	long int m, long int n){
	int i, j;
	/* create new m by n vector */
	vec2d_b* x = vec2d_bNew(m,n);

	if(x == NULL){
		return NULL;
	}

	/* iterate over column, setting all values equal to
 	* linear spacing over m */ 	
	for(i = 0; i < m; i++){
		double val = a + (b-a)/(m-1)*i;
		for(j = 0; j < n; j++){
			x->data[(i*n) +j] = val;
		}
	}
	return x;
}

/* Randomly populated vector constructor */
vec2d_b* vec2d_bRandom(long int m, long int n){
	/*create vector*/
	vec2d_b* x = vec2d_bNew(m,n);

	int i,j;
	
	if(x->data == NULL){
		fprintf(stderr, "Rand Vec not initialized!");
	}
	double norm = pow(2.0,31.0) - 1.0;

	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			x->data[(i*n)+j] = random()/norm;
		}	
	}	
	return x;
}

/* Destructor, deallocates memory */
void vec2d_bDestroy(vec2d_b* inVec){
	int i;
	/* Dealloc outerlevel */
	free(inVec->data);
	/* Zero out values */
	inVec->m = 0;
	inVec->n = 0;
}

/* Create vector from input array and sizes */
vec2d_b* vec2d_bVectorize(double * inData, int size, int mlen, int nlen){
	vec2d_b * x;
	int i, j;

	if(mlen * nlen != size){
		fprintf(stderr, "Vectorize error, unmatched size");
		return NULL;
	}

	if(inData == NULL){
		fprintf(stderr, "Vectorize error, data array empty");
		return NULL;
	}

	/* Create new zero vector */
	x = vec2d_bNew(mlen, nlen);

	/* Assign values of array to vector
 	* This is done by casting array indicies in form [i*n + j] 
 	* one to one in vector object */
	for(i = 0; i < mlen; i ++){
		for(j = 0; j < nlen; j++){
			x->data[i*nlen + j] = inData[i*nlen + j];
		}
	}

	return x;

}


/* I/O functions */
/* Write to stdout*/
int vec2d_bWrite(vec2d_b* x){
	int i, j;
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			printf("%f ", x->data[(i*x->n)+j]);
		}
		printf("\n");
	}
	return 0;
}

/* Write to file */
int vec2d_bWriteFile(vec2d_b* x, const char* outfile){
	long int i, j;
	FILE * fptr = NULL;

	/* quite if vector has no alloced array*/
	if (x->data == NULL){
		fprintf(stderr, "Write to file error, no array\n");
		return 1;
	}

	/* return if no file string specified */
	if (strlen(outfile) < 1){
		fprintf(stderr, "Write for file error, no name\n");
		return 1;
	}

	/* open file and verify valid */
	fptr = fopen(outfile, "w");
	if(fptr == NULL){
		fprintf(stderr, "Write to file error, unable to open %s \n", outfile);
		return 1;
	}

	/* print to file */
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			fprintf(fptr, "%f ", x->data[(i*x->n)+j]);
		}
		fprintf(fptr, "\n");
	}

	/* close file */
	fclose(fptr);
	return 0;
} 

/* Read in from file */
/* subsets of data can be selected from file, but will be linear from top */
vec2d_b* vec2d_bReadFile(const char * infile, int m, int n){
	
	/* declare vars */
	FILE * input;
	int i, j;

	/* create new zeroed out vector of m by n */
	vec2d_b * x = vec2d_bNew(m, n);

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

	/* open file and check success */
	input = fopen(infile, "r");

	if(input == NULL){
		fprintf(stderr, "Read File error, file not opened\n");
		return NULL;
	}

	double tmp;
	char nwln;

	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			fscanf(input, "%lf", &tmp);
			x->data[i*n + j] = tmp;

			/* check if eof reached prematurely */
			if (feof(input)){
				fprintf(stderr, "File read error, eof reached early\n");
				return NULL;
			}
		}
		
		/* grab newline at end of each line */
		fscanf(input, "%c\n", &nwln); 
	}	

	fclose(input);

	return x;
}

/* Linear sum function */
/* x = a*y + b*z */
int vec2d_bLinearSum(vec2d_b* x, double a, vec2d_b* y,
		double b, vec2d_b* z ){
	int i, j;
	
	if(x->data == NULL || y->data == NULL || z->data ==NULL){
		fprintf(stderr, "vec2dLinsum error, null data\n");
		return 1;
	}

	if(x->m != y->m || x->m != z->m || x->n != y->n || x->n != z->n){
		fprintf(stderr, "vec2dLinsum error, diff sized arrays\n");
		return 1;
	}
	
	/* calculate values for each iteration*/
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			x->data[(i*x->n)+j] = y->data[(i*x->n)+j] * a + z->data[(i*x->n)+j] * b;
		}
	}
	return 0;
}

/* Scale function */
/* x = a*x */
int vec2d_bScale(vec2d_b* x, double a){
	int i,j;

	if(x->data == NULL){
		fprintf(stderr, "vec2dScale error, null data\n");
		return 1;
	}
	
	/*iterate over matix, scale each entry*/
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			x->data[(i*x->n)+j] = x->data[(i*x->n)+j] * a;
		}
	}
	return 0;
}

/* Vector copy */
/* x = y */
int vec2d_bCopy(vec2d_b* x, vec2d_b* y){
	int i, j;

	/*return with error for null data and diff sizes*/
	if(x->data == NULL || y-> data == NULL){
		fprintf(stderr, "vec2dCopy error, null data array");
		return 1;
	}
	
	if(x->m != y->m || x->n != y->n){
		fprintf(stderr, "vec2dCopy error, diffe sized arrays");
		return 1;
	}

	/*copy elements */

	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			x->data[(i*x->n)+j] = y->data[(i*x->n)+j];
		}
	}
	return 0;
}

/* Vector Constant */
/* sets all entries to a */
int vec2d_bConstant(vec2d_b* x, double a){
	int i,j;
	
	/*error if null data*/
	if(x->data == NULL){
		fprintf(stderr, "vec2dConstant error, null data array");
		return 1;
	}

	/* iterate over vector, set all entries to a */
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			x->data[(i*x->n)+j] = a;
		}
	}
	return 0;
}

/* Scalar values */

/*Vector Minimum*/
double vec2d_bMin(vec2d_b* x){
	int i, j;
	double min = x->data[0];

	/*Reset min if lower value is found*/
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			if(x->data[(i*x->n)+j] < min){
				min  = x->data[(i*x->n)+j];
			}
		}
	}
	return min;
}

/*Vector Maximum*/
double vec2d_bMax(vec2d_b* x){
	int i, j;
	double max = x->data[0];

	/*Reset max if higher value is found*/
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			if(x->data[(i*x->n)+j] > max){
				max  = x->data[(i*x->n)+j];
			}
		}
	}
	return max;
}

/*Vector dot product*/
double vec2d_bDot(vec2d_b* x, vec2d_b* y){
	int i, j;
	double sum = 0.0;

	if(x->m != y->m || x->n != y->n){
		fprintf(stderr, "vec2dDot error, different sized vecs\n");
		sum = 0.0;
		return sum;
	}

	/* iterate over entries, multiply pairs, sum*/	
	for(i = 0; i < x->m; i++){
		for(j = 0; j < x->n; j++){
			sum = sum + (x->data[(i*x->n)+j] * y->data[(i*x->n)+j]);
		}
	}	
	return sum;
}

/*Vector Two Norm*/
double vec2d_bTwoNorm(vec2d_b* x){
	double dot = vec2d_bDot(x,x);
	return sqrt((double)dot);
}

/*Root Mean Square Norm*/
double vec2d_bRmsNorm(vec2d_b* x){
	double dot = vec2d_bDot(x,x);
	return sqrt((double)dot/(double)(x->n * x->n));
}

/*Max norm/inf norm*/
double vec2d_bMaxNorm(vec2d_b* x){
	int i, j;	
	double max = 0.0;
	for(i = 0; i < x->n; i++){
		for(j = 0; j < x->n; j++){
			max = (max > fabs(x->data[(i*x->n)+j])) ? max : fabs(x->data[(i*x->n)+j]);
		}
	}
	return max;
}
