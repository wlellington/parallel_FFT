# Makefile for FFT serial driver implementation
#
# Wesley Ellington
# Math 4370
# March 2, 2017
#
# This file builds a serial implementaion of an FFT algorithm

# Compiler flags
MPICC = mpicc
CFLAGS = -O0 

# Build targets

all : driver_mpi.exe


driver_mpi.exe : driver_mpi.c vec2d_b.o cvec1d.o fft_mpi.o
	$(MPICC) $(CFLAGS) $^ -lm -g -o $@

vec2d_b.o : vec2d_b.c vec2d_b.h
	$(MPICC) $(CFLAGS) -lm -c $< -g -o $@

cvec1d.o : cvec1d.c cvec1d.h
	$(MPICC) $(CFLAGS) -lm -c $< -g -o $@

fft_mpi.o : fft_mpi.c fft_mpi.h
	$(MPICC) $(CFLAGS) -lm -c $< -g -o $@

clean :
	\rm -f *.o *.out

realclean : clean 
	\rm -f *.exe *~

# End Makefile
