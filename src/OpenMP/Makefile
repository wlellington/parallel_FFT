# Makefile for FFT serial driver implementation
#
# Wesley Ellington
# Math 4370
# March 2, 2017
#
# This file builds a serial implementaion of an FFT algorithm

# Compiler flags
CC = gcc
CFLAGS = -O0 -fopenmp 

# Build targets

all : driver.exe


driver.exe : driver.c vec2d_b.o cvec1d.o fft.o get_time.o
	$(CC) $(CFLAGS) $^ -lm -o $@

vec2d_b.o : vec2d_b.c vec2d_b.h
	$(CC) $(CFLAGS) -lm -c $< -o $@

cvec1d.o : cvec1d.c cvec1d.h
	$(CC) $(CFLAGS) -lm -c $< -o $@

fft.o : fft.c fft.h vec2d_b.h cvec1d.h
	$(CC) $(CFLAGS) -lm -c $< -o $@

get_time.o : get_time.c get_time.h
	$(CC) $(CFLAGS) -lm -c $< -o $@

clean :
	\rm -f *.o *.out

realclean : clean 
	\rm -f *.exe *~

# End Makefile
