# Parallel FFT
A set of Fast Fourier Transform programs written in serial, using OpenMP, and MPI

Wesley Ellington

Daniel Reynolds

Math 4370

## Overview
This set of tools implements a simple FFT tool using various parallelization tools. It was designed mostly with MPI in mind,
so most performace benefits will be seen in that version on a large scale. It uses an iterative solution rather than the
more common recursive method, in an attempt to deal with larger problems that may lead to memory issues.

It is important to note that there is a specific format each program is expecting its input data to be in. Each folder 
contains a simple Python script to create input of this type and to serve as an example. In its current state, it must 
read in a text file (either all real, or complex will work) as input. Future plans include adding image filetypes.

Each different build will require special dependancies to construct. A Makefile exists in each version for easy compilation,
but may require a module loaded, or some minor tweaks.

All versions were tested and run on Southern Methodist University's ManeFrame supercomputer

## Serial
This version relies heavily on the use of a complex vector "class" for operation. All functions are written for this type of
object and must be supplied with the proper input format to run properly. In this stage, both 1D and 2D signals can be 
processed, but do assume problems sizes of 2^n nodes or 2^n by 2^n for 1D and 2D respectivly. 

This complex vector struct can handle most of its own I/O, meaning that it can be printed to console, read from a file, or 
dumped into a file at command.

## OpenMP
In order to build this version, it must be compiled on a machine that is OpenMP compatable. Otherwise, it will be no different
from the serial implementation, but should still compile.

Since very little can be done to parallelize the inner loops of a 1D FFT (at least it was of very little benefit), so most of
the performance benefit will be seen in the 2D version of the tool. Here, the whole signal was decomposed into strips and
processed individually by each thread.

## MPI
In terms of large scale computing, this is by far the most effective choice. It should be noted that this version requires a
power of two MPI tasks to run, but can scale to any square, power of two signal.

Unlike previous implementations, the vector structs have been removed in calculation stages
	This is to ensure that no index issues arise during data reordering from all to all

The overall proccess occurs as follows:

1. MPI is initialized, spawning N processes

2. All processes get their Id and number of procs in rank
	They will also know the problem size at this point

3. Process 0 loads all data for transform from file
	This will be done without a struct, dealing out data in read time as complex values
	0 will read in 1/N data as several rows, concatonate as a long array
	and push to respective procs


	Each process will receive 1/N the work for the comm group
	This includes proc 0, all processes will know total problem size

4. Processes will compute 1D transforms for each of their rows
	For example, if 4 procs exist in a 16 X 16 problem, each will do 4 rows

5. Each proc will "weave" data together to create correct distribution pattern for data
	EX: The first N values in the send buffer must be the first item in each resultant array,
	so the receving process must collect information and read it as the first value of the array
	then do so with each subsiquent value from the chunked data (each prox does more than one row)

6. Data is sent to new threads, decoded into proper indecies for iteration

7. FFTs are run over each chunk

8. All processes call back to 0, returning calculated values to base process

9. Information is stored in mother array

10. Data is output to file.

As of now, the MPI version cannot be called from outside programs, but must be run with input throught its driver.
Future versions will not have this limitation.
