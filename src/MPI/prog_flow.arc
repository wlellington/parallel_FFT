# Wesley Ellington
# Math 4370

# MPI parallelized FFT

This file contains a basic explaination of parallelized 
workflow for this project in MPI

Unlike previous implementations, the vector structs have been removed in calculation stages
	This is to ensure that no index issues arise during data reordering from all to all

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
