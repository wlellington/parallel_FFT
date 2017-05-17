#!/bin/bash

echo "Preparing Build"

# if not loaded, load mpavich
if [ "$MPIHOME" == "/grid/software/mvapich2/2.0/gcc-4.9.1" ]
then
	echo "MPI already loaded!"
else
	echo "Loading MPI"

	eval "module load mvapich2/2.0-gcc-4.9.1"
	export $MPIHOME
	export 
fi

eval "env | grep \"MPI\""

#clean from previous build
make clean

#build new target
make

echo "Build Complete"
