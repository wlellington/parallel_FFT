#!/usr/bin/env python

# Wesley Ellington
# Math 4370
# Serial FFT plotting functions

# Import plotting library
import matplotlib.pyplot as plt
import math
import re

############ Plot original Function #############

numbers = []

infile = open("samples.txt", "r")

#import entries from signal input file
for line in infile:
	numbers.append(float(line))

infile.close()

plt.figure()
plt.plot(numbers,"r-")
plt.axis([0, len(numbers), min(numbers), max(numbers)])
plt.title("Generated Function")
plt.xlabel("Sample")
plt.ylabel("Value")
plt.show()

############## Plot power Spectrum ##############

#master number list of data points entries
numbers = []

infile = open("powerSpec.txt", "r")

# Import entries from power spectrum output file
for line in infile:
	numbers.append(float(line) + 0.0000000000001)

infile.close()

plt.figure()
plt.plot(numbers,"b-")
plt.yscale("log")
plt.axis([0, 120, 0, (max(numbers))])
plt.title("Normalized Power Spectrum")
plt.xlabel("Frequency")
plt.ylabel("Relative power")
plt.show()

############# Plot reconstructed funciton ########

numbers = []

infile = open("rSignal.txt", "r")

# Import entries from fft output file
for line in infile:
	# Becuase the output files are now complex, we must use regex
	numbers.append(float(re.findall(r"-*\d+.\d*",line)[0]))

infile.close()

plt.figure()
plt.plot(numbers,"g-")
plt.yscale("linear")
plt.axis([0, len(numbers), min(numbers), max(numbers)])
plt.title("Reconstructed Signal")
plt.xlabel("Sample")
plt.ylabel("Value")
plt.show()
