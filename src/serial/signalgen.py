#!/usr/bin/env python

# Wesley Ellington
# Math 4370
# FFT signal Generator

import cmath
import math
import numpy

pi = cmath.pi
# List of frequencies to generate data with
#coeffs = [1,7,20]

# Secant?
coeffs = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

# Square Wave
#coeffs = [1.0,3.0,5.0,7.0,9.0, 11.0, 13.0 ]

# List of scalars for freqs, must be one to one
#mags = [1,.75,.5]

# Secant?
mags = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# Square Wave
#mags = [4/(1.0*pi), 4/(3.0*pi), 4/(5.0*pi), 4/(7.0*pi), 4/(9.0*pi), 4/(11.0*pi), 4/(13.0*pi)]

# List of data points created from Fourier series
samples = []

# Seconds of signal
n = 1

# Samples per second
Fs = 8192

print ("Generating signal of size: " + str(n*Fs))

########### Signal Creation Phase ##############
for x in range(0,n):
	for step in range(0,Fs):
		val = float(0.0)
		t = float(x) + float(step)/float(Fs)
		for i in range(0, len(coeffs)):
			# calculate contribution from one frequency element and magnitude
			val += float(mags[i])*math.sin(2.0*math.pi*float(coeffs[i])*float(t))
		# Set sample value for time step
		samples.append(float(val))

#============ Noise Addition =============#
noise = numpy.random.normal(0,.1,len(samples))

# Comment out next line for pure signal
#samples = samples + noise

#output values to file
outfile = open("samples.txt", "w")

#write value to file
for value in samples:
	outfile.write(str(str(value) + "\n"))

outfile.close()


# Create perfect square wave
outfile = open("rectangle.txt", "w")

for value in range(0,2048):
	outfile.write(str(str(1) + "\n"))
for value in range(0,2048):
	outfile.write(str(str(-1) + "\n"))
for value in range(0,2048):
	outfile.write(str(str(1) + "\n"))
for value in range(0,2048):
	outfile.write(str(str(-1) + "\n"))

outfile.close()
