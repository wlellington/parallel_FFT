#!/usr/bin/env python

# Wesley Ellington
# Math 4370
# FFT 2d Signal generator

import cmath
import math
import numpy

pi = cmath.pi

# Number of nodes per direction, n by n
n = 4096

#Coefficients in x
xcoeffs = []

#Coefficients in y
ycoeffs = []

print ("Generating mesh of " + str(n ** 2) + " nodes")

output = open("2Dlarge.txt", "w")

for y in range(0, n):
	for x in range(0, n):
		val = 0;
		t = float(x) / float(n)
		s = float(y) / float(n)
		val = math.sin(2*pi*t)
		output.write('{0:.16f}'.format(val) + " ")

	output.write("\n")

output.close()

output = open("2Dcount.txt", "w")

for i in range(0, n ** 2):
	if i % n == 0:
		output.write("\n")
	output.write(str(i) + " ")

output.close()

