#!/usr/bin/env python2
"""
This script plots the results of plottillpressure.c.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *

"""
data1 = loadtxt('out.old.dat',skiprows=3)
data2 = loadtxt('out.dat',skiprows=3)
"""

# Tillotson material
data = loadtxt('plottillpressure.txt')

print where(fabs(data) < 1e-10)

figure(1)
imshow(data)

colorbar()

show()


