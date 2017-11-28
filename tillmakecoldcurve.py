#!/usr/bin/env python2
"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *

# Load the lookup table for the cold curve
data1 = loadtxt('lookup.txt')

# Load data from direct integration
data2 = loadtxt('out')

"""
Load data.
"""
rho1 = data1[:,0]
u1   = data1[:,1]

rho2 = data2[:,0]
u2   = data2[:,1]

"""
xmax = 25
ymax = 25
xlim(0, xmax)
ylim(0, ymax)
"""

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

plot(rho1, u1, '-', color='red', label='Cold curve')
plot(rho2, u2, '--', color='blue', label='Cold curve')

savefig('tillmakecoldcurve.png')
show()

