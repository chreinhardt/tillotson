#!/usr/bin/env python2
"""
This script compares the pressure obtained from tilleos.c with the one from tillotson.c.
obtained from gendatotipsy.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *

data1 = loadtxt('soundspeed_cutoff1.txt')
data2 = loadtxt('soundspeed_cutoff2.txt')

# Values for granite (in code units!)
us = 3.5
us2 = 18.0
rho0 = 7.33

print where(fabs(data1-data2) > 1e-30)
print

figure(1)
imshow(data1-data2)
#xlim(0,20)
#ylim(0,7.33)
xlabel('Int. energy')
ylabel('Density')
colorbar()

figure(2)
imshow(data1)
xlabel('Int. energy')
ylabel('Density')
colorbar()

figure(3)
imshow(data2)
xlabel('Int. energy')
ylabel('Density')
colorbar()
show()

