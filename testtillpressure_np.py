#!/usr/bin/env python2
"""
This script compares pressure calculated in tillPressureSoundNP with the results from tillPressure.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *

"""
Setup the plot.
"""
# Set a font
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 10.0

# Legend
# mpl.rcParams['legend.handlelength']  = 2.9
rcParams['legend.handlelength']  = 0.5
rcParams['legend.frameon']       = False
rcParams['legend.numpoints']     = 1
rcParams['legend.scatterpoints'] = 1

# Adjust axes line width
rcParams['axes.linewidth']   = 0.5

# Adjust ticks
rcParams['xtick.major.size'] = 4
rcParams['xtick.minor.size'] = 2
rcParams['ytick.major.size'] = 4
rcParams['ytick.minor.size'] = 2

# Adjust Font Size
rcParams['xtick.labelsize']  = 'x-small'
rcParams['ytick.labelsize']  = 'x-small'
rcParams['axes.labelsize']   = 'small'

# Set Up Figure, Single Column MNRAS
fig = gcf()
ax = gca()
fig, ax = subplots(1,1)
fig.set_size_inches(8.27*0.39,8.27*(6./8.)*0.39)

"""
data1 = loadtxt('out.old.dat',skiprows=3)
data2 = loadtxt('out.dat',skiprows=3)
"""

# Tillotson material
data1 = loadtxt('testtillpressure_np.txt')
data2 = loadtxt('testtillpressure.txt')
data3 = loadtxt('testtillpressure_np_diff.txt')

print where(fabs(data1-data2) > 1e-30)
print
print where(fabs(data3) > 1e-30);

figure(1)
imshow(data1-data2)

figure(2)
imshow(data1)

figure(3)
imshow(data2)

figure(4)
imshow(data3)

show()

