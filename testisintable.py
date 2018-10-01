#!/usr/bin/env python2
"""
This script plots the data from testisintable.c.
"""
from matplotlib import *
from matplotlib.pyplot import *
from numpy import *
from sys import *

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

# Load the whole lookup table
data = loadtxt('lookup.txt')

# Load data points
data2 = loadtxt('testisintable.txt')
data3 = loadtxt('testisintable2.txt')

rho = data2[:,0]
u = data2[:,1]

# Number of isentropes
v_n = size(data[0, :])

# Plot the lookup table
for i in range(1, v_n, 1):
		plot(data[:,0], data[:,i], '-', color='red', markersize=1, label='Table')

scatter(rho, u, s=16, color='blue')


xmax = 25
ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

savefig('testisintable.png', dpi=300, bbox_inches='tight')

figg = gcf()

fig.clear()

rho = data3[:,0]
u = data3[:,1]
error = data3[:,2]

# Plot the lookup table
for i in range(1, v_n, 1):
		plot(data[:,0], data[:,i], '-', color='red', markersize=1, label='Table')

scatter(rho, u, s=16, c=error, linewidth=0)

colorbar()

xmax = 25
ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

savefig('testisintable2.png', dpi=300, bbox_inches='tight')

show()

