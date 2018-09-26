#!/usr/bin/env python2
"""
This script produces a plot of the cold curve for the Tillotson EOS.
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

# Load the cold curve
data1 = loadtxt('coldcurve.txt')

# Load interpolated values
data2 = loadtxt('testtillcoldu.txt')

rho_cold = data1[:,0]
u_cold   = data1[:,1]

rho      = data2[:,0]
u        = data2[:,1]

# Plot the lookup table
plot(rho_cold, u_cold, '-', color='red', markersize=1, label='Table')

# Plot the interpolated data
plot(rho, u, '--', color='blue', markersize=1, label='Table')

xmax = 25
ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

savefig('testtillcoldu.png', dpi=300, bbox_inches='tight')
show()
