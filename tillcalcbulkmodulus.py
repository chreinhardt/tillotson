#!/usr/bin/env python2
"""
Plot the bulk modulus.
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

# Load the data
data1 = loadtxt('bulkmodulus_granite.txt')
data2 = loadtxt('bulkmodulus_ice.txt')
data1 = loadtxt('bulkmodulus_ice.txt')
data2 = loadtxt('bulkmodulus_ice_alpha_beta.txt')
# Plot the bulk modulus
for i in range(1,size(data1[0,:]),1):
		plot(data1[:,0],data1[:,i],'-',color='red',markersize=1)
#		semilogx(data[:,0],data[:,i],'-',color='red',markersize=1,label='Table')

for i in range(1,size(data2[0,:]),1):
		plot(data2[:,0],data2[:,i],'--',color='blue',markersize=1)
#		semilogx(data[:,0],data[:,i],'-',color='red',markersize=1,label='Table')

xmax = 25
ymax = 25
#xlim(0,xmax)
ylim(0,1e5)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Bulk modulus')

savefig('tillcalcbulkmodulus.png', dpi=300, bbox_inches='tight')
show()
