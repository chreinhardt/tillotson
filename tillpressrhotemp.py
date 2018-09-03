#!/usr/bin/env python2
"""
Plot a ballic model.
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

"""
if len(argv) != 3:
    print "Usage: profile-plot.py <profile> <model>"
    exit(1)

profile = argv[1]
"""

"""
Units in cgs assuming Lunit = 1RE and vunit = 1km/s.
"""
ErgPerGmUnit = 9.99823e+09  # erg/g
GmPerCcUnit = 0.368477      # g/cm^3
SecUnit = 6378.69           # s
Lunit = 6.37813e+08
Munit = 9.56072e+25

data = loadtxt('press_rho_temp.txt')

rho = data[:,0]

for i in range(1,size(data[:,0]),1):
		plot(data[:,0],data[:,i],'-',color='red',markersize=1,label='T')

"""
xlim(0, max(max(R), max(Rmodel)*1.05))
ylim(0, max(max(rho), max(rhomodel)*1.05))
"""

xlabel(r'Density [code units]')
ylabel(r'Pressure [code units]')
#ylabel(r'Density [g/cm$^3$]')
savefig('press_rho_temp.png', dpi=300, bbox_inches='tight')

