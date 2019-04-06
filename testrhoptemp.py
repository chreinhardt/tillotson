#!/usr/bin/env python2
"""
This script plots the results of testrhoptemp.c.
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
rcParams['xtick.labelsize']  = 'xx-small'
rcParams['ytick.labelsize']  = 'xx-small'
rcParams['axes.labelsize']   = 'x-small'

# Set Up Figure, Single Column MNRAS
fig = gcf()
ax = gca()
fig, ax = subplots(1,1)
fig.set_size_inches(8.27*0.39,8.27*(6./8.)*0.39)

# Physical constants used for unit convertion (cgs)
KBOLTZ = 1.38e-16		# bolzman constant in cgs
MHYDR  = 1.67e-24		# mass of hydrogen atom in grams
MSOLG  = 1.99e33		# solar mass in grams
GCGS   = 6.67e-8		# G in cgs
KPCCM  = 3.085678e21		# kiloparsec in centimeters

# These two values define the unit system we use
dKpcUnit  = 2.06701e-13	        # kiloparsec in code units
dMsolUnit = 4.80438e-08	        # solar mass in code units

# Mass of Earth
MEarth = 5.98e27		# g
MSun   = 1.989e33	        # g

"""
Convert kboltz/mhydrogen to system units, assuming that
G == 1.
"""
# code energy per unit mass --> erg per g
dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM)
# code density --> g per cc
dGmPerCcUnit = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0)
# code time --> seconds
dSecUnit = sqrt(1/(dGmPerCcUnit*GCGS))

# The units are not nescessary but can be derived from the quanities above
Lunit = dKpcUnit*KPCCM
Munit = dMsolUnit*MSOLG

# Reference density (material dependent)
rho0 = 2.7
us   = 3.5e10
us2  = 1.8e11

"""
rho_min = 1e-4
rho_max = 10.0
u_min   = 0.0
u_max   = 25.0

rho_min = 1e-4
rho_max = 8.0602
u_min   = 0.0
u_max   = 19.8035
"""
P_min = 0.0
P_max = 1e3

T_min = 0.0
T_max = 1e3

"""
# Convert to cgs
rho_min *= dGmPerCcUnit
rho_max *= dGmPerCcUnit

u_min   *= dErgPerGmUnit
u_max   *= dErgPerGmUnit

print dGmPerCcUnit
print dErgPerGmUnit

print "rho_min=", rho_min, "rho_max=", rho_max
print "u_min  =", u_min, "u_max  =", u_max
"""

"""
Plot rho(P, T).
"""
data = loadtxt('testtillpressure_np_diff.txt')

print where(fabs(data) < 1e-10)

imshow(data, origin='lower', extent=(rho_min, rho_max, u_min, u_max), aspect='auto')

xlabel("Pressure [code units]")
ylabel("Temperature [K]")

colorbar()

savefig('testrhoptemp.png', dpi=300, bbox_inches='tight')

fig.clear()



