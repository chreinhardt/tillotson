#!/usr/bin/env python2
"""
This script plots the results of plottillpressure.c.
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

"""
data1 = loadtxt('out.old.dat',skiprows=3)
data2 = loadtxt('out.dat',skiprows=3)
"""

"""
# Load values for rho and u
data = loadtxt('testtillpressure_np_diff.txt')

rho = data[:, 0]
u   = data[:, 1]
T   = data[:, 2]

print "rho_min=", min(rho), "rho_max=", max(rho)
print "u_min  =", min(u), "u_max  =", max(u)

# Convert to cgs
rho /= dGmPerCcUnit
u   /= dErgPerGmUnit

#data /= (dErgPerGmUnit*dGmPerCcUnit)
#data /= (dGmPerCcUnit*dErgPerGmUnit)

rho_min = min(rho)
rho_max = max(rho)

u_min   = min(u)
u_max   = max(u)

T_min   = min(T)
T_max   = max(T)

print "rho_min=", min(rho), "rho_max=", max(rho)
print "u_min  =", min(u), "u_max  =", max(u)
"""
# Reference density (material dependent)
rho0 = 2.7
us   = 3.5e10
us2  = 1.8e11

rho_min = 1e-4
rho_max = 10.0
u_min   = 0.0
u_max   = 25.0

rho_min = 1e-4
rho_max = 8.0602
u_min   = 0.0
u_max   = 19.8035

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
Plot the difference between tillPressureSoundNP and tillPressure.
"""
data = loadtxt('testtillpressure_np_diff.txt')

print where(fabs(data) < 1e-10)

imshow(data, origin='lower', extent=(rho_min, rho_max, u_min, u_max), aspect='auto')

plot([rho0, rho0], [u_min, u_max], '--', color='red')
plot([rho_min, rho_max], [us, us], '--', color='red')
plot([rho_min, rho_max], [us2, us2], '--', color='red')

xlabel("Density [code units]")
ylabel("Int. energy [code units]")

colorbar()

savefig('testtillpressure_np_diff.png', dpi=300, bbox_inches='tight')

fig.clear()

"""
Plot where the pressure is negative.
"""
data = loadtxt('testtillpressure_np_region.txt')

imshow(data, origin='lower', extent=(rho_min, rho_max, u_min, u_max), aspect='auto')

plot([rho0, rho0], [u_min, u_max], '--', color='red')
plot([rho_min, rho_max], [us, us], '--', color='red')
plot([rho_min, rho_max], [us2, us2], '--', color='red')

xlabel("Density [code units]")
ylabel("Int. energy [code units]")

colorbar()

savefig('testtillpressure_np_region.png', dpi=300, bbox_inches='tight')

fig.clear()

"""
Plot P(rho, T).
"""
data = loadtxt('testtillpressure_np_rhot.txt')

rho = data[:,0]

for i in range(1,size(data[1,:]),1):
    plot(data[:,0],data[:,i],'-',color='red',markersize=1,label='T')

rho_min = min(rho)
rho_max = max(rho)

xlim(rho_min, rho_max)
#ylim(0, 4)

xlabel("Density [code units]")
ylabel("Pressure [code units]")

savefig('testtillpressure_np_rhot.png', dpi=300, bbox_inches='tight')

show()


