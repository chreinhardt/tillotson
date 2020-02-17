#!/usr/bin/env python2
"""
This script plots the results of testtillhugoniotpressure.c.
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
rcParams['legend.fontsize']      = 'xx-small'

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

# Physical constants used for unit convertion (cgs)
KBOLTZ = 1.38e-16		# bolzman constant in cgs
MHYDR  = 1.67e-24		# mass of hydrogen atom in grams
MSOLG  = 1.99e33		# solar mass in grams
GCGS   = 6.67e-8		# G in cgs
KPCCM  = 3.085678e21		# kiloparsec in centimeters

# These two values define the unit system we use
dKpcUnit  = 2.06701e-13 # kiloparsec in code units
dMsolUnit = 4.80438e-08 # solar mass in code units

# Mass of Earth
MEarth = 5.98e27		# g
MSun   = 1.989e33	    # g

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

data = loadtxt("testtillhugoniotpressure.txt")

rho = data[:, 0]
P_H = data[:, 1]
P_S = data[:, 2]

# Convert to cgs
rho *= dGmPerCcUnit

P_H *= (dErgPerGmUnit*dGmPerCcUnit)
P_S *= (dErgPerGmUnit*dGmPerCcUnit)

rho_min = min(rho)
rho_max = max(rho)

P_min   = min(min(P_H), min(P_S))
P_max   = max(max(P_H), max(P_S))

print "rho_min=", rho_min, "rho_max=", rho_max
print "P_min  =", P_min, "P_max  =", P_max

# Material parameters for granite
a     = 0.5;
b     = 1.3;
u0    = 1.6e11;
rho0  = 2.7;
A     = 1.8e11;
B     = 1.8e11;
us    = 3.5e10;
us2   = 1.8e11;
alpha = 5.0;
beta  = 5.0;

c1 = '#66c2a5'
c2 = '#fc8d62'
c3 = '#8da0cb'

# Calculate the linear Hugoniot
Gamma0 = a + b
s = 0.5*(1 + B/A + 0.5*(a+b))

rho_max = rho0*s/(s-1.0)

print "rho_max=", rho_max

print rho_max/rho0
rho_max = 2*rho0

# Limit the density
i = where(rho < rho_max)

rho = rho[i]
P_H = P_H[i]
P_S = P_S[i]

rho_min = min(rho)
rho_max = max(rho)

P_min   = min(min(P_H), min(P_S))
P_max   = max(max(P_H), max(P_S))

# Add the reference point to the plot
P0 = 0.0

rho = numpy.append(rho0, rho)
P_H = numpy.append(P0, P_H)
P_S = numpy.append(P0, P_S)

# Mark the reference state
#scatter(rho0, P0, s=10, color=c3)

"""
Plot the Hugoniot.
"""
plot(rho, P_H, '-', color=c1)
plot(rho, P_S, '--', color=c2)

xlim(rho0, rho_max)
ylim(-1e10, P_max)

i = where(abs(rho-1.8*rho0)<5e-2)

rho1 = min(rho[i])
P1   = min(P_H[i])

# Annotate the Hugoniot and the cold cruve
rho_min = 0.0
P_min = 0.0

"""
i_min = where(abs(rho-4.1)<5e-2)
print i_min

i_max = where(abs(rho-4.1)<5e-2)
print i_max

x_text_min = min(rho[i_min])

x_text_max = 4.7

y_text_min = 0.35e12
y_text_max = 4.7
"""

annotate("Hugoniot", xy=(rho1, P1), xytext=(0.65, 0.45), textcoords='axes fraction', horizontalalignment='center', fontsize='x-small', rotation=41) 
#annotate("Hugoniot", xy=(0,0), xytext=(rho_min+0.5*(rho_max-rho_min), P_min+0.5*(P_max-P_min))) 

annotate("Isentrope", xy=(rho1, P1), xytext=(0.65, 0.245), textcoords='axes fraction', horizontalalignment='center', fontsize='x-small', rotation=25.5) 
print rho1, P1

#scatter(rho1, P1, s=25, color='blue')
#plot([rho0, rho1], [P0, P1], '-', color='red')

xlabel("Density [g]")
ylabel("Pressure [erg cm$^3$]")

savefig('testtillhugoniotpressure.png', dpi=300, bbox_inches='tight')
savefig('testtillhugoniotpressure.pdf', dpi=300, bbox_inches='tight')
#show()

