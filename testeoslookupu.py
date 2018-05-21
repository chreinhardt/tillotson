#!/usr/bin/env python2
"""
A script to debug eosLookupU().
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
Plot numerical solution.
"""
data = loadtxt('testeoslookupu.txt')

nU = 10

# Store the different isentropes in rho[0], ... , rho[nU-1]
rho_num = numpy.split(data[:, 0], nU)
u_num   = numpy.split(data[:, 1], nU)

for i in range(nU):
    plot(rho_num[i], u_num[i])


"""
Plot analytic solution.
"""
data = loadtxt('testeoslookupu2.txt')

# Store the different isentropes in rho[0], ... , rho[nU-1]
rho_ana = numpy.split(data[:, 0], nU)
u_ana   = numpy.split(data[:, 1], nU)

for i in range(nU):
    plot(rho_ana[i], u_ana[i], '--')

xlabel(r'Density [g/cm^3]')
ylabel(r'Internal energy [erg/g]')

#ylim(0, min([max(u_num), max(u_ana)]))

#xlim(min(Vr), max(Vr))
#xlim(0.4, 3.0)
#ylim(0, 1.6)

title('mu = 23.0')
legend(loc='best')

savefig('testeoslookupu.png', dpi=300, bbox_inches='tight')

show()

