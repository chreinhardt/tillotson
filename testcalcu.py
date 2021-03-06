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

# Load the whole lookup table
data = loadtxt('lookup.txt')

# Load interpolated values
data2 = loadtxt('testsplint.txt')

# Load particles
data3 = loadtxt('testcalcu.txt')

rho1 = data3[:,0]
u1   = data3[:,1]

rho2 = data3[:,2]
u2   = data3[:,3]

# Plot the lookup table
for i in range(1,size(data[0,:]),1):
		plot(data[:,0],data[:,i],'-',color='red',markersize=1,label='Table')
#		semilogx(data[:,0],data[:,i],'-',color='red',markersize=1,label='Table')

# Plot the interpolated data
for i in range(1,size(data2[0,:]),1):
		plot(data2[:,0],data2[:,i],'-',color='green',markersize=1,label='Lookup')
#		semilogx(data2[:,0],data2[:,i],'-',color='green',markersize=1,label='Lookup')

scatter(rho1, u1, color='blue', s=5)
scatter(rho2, u2, color='blue', s=5)

xmax = 25
ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

#plot(rho,u,'.',color='blue',markersize=1,label='Lookup')

#plot(rho,u,'-',color='blue',linewidth=1,label='Lookup table')

#fill_between(rhocold,ucold,color='orange')

#plot([0,rho0],[us,us],'r--')
#plot([0,rho0],[us2,us2],'r--')
#plot([rho0,rho0],[0,25],'r--')

#xticks([0,rho0,2*rho0,3*rho0],[r'$0$',r'$\rho_0$',r'$2\rho_0$',r'$3\rho_0$'])
#xticks([0,rho0],[r'$0$',r'$\rho_0$'],size='large')
#yticks([0,us,us2],[r'$0$',r'$u_{IV}$',r'$u_{CV}$'],size='large')

savefig('testsplint.png', dpi=300, bbox_inches='tight')
show()
