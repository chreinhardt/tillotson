#!/usr/bin/env python2
"""
This script compares the density obtained from tillRhoPU() with the expected value.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *

"""
data1 = loadtxt('out.old.dat',skiprows=3)
data2 = loadtxt('out.dat',skiprows=3)
"""

# Load data
P1 = loadtxt('testtillidealgas.pressure1.txt')
P2 = loadtxt('testtillidealgas.pressure2.txt')

cs1 = loadtxt('testtillidealgas.soundspeed1.txt')
cs2 = loadtxt('testtillidealgas.soundspeed2.txt')

dPdrho1 = loadtxt('testtillidealgas.dPdrho1.txt')
dPdrho2 = loadtxt('testtillidealgas.dPdrho2.txt')

dPdu1 = loadtxt('testtillidealgas.dPdu1.txt')
dPdu2 = loadtxt('testtillidealgas.dPdu2.txt')

"""
P(rho, u)
"""
print "P1-P2:"
print where(fabs(P1-P2) > 1e-30)
print

figure(1)
imshow(P1-P2)
#xlim(0,20)
#ylim(0,7.33)
xlabel('Int. energy')
ylabel('Density')
title('P1-P2')
colorbar()

"""
cs2(rho, u)
"""
print "cs1-cs2:"
print where(fabs(cs1-cs2) > 1e-30)
print

figure(2)
imshow(cs1-cs2)
xlabel('Int. energy')
ylabel('Density')
title('cs1-cs2')
colorbar()

"""
dPdrho(rho, u)
"""
print "dPdrho1-dPdrho2:"
print where(fabs(dPdrho1-dPdrho2) > 1e-30)
print

figure(3)
imshow(dPdrho1-dPdrho2)
xlabel('Int. energy')
ylabel('Density')
title('dPdrho1-dPdrho2')
colorbar()

"""
dPdu(rho, u)
"""
print "dPdu1-dPdu2:"
print where(fabs(dPdu1-dPdu2) > 1e-30)
print

figure(4)
imshow(dPdu1-dPdu2)
xlabel('Int. energy')
ylabel('Density')
title('dPdu1-dPdu2')
colorbar()

show()

