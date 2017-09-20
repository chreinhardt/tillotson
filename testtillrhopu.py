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

# Tillotson material
data1 = loadtxt('testtillrhopu1.txt')
data2 = loadtxt('testtillrhopu2.txt')

print "data1-data2:"
print where(fabs(data1-data2) > 1e-30)
print

"""
index = where(fabs(data1-data2) > 1e-30)
print
print "data1-data2:"
print data1[index]
print
print data2[index]
"""
figure(1)
imshow(data1-data2)
xlim(0,20)
ylim(0,7.33)
xlabel('Int. energy')
ylabel('Density')
colorbar()

figure(2)
imshow(data1-data2)
xlabel('Int. energy')
ylabel('Density')
colorbar()

figure(3)
imshow(data1)
xlabel('Int. energy')
ylabel('Density')
colorbar()

figure(4)
imshow(data2)
xlabel('Int. energy')
ylabel('Density')
colorbar()

data3 = loadtxt('testtillrhopu.txt')

rho = data3[:,0]
u = data3[:,1]
P = data3[:,2]

#
# Plot the points where the pressure differ in the u(rho) plane.
#
figure(5)

"""
The region where the pressure is zero.
"""
np = loadtxt('till_granite_np.txt')
rho_np = np[:,0]
u_np = np[:,1]
print size(np)

# The cold curve for granite.
coldcurve = loadtxt('/home/ics/creinh/code/ballic-test/tools/cold_granite.txt')

# Values for granite (in code units!)
us = 3.5
us2 = 18.0
rho0 = 7.33

rhocold = coldcurve[:,0]
ucold = coldcurve[:,1]


scatter(rho,u,color='blue')
#plot(rhocold,ucold,color='red')
xmax = max(rho)*1.1
ymax = max(u)*1.1

#plot(rho_np,u_np,color='red')

xlim(0,xmax)
ylim(0,ymax)

plot([0,rho0],[us,us],'r--')
plot([0,rho0],[us2,us2],'r--')
plot([rho0,rho0],[0,ymax],'r--')

xlabel('Density')
ylabel('Int. energy')

show()

