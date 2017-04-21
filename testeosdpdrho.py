"""
This script compares the pressure obtained from tilleos.c with the one from tillotson.c.
obtained from gendatotipsy.
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
data1 = loadtxt('testeosdpdrho1.txt')
data2 = loadtxt('testeosdpdrho2.txt')

# Ideal gas
data3 = loadtxt('testeosdpdrho3.txt')

print where(fabs(data1-data2) > 1e-30)
print
print where(fabs(data3) > 1e-30);

figure(1)
imshow(data1-data2)

figure(2)
imshow(data1)

figure(3)
imshow(data2)

figure(4)
imshow(data3)

show()

