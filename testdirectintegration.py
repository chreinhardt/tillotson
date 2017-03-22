"""
This script plot the lookup table and one isentrope, that was directly integrated using tillCalcU().
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# This might be useful later to debug our look up table.
lookup = loadtxt('lookup.txt')
data = loadtxt('testdirectintegration.txt')

# Cold curve
rhocold = lookup[:,0]
ucold = lookup[:,1]

rho = data[:,0]
u = data[:,1]
u2 = data[:,2]

#rho2 = data[:,2]
#u = data[:,1]

# Values for granite
us = 3.5
us2 = 18.0
rho0 = 7.33

"""
xmax = 25
ymax = 25
"""

xmax = min(max(rho)*1.1,max(rhocold))
ymax = max(u)*1.1

xlim(0,xmax)
ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

print "rho=",[data[0],data[2]],"u=",[data[1],data[3]]
#plot(rho,u,'.',color='green',markersize=0.1,label='Lookup table')
plot(rhocold,ucold,'-',color='red',linewidth=2,label='Cold curve (T=0)')

"""
plot([data[0],data[2]],[data[1],data[3]],'-',color='blue')
scatter(data[0],data[1],s=10,color='blue')
"""
plot(rho,u,'-',color='blue',linewidth=2,label='Direct integration')
plot(rho,u2,'--',color='red',linewidth=2,label='Lookup')
# Plot different isentropes (Lookup(i, j) = Lookup(rho,v))
for i in arange(1,1000,10):
		plot(rhocold,lookup[:,i],'-',color='green',markersize=1)

#plot(rho,u,'-',color='blue',linewidth=1,label='Lookup table')

#fill_between(rhocold,ucold,color='orange')

plot([0,rho0],[us,us],'r--')
plot([0,rho0],[us2,us2],'r--')
plot([rho0,rho0],[0,ymax],'r--')

#xticks([0,rho0,2*rho0,3*rho0],[r'$0$',r'$\rho_0$',r'$2\rho_0$',r'$3\rho_0$'])
#xticks([0,rho0],[r'$0$',r'$\rho_0$'],size='large')
#yticks([0,us,us2],[r'$0$',r'$u_{IV}$',r'$u_{CV}$'],size='large')

show()

savefig('testdirectintegration.png')

