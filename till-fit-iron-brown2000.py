#!/usr/bin/env python2
"""
Calculate the Gruneisen parameter of the Tillotson EOS.
"""
from matplotlib import *
from matplotlib.pyplot import *
import numpy
from scipy.integrate import odeint
from scipy import optimize
import sys


def Gamma(rho, u, rho0, u0, a, b):
    """
    Calculate the generalized Gruneisen parameter of the Tillotson EOS.

    Paramters:
    ----------
    rho:    Density
    u:      Internal energy

    rho0, u0, a, b are Tillotson EOS material parameters     

    Return:
    -------
    Gamma:  Generalized Gruneisen parameter
    """
    eta = numpy.array(rho/rho0) 
    eta2 = eta*eta

    omega0 = u/(u0*eta2)+1.0

    Gamma = a + b/omega0

    return numpy.array(Gamma)


def pressure(rho, u, rho0, u0, a, b, A, B):
    """
    Calculate the pressure in the condensed states of the Tillotson EOS.

    Paramters:
    ----------
    rho:    Density
    u:      Internal energy

    rho0, u0, a, b, A, B are Tillotson EOS material parameters     

    Return:
    -------
    P:      Pressure
    """
    eta = rho/rho0
    mu = eta-1.0

    P = Gamma(rho, u, rho0, u0, a, b)*rho*u + A*mu + B*mu**2

    return P


def pressure_hugoniot_lin(rho, C0, s, rho0, u0, a, b, A, B):
    """
    Calculate the pressure along the principal hugoniot of the Tillotson EOS assuming a linear
    relation
    
    Us = C0 + s*Up

    between Us and Up.

    Paramters:
    ----------
    rho:    Density
    C0:     Intercept
    s:      Slope

    rho0, u0, a, b, A, B are Tillotson EOS material parameters     

    Return:
    -------
    P:      Pressure
    """

    """
    # Us = C0 + s*Up + q*Up^2
    C0 = 3.691
    s  = 1.788
    q  = -0.038
    """
    
    # Calculate u(rho, P) from the Rankine-Hugoniot equations 
    X = 1.0 - rho0/rho
    P_H = rho0*X*C0**2/(1.0 - s*X)**2

    u = 0.5*P_H*(1.0/rho0 - 1.0/rho)

    P = pressure(rho, u, rho0, u0, a, b, A, B)

    return P


def pressure_hugoniot_wrapper_lin(C0, s, rho0, a, b, A):
    """
    A wrapper of pressure_hugoniot_lin for optimize.curve_fit.

    Paramters:
    ----------
    C0:     Intercept
    s:      Slope
    rho0, a, b, A, re Tillotson EOS material parameters     

    Return:
    -------
    func:   A function that returns P(rho) for optimize.cuve_fit
    """
    def func(rho, u0, B):
        
        print "u0=", u0, "B=", B
        P = pressure_hugoniot_lin(rho, C0, s, rho0, u0, a, b, A, B)
        return P

    return func


def pressure_hugoniot(rho, C0, s, q, rho0, u0, a, b, A, B):
    """
    Calculate the pressure along the principal hugoniot of the Tillotson EOS assuming
    
    Us = C0 + s*Up + q*Up^2

    between Us and Up.

    Paramters:
    ----------
    rho:        Density
    C0, s, q:   Fitting parameters

    rho0, u0, a, b, A, B are Tillotson EOS material parameters     

    Return:
    -------
    P:          Pressure
    """

    """
    # Us = C0 + s*Up + q*Up^2
    C0 = 3.691
    s  = 1.788
    q  = -0.038
    """
    
    Us = C0 + s*Up + q*Up**2

    # Calculate u(rho, P) from the Rankine-Hugoniot equations 
    P_H = rho0*Us*Up

    u = 0.5*P_H*(1.0/rho0 - 1.0/rho)

    P = pressure(rho, u, rho0, u0, a, b, A, B)

    return P


def pressure_hugoniot_wrapper(C0, s, rho0, a, b, A):
    """
    A wrapper of pressure_hugoniot for optimize.curve_fit.

    Paramters:
    ----------
    C0, s, q:   Fitting parameters

    rho0, a, b, A, re Tillotson EOS material parameters     

    Return:
    -------
    func:   A function that returns P(rho) for optimize.cuve_fit
    """
    def func(rho, u0, B):

        P = pressure_hugoniot(rho, C0, s, q, rho0, u0, a, b, A, B)
        return P

    return func



def dudrho(u, rho, param):
    """
    Calculate

    du/drho = P/rho^2

    from the Tillotson EOS for odeint.

    Paramters:
    ----------
    rho:    Density
    u:      Internal energy
    param:  Material parameters

    Return:
    -------
    dudrho: Derivative the internal energy at constant entropy

    Note: scipy.integrate.odeint requires the function to be of the form func(y, t)
    """
    rho0 = param[0]
    u0   = param[1]
    a    = param[2]
    b    = param[3]
    A    = param[4]
    B    = param[5]

    return pressure(rho, u, rho0, u0, a, b, A, B)/rho**2


def main():
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
    Material parameters of iron from Benz et al. (1986).
    """
    a     = 0.5
    b     = 1.5
    u0    = 9.5e10
    rho0  = 7.86
    A     = 1.279e12
    B     = 1.05e12
    us    = 1.42e10
    us2   = 8.45e10
    alpha = 5.0
    beta  = 5.0

    """
    Material parameters of iron from Wissing et al. (2020).
    """
    a_new    = 0.5
    b_new    = 1.28
    u0_new   = 1.425e11
    rho0_new = 7.85
    A_new    = 1.28e12
    B_new    = 1.815e12

    num_rho = 1000
    num_u   = 100

    """
    Set material constants for the Tillotson EOS. The fitting parameters u0 and B are obtained by
    fitting to Hugoniot data from Brown et al. (2000).
    """
    a     = 0.5         # a = 0.5 provides the best asympthotic fit (Tillotson 1962)
    b     = 1.5         # Gamma0 = a + b, again from Tillotson (1962)
    u0    = 0.0
    rho0  = 7.85        # Reference density at zero compression from Brown et al (2000)
    A     = 1.28e12     # Bulk modulus at the reference state (Tillotson 1962)
    B     = 0.0
    us    = 0.0
    us2   = 0.0
    alpha = 5.0         # These values are used for most solids, why they were used originally
    beta  = 5.0         # by Tillotson is unclear from the report

    """
    Hugoniot data for iron from Brown et al. (2000).
    """
    data = numpy.loadtxt("hugoniot_iron_brown2000.txt")

    Up    = data[:,0]  # km/s
    Us    = data[:,1]  # km/s
    P_H   = data[:,2]  # GPa
    rho_H = data[:,3]  # g/cm^3

    # Convert to cgs
    Up  *= 1e5
    Us  *= 1e5
    P_H *= 1e10         # 1 GPa = 1e10 erg/cm^3

    """
    Fitting parameters from Brown et al. (2000).
    """
    #C0 = 3.691
    #s  = 1.788
    #q  = -0.038

    # Linear Hugoniot from Brown et al. (2000)
    C0 = 3.935
    s  = 1.578
    
    # Convert to cgs (note that s is not converted)
    C0 *= 1e5
    
    print "Parameters:"
    print "rho0=", rho0
    print "C0=  ", C0
    print "s=   ", s

    rho_max = numpy.max(rho_H)
    rho = numpy.linspace(rho0, rho_max, 1000)
    
    #X = 1.0 - rho0/rho
    #P = rho0*X*C0**2/(1.0 - s*X)**2

    """
    X = 1.0 - rho0/rho
    a1 = q*X**2
    a2 = (s*X-1.0)
    a3 = C0

    plot(X, a1*X**2+a2*X+a3, color='red')
    show()
    exit(1)
    Us = C0 + s*Up + q*Up**2
    """

    """
    # Plot the data
    scatter(rho_H, P_H, color='red', marker='x', label="Hugoniot")
    
    P = pressure_hugoniot_lin(rho, C0, s, rho0, u0, a, b, A, B)
    plot(rho, P, color='blue', label="Pressure")
    
    show()
    exit(1)
    """

    popt, pcov = optimize.curve_fit(pressure_hugoniot_wrapper_lin(C0, s, rho0, a, b, A), rho_H, P_H)
    #popt, pcov = optimize.curve_fit(pressure_hugoniot_wrapper(rho0, a, b, A), numpy.array([rho_H, Us, Up]), P_H, bounds=(0, numpy.inf))

    u0 = popt[0]
    B  = popt[1]

    print popt
    print pcov
    exit(1)


    rho_max = 1e5*rho0
    u_max   = 10.0*us2
    u_initial = numpy.linspace(0.0, u_max, num_u)
    #u_initial = numpy.logspace(0.0, numpy.log(u_max), num_u, base=numpy.e)

    print "Solving isentropes:"
    print "u_initial=", u_initial
    print

    """
    Calculate an isentrope.
    """
    rho = numpy.linspace(rho0, rho_max, num_rho)
    u = odeint(dudrho, u_initial, rho, args=([rho0, u0, a, b, A, B], ))

    # odeint returns an array for each initial value
    for i in range(0, num_u):
        G = Gamma(rho, u[:,i])
        #print G
        #plot(rho/rho0, G, linewidth=1)
        semilogx(rho/rho0, G, linewidth=1)
        #plot(rho/rho0, u[:, i], color='red', linewidth=1)

    # Mark the limits
    plot([1.0, rho_max/rho0], [0.5, 0.5], '--', color='red')
    plot([1.0, rho_max/rho0], [1.8, 1.8], '--', color='red')

    ylim(0.4, 1.9)

    xlabel(r'Normalized Density ($\rho$/$\rho_0$)')
    ylabel(r'Generalized Gruneisen Parameter')

    legend(loc='best')

    savefig('tillotson-gamma-isentrope.png', dpi=300, bbox_inches='tight')
    
    fig = gcf()
    fig.clear()

    """
    Now plot Gamma for different densities.
    """
    num_rho = 10
    num_u = 1000

    rho_max = 5*rho0
    u_max   = 1e3*us2

    u = numpy.linspace(0.0, u_max, num_u)
    
    colors = iter(matplotlib.pyplot.cm.hot(numpy.linspace(0, 1, num_rho)))
    colors = iter(matplotlib.pyplot.cm.Set1(numpy.linspace(0, 1, num_rho)))
    
    i = 0
    for rho in numpy.linspace(rho0, rho_max, num_rho):
        #rho = rho0
        G = Gamma(rho, u)
        c = next(colors)
        #semilogx(u/us2, G, linewidth=1)
        print "i=", i, "rho=", rho
        semilogx(u/us2, G, color=c, linewidth=1)
        #plot(u/us2, G, color=c, linewidth=1)
        i += 1

    # Mark the limits
    plot([0.0, u_max/u0], [0.5, 0.5], '--', color='red')
    plot([0.0, u_max/u0], [1.8, 1.8], '--', color='red')

    xlim(0, u_max/us2)
    ylim(0.4, 1.9)

    xlabel(r'Normalized Internal energy (u/u$_{CV}$)')
    #xlabel(r'Internal energy ($\rho$/$\rho_0$)')
    ylabel(r'Generalized Gruneisen Parameter')

    legend(loc='best')

    savefig('tillotson-gamma-u.png', dpi=300, bbox_inches='tight')


    exit(0)

if __name__ == '__main__':
    main()

