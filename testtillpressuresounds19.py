#!/usr/bin/env python3
"""
Show how the two different ways to calculate the pressure and sound speed differ.
"""
from matplotlib import *
from matplotlib.pyplot import *
import numpy

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

    # Restore classic font used for math
    rcParams['mathtext.fontset'] = 'cm'
    rcParams['mathtext.rm']      = 'serif'

    # Set Up Figure, Single Column MNRAS
    fig = gcf()
    ax = gca()
    fig, ax = subplots(1,1)
    fig.set_size_inches(8.27*0.39,8.27*(6./8.)*0.39)

    ErgPerGmUnit = 9998228982.69
    GmPerCcUnit  = 0.368477421278

    """
    Material constants for granite from Benz et al. (1986).
    """
    rho0 = 2.7/GmPerCcUnit;
    us   = 3.5e10/ErgPerGmUnit
    us2  = 1.8e11/ErgPerGmUnit 
    A    = 1.8e11/ErgPerGmUnit

    cs_min = (A/rho0)**0.5

    """
    Load the rho and u axis.
    """
    rho = numpy.loadtxt("testtillpressuresounds19_rhoaxis.txt")
    u   = numpy.loadtxt("testtillpressuresounds19_uaxis.txt")

    rho_min = numpy.min(rho)
    rho_max = numpy.max(rho)

    u_min = numpy.min(u)
    u_max = numpy.max(u)

    """
    Plot the relative error in the pressure.
    """
    data = numpy.loadtxt("testtillpressures19.txt")

    """
    data1 = numpy.zeros(data.shape)

    data1[numpy.where(numpy.abs(data) < 1e-2)] = 4
    data1[numpy.where((numpy.abs(data) >= 1e-2) & (numpy.abs(data) < 1e-1))] = 3
    data1[numpy.where((numpy.abs(data) >= 1e-1) & (numpy.abs(data) < 0.5))] = 2
    data1[numpy.where((numpy.abs(data) >= 0.5) & (numpy.abs(data) < 1e1))] = 1

    print numpy.where(data1 <= 0.0)

    imshow(data1)
    colorbar()
    """
    imshow(data, origin='lower', extent=[rho_min, rho_max, u_min, u_max])
    colorbar()

    title("Pressure")
    xlabel("Density")
    ylabel("Internal energy")

    plot([rho0, rho0], [u_min, u_max], '--', color='red')
    plot([rho_min, rho0], [us, us], '--', color='red')
    plot([rho_min, rho0], [us2, us2], '--', color='red')

    savefig('testtillpressures19.png', dpi=300, bbox_inches='tight')
    show()

    fig = gcf()
    fig.clear()

    """
    Plot the relative error in the sound speed.
    """
    data = numpy.loadtxt("testtillsounds19.txt")

    imshow(data, origin='lower', extent=[rho_min, rho_max, u_min, u_max])
    colorbar()

    title("Sound speed")
    xlabel("Density")
    ylabel("Internal energy")

    plot([rho0, rho0], [u_min, u_max], '--', color='red')
    plot([rho_min, rho0], [us, us], '--', color='red')
    plot([rho_min, rho0], [us2, us2], '--', color='red')

    savefig('testtillsounds19.png', dpi=300, bbox_inches='tight')
    show()

    fig = gcf()
    fig.clear()

    # This doesnt work yet.
    """
    Plot where the minimum sound speed is set.
    """
    data = numpy.loadtxt("testtillsounds19.txt")

    data1 = numpy.zeros(data.shape)

    data1[numpy.where(data < cs_min)] = 1
    data1[numpy.where(data > cs_min)] = 2

    imshow(data1, origin='lower', extent=[rho_min, rho_max, u_min, u_max])
    colorbar()

    title("Min. Sound speed")
    xlabel("Density")
    ylabel("Internal energy")

    plot([rho0, rho0], [u_min, u_max], '--', color='red')
    plot([rho_min, rho0], [us, us], '--', color='red')
    plot([rho_min, rho0], [us2, us2], '--', color='red')

    savefig('testtillsounds19.png', dpi=300, bbox_inches='tight')
    show()

    exit(0)

if __name__ == '__main__':
    main()

