/*
 * Copyright (c) 2018 Christian Reinhardt.
 * 
 * This file provides all the functions needed to implement an ideal gas EOS
 * into a hydro code.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "ig-eos.h"

/*
 * Basic functions:
 *
 * tillInit: initialise the library (contains all materials and converts to a given unit system).
 *
 * tillPressureSound: calculate the pressure and soundspeed for a given density and internal enegry using the Tillotson EOS.
 *
 * tillPressure: calculate the pressure for a given rho and u from the Tillotson EOS (uses tillPressureSound).
 *
 * tillSoundSpeed: calculate the sound speed for a given rho and u from the Tillotson EOS (uses tillPressureSound).
 *
 * tillFinalize: free memory.
 */


IGMAT *igInitMaterial(double dConstGamma, double dMeanMolMass, double b, double dKpcUnit, double dMsolUnit)
{
	/*
	 * Initialise a material
	 *
	 * We do:
	 * Initialize variables
	 * Convert quantities to code units
	 */
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    const double NA = 6.022e23;          /* Avogadro's number */
    IGMAT *material;
	int i;
	 
    material = calloc(1, sizeof(IGMAT));
    assert(material != NULL);

    /* There is only one material id. */
	material->iMaterial = IDEALGAS;

	/* This two parameters define the unit system we use */
    material->dKpcUnit = dKpcUnit;
    material->dMsolUnit = dMsolUnit;

	if (dKpcUnit <= 0.0 && dMsolUnit <= 0.0)
	{
		/* 
         * In this case units are not converted, so the code units are cgs.
         * Keep in mind that dGasConstant is not R but seems to be  kb/mp
         * (where mp is the mean particle mass) in Gasoline.
         */
		material->dGasConst = KBOLTZ/dMeanMolMass/MHYDR;
		material->dErgPerGmUnit = 1.0;
		material->dGmPerCcUnit = 1.0;
		material->dSecUnit = 1.0;
	} else {
		/*
		** Convert kboltz/mhydrogen to system units, assuming that
		** G == 1.
		*/
		material->dGasConst = material->dKpcUnit*KPCCM*KBOLTZ
		/MHYDR/GCGS/material->dMsolUnit/MSOLG;
		/* code energy per unit mass --> erg per g */
		material->dErgPerGmUnit = GCGS*material->dMsolUnit*MSOLG/(material->dKpcUnit*KPCCM);
		/* code density --> g per cc */
		material->dGmPerCcUnit = (material->dMsolUnit*MSOLG)/pow(material->dKpcUnit*KPCCM,3.0);
		/* code time --> seconds */
		material->dSecUnit = sqrt(1/(material->dGmPerCcUnit*GCGS));
	}

    material->dConstGamma = dConstGamma;
    material->dMeanMolMass = dMeanMolMass;

    /*
     * Add a finite volume to each gas particle to avoid crazy
     * densities at large pressure but neglect self-interaction.
     *
     * Hydrogen: b = 26.6 cm^3/mol (Wikipedia)
     *
     * The code uses b'=b/dMeanMolarMass instead, so be careful when
     * converting the parameter ([b] = cm^3/(g*mol)).
     *
     * For b=0 the ideal gas EOS is obtained.
     *
     * NOTE:    Introducing the parameter b defines a maximum density
     *          and the code must assert that rho < rho_max.
     */
    material->b = b/(material->dMeanMolMass*MHYDR*NA);

    /* cv = kb/mp */
    material->cv = KBOLTZ/((material->dConstGamma-1.0)*MHYDR*material->dMeanMolMass);
    material->rho0 = 0.001;

    /* Set a maximum density if b>0. */
    if (material-> b > 0.0) material->rhomax = 1.0/material->b;

    /*
     * Convert to code units!
     */
 	material->cv /= material->dErgPerGmUnit;
    material->b *=material->dGmPerCcUnit;

    return(material);
}

void tillFinalizeMaterial(TILLMATERIAL *material)
{
    /* Free the memory */
    if (material != NULL) free(material);
}

/*
 * This function prints the material constants for a given material.
 */
void igPrintMat(IGMAT *material)
{
    char MatName[256];

    assert(material != NULL);

    sprintf(MatName, "IDEAL_GAS");

    fprintf(stderr,"Material: %i (%s)\n", material->iMaterial, MatName);

    fprintf(stderr,"dConstGamma: %g\n", material->dConstGamma);
    fprintf(stderr,"dMeanMolMass: %g\n", material->dMeanMolMass);    
    fprintf(stderr,"rho0: %g\n", material->rho0);
    fprintf(stderr,"rho_max: %g\n", material->rho_max);
    fprintf(stderr,"cv: %g\n", material->cv);
}

/*
 * Calculate the pressure and sound speed.
 */
double igPressureSound(IGMAT *material, double rho, double u, double *pcSound)
{
    /*
     * (CR) 25.12.17: Generalized the ideal gas EOS by introducing a volume
     * to each gas particle. In the limit b=0 an ideal gas is obtained.
     *
     * NOTE: For densities larger than 1/b the pressure is negative. This
     *       must be treated properly or the code will crash.
     */
    if ((material->b > 0.0) && (rho >= (1.0/material->b)*0.99))
        rho = (1.0/material->b)*0.99;

    if (pcSound != NULL)
        *pcSound = material->dConstGamma*(material->dConstGamma-1.0)*u/pow(1.0-material->b*rho, 2.0);

    return ((material->dConstGamma-1.0)*rho*u/(1.0-material->b*rho));
}

double igPressure(IGMAT *material, double rho, double u)
{
    return (igPressureSound(material, rho, u, NULL));
}

double igdPdrho(IGMAT *material, double rho, double u)
{
    return ((material->dConstGamma-1.0)*u/(pow(1.0-material->b*rho, 2.0)));
}

double igdPdu(TILLMATERIAL *material, double rho, double u)
{
    return ((material->dConstGamma-1.0)*rho/(1.0-material->b*rho));
}

double igTempRhoU(TILLMATERIAL *material, double rho, double u)
{
    /*
     * Calculate T(rho,u) for a material. For an ideal gas the temperature is
     * independent of the density.
     */
    return(u/material->cv);
}

/*
 * Calculate rho(P, u) for an ideal gas EOS.
 */
double igRhoPU(TILLMATERIAL *material, double P, double u)
{
    return (P/((material->dConstGamma-1.0)*u+P*material->b));
}

/*
 * Calculate u(rho, P) for an ideal gas EOS.
 */
double igURhoP(TILLMATERIAL *material, double rho, double P)
{
    double a, b, c, Pa, Pb, Pc;

    // umin = 0.0
    a = 0.0;
    Pa = eosPressure(material, rho, a);

    b = material->vmax;
    Pb = eosPressure(material, rho, b);

    /*
     * Make sure the root is bracketed.
     */
    while (Pb <= P)
    {
        b = 2.0*b;
        Pb = eosPressure(material, rho, b);
    }

    fprintf(stderr, "rho= %g P= %g: a= %g, Pa= %g, b= %g, Pb= %g\n", rho, P, a, Pa, b, Pb);

    assert(Pa<P && Pb>P);

    while ((Pb-Pa) > 1e-10*Pc)
    {
        c = 0.5*(a + b);
        Pc = eosPressure(material, rho, c);

        if (Pc < P)
        {
            a = c;
            Pa = Pc;
        } else {
            b = c;
            Pb = Pc;
        }
    }

    return(c);
}



