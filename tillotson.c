/*
 * Copyright (c) 2014-2018 Christian Reinhardt and Joachim Stadel.
 *
 * This file provides all the functions for the Tillotson EOS library.
 * The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 * and convenient to use equation of state that can describe matter over
 * a large range of pressures, densities and internal energies.
 *
 * tillotson.c:
 *
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
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"

/* This will cut the pressure in the cold expanded states for rho/rho0 < 0.8 as suggested in Melosh1989. */
//#define TILL_PRESS_MELOSH

/*
 * Initialize a material from the Tillotson library
 *
 * We do:
 * Initialize variables
 * Convert quantities to code units
 * The memory for the look up table is allocated in tillInitLookup()
 */
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, int nTableRho, int nTableV, double rhomax, double vmax)
{
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    const double NA = 6.022e23;          /* Avogadro's number */
    TILLMATERIAL *material;
	int i;
	 
    material = calloc(1, sizeof(TILLMATERIAL));
    assert(material != NULL);

	material->iMaterial = iMaterial;

	/* This two parameters define the unit system we use. */
    material->dKpcUnit = dKpcUnit;
    material->dMsolUnit = dMsolUnit;

    /* Min and max values for the lookup table (in code units). */
	material->rhomax = rhomax;
	material->vmax = vmax;

    assert(material->rhomax > 0.0);
    assert(material->vmax > 0.0);

	if (dKpcUnit <= 0.0 && dMsolUnit <= 0.0)
	{
		/* In this case units are not converted, so the code units are cgs. */
		material->dGasConst = KBOLTZ;   // dGasConst is NOT KBOLZ !!!!!!!
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

	/*
	 * Set the Tillotson parameters for the material.
	 */
	switch(iMaterial)
	{
        case IDEALGAS:
            /*
             * Ideal gas EOS. Currently we are limited to monoatomic gases.
             * Keep in mind that dGasConstant is not R but seems to be kb/mp
             * (where mp is the mean particle mass) in Gasoline.
             */
            material->dConstGamma = 5.0/3.0;
            material->dMeanMolMass = 1.0;

//            material->dMeanMolMass = 23.0; // 10x solar value (mu=2.3)
//          material->dMeanMolMass = 11.5; // 5x solar value (mu=2.3)
//          material->dMeanMolMass = 17.25; // 7.5x solar value (mu=2.3)
#if 0
            /*
             * This doesnt work as cv is converted to code units below.
             */
            material->cv = material->dGasConst/((material->dConstGamma-1.0)*material->dMeanMolMass);
#endif
            // cv = kb/mp
            material->cv = KBOLTZ/((material->dConstGamma-1.0)*MHYDR*material->dMeanMolMass);
            material->rho0 = 0.001;

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
//            material->b = 26.6/(material->dMeanMolMass*MHYDR*NA);
            material->b = 0.0; 
            material->a = 0.0;
            fprintf(stderr, "Modified ideal gas: b=  %g [cm^3/g]\n", material->b);
            fprintf(stderr, "Modified ideal gas: mu= %g [1/m_H]\n", material->dMeanMolMass);
			break;
		case GRANITE:
			/*
			 * Material parameters from Benz et al. (1986).
			 */
			material->a = 0.5;
			material->b = 1.3;
			material->u0 = 1.6e11;		/* in ergs/g */
			material->rho0 = 2.7;		/* g/cc */
			material->A = 1.8e11;		/* ergs/cc */
			material->B = 1.8e11;		/* ergs/cc */
			material->us = 3.5e10;		/* ergs/g */
			material->us2 = 1.8e11;		/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.79e7;		/* ergs/g K (or 790 J/kg K) */ 
			break;
		case IRON:
			/*
			 * Material parameters from Benz et al. (1987).
			 */
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 9.5e10;		/* in ergs/g */
			material->rho0 = 7.86;		/* g/cc */
			material->A = 1.279e12;		/* ergs/cc */
			material->B = 1.05e12;		/* ergs/cc */
			material->us = 1.42e10;		/* ergs/g */
			material->us2 = 8.45e10;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.449e7;		/* ergs/g K */ 
			break;
		case BASALT:
			/*
			 ** Material parameters from Benz & Asphaug (1999).
			 */
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 4.87e12;		/* in ergs/g */
			material->rho0 = 2.7;		/* g/cc */
			material->A = 2.67e11;		/* ergs/cc */
			material->B = 2.67e11;		/* ergs/cc */
			material->us = 4.72e10;		/* ergs/g */
			material->us2 = 1.82e11;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.84e7;		/* ergs/g K */ 
			break;
		case ICE:
			/*
			 * Material parameters from Benz & Asphaug (1999).
			 */
			material->a = 0.3;
			material->b = 0.1;
			material->u0 = 1.0e11;		/* in ergs/g */
			material->rho0 = 0.917;		/* g/cc */
			material->A = 9.47e10;		/* ergs/cc */
			material->B = 9.47e10;		/* ergs/cc */
			material->us = 7.73e9;		/* ergs/g */
			material->us2 = 3.04e10;	/* ergs/g */
			/*
			 * Note that the expression for Pe provided in Benz et al. (1986) does not agree
			 * in alpha and beta with the one provided in Melosh (1989) or Benz & Asphaug (1999).
			 */
			material->alpha = 10.0;
			material->beta = 5.0;
			// Have to look up in more details?
			material->cv = 1.0e7;		/* ergs/g K */ 
			break;
		case WATER:
			/*
			 * Material parameters from Woolfson (2007).
			 */
			material->a = 0.5;
			material->b = 0.9;
			material->u0 = 2.0e10;		/* in ergs/g */
			material->rho0 = 1.00;		/* g/cc */
			material->A = 2.00e11;		/* ergs/cc */
			material->B = 1.00e11;		/* ergs/cc */
			material->us = 4.00e9;		/* ergs/g */
			material->us2 = 2.04e10;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			// Have to look up in more details?
			material->cv = 4.1814e7;	/* ergs/g K */ 
			break;
		case DUNITE:
			/*
			 * Material parameters that imitate dunite in ANEOS (Thompson 1972).
			 */
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 4.87e12;		/* in ergs/g */
            // Key difference is the reference density
			material->rho0 = 3.32;		/* g/cc */
			material->A = 2.67e11;		/* ergs/cc */
			material->B = 2.67e11;		/* ergs/cc */
			material->us = 4.72e10;		/* ergs/g */
			material->us2 = 1.82e11;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.84e7;		/* ergs/g K */ 
			break;
		default:
			/* Unknown material */
			assert(0);
	}

    /*
     * Convert energies and densities to code units!
     */
    material->u0 /= material->dErgPerGmUnit;
    material->us /= material->dErgPerGmUnit;
    material->us2 /= material->dErgPerGmUnit;
    material->rho0 /= material->dGmPerCcUnit;
    material->A /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
    material->B /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
 
 	material->cv /= material->dErgPerGmUnit;

    /*
     * In case of an ideal gas the parameters a and b have to be converted to
     * code units too. Note that [b] = cm^3/(g*mol).
     */
    if (iMaterial == IDEALGAS)
    {
        material->b *=material->dGmPerCcUnit;
//        fprintf(stderr, "b= %g [RE^3/Munit]\n", material->b);
    }

#if 0
    if (material->iMaterial == IDEALGAS)
    {
        fprintf(stderr, "Ideal gas: cv= %g\n in code units.\n", material->cv);
    }
#endif

	/* Number of grid points for the look up table. */
	material->nTableRho = nTableRho;
	material->nTableV = nTableV;

    /* The stuff below does not have to be done for the ideal gas as there is no lookup table. */
    if (material->iMaterial == IDEALGAS)
    {
        material->rhomin = 0.0;
        material->n = 0;
        /* rhomax is set already. */
        material->dlogrho = 0.0;
    } else {
        /* Set rhomin */
        material->rhomin = TILL_RHO_MIN;

        /* rhomin has to be larger than zero otherwise the logarithmic spacing does not work. */
        assert(material->rhomin > 0.0 && material->rhomin < material->rhomax);
        
        /* Set dlogrho so that log(rho0) lies on the grid. */
        material->n = floor((log(material->rho0)-log(material->rhomin))/(log(material->rhomax)-log(material->rhomin))*material->nTableRho);
        material->dlogrho = (log(material->rho0)-log(material->rhomin))/material->n;

        /* Set the actual rhomax (note that rhomax can differ more after the correction when using log(rho) as variable. */ 
        material->rhomax = tillLookupRho(material, material->nTableRho-1);
    }

//#ifdef TILL_VERBOSE
    fprintf(stderr, "tillInitialize: iMat= %i n= %i dlogrho= %g dv= %g.\n", material->iMaterial, material->n, material->dlogrho, material->dv);
    fprintf(stderr, "tillInitialize: nTableRho= %i nTableV= %i.\n", material->nTableRho, material->nTableV);
//#endif


    /* This is the same for all materials. */
    material->dv = material->vmax/(material->nTableV-1);

    return material;
}

void tillFinalizeMaterial(TILLMATERIAL *material)
{
    /* Free the memory */
    if (material->Lookup != NULL) free(material->Lookup);
    //	if (material->cold != NULL) free(material->cold);

    free(material);
}

/*
 * This function stores the materials name in a string string.
 */
void tilliMatString(TILLMATERIAL *material, char *MatName)
{
    assert(MatName != NULL);

    switch(material->iMaterial)
    {
        case IDEALGAS:
            sprintf(MatName, "IDEAL_GAS");
            break;
        case GRANITE:
            sprintf(MatName, "GRANITE");
            break;
        case IRON:
            sprintf(MatName, "IRON");
            break;
        case BASALT:
            sprintf(MatName, "BASALT");
            break;
        case ICE:
            sprintf(MatName, "ICE");
            break;
        case WATER:
            sprintf(MatName, "WATER");
            break;
        default:
            /* Unknown material */
            assert(0);
    } 
}

/*
 * This function returns an error message for each error code.
 */
void tillErrorString(int iError, char *ErrorString)
{
    assert(ErrorString != NULL);

    switch(iError)
    {
        case TILL_LOOKUP_SUCCESS:
            sprintf(ErrorString, "TILL_LOOKUP_SUCCESS");
            break;
        case TILL_LOOKUP_OUTSIDE_RHOMIN:
            sprintf(ErrorString, "TILL_LOOKUP_OUTSIDE_RHOMIN");
            break;
        case TILL_LOOKUP_OUTSIDE_RHOMAX:
            sprintf(ErrorString, "TILL_LOOKUP_OUTSIDE_RHOMAX");
            break;
        case TILL_LOOKUP_OUTSIDE_VMIN:
            sprintf(ErrorString, "TILL_LOOKUP_OUTSIDE_VMIN");
            break;
        case TILL_LOOKUP_OUTSIDE_VMAX:
            sprintf(ErrorString, "TILL_LOOKUP_OUTSIDE_VMAX");
            break;
        default:
            sprintf(ErrorString, "UNKNOWN ERROR");
    } 
}

/*
 * This function prints the material constants for a given material.
 */
void tillPrintMat(TILLMATERIAL *material)
{
    char MatName[256];

    assert(material != NULL);
    tilliMatString(material, MatName);

    fprintf(stderr,"Material: %i (%s)\n", material->iMaterial, MatName);

    /*
     * Currently the ideal gas is treated differently.
     */
    if (material->iMaterial == IDEALGAS)
    {
        fprintf(stderr,"dConstGamma: %g\n", material->dConstGamma);
        fprintf(stderr,"dMeanMolMass: %g\n", material->dMeanMolMass);    
        fprintf(stderr,"rho0: %g\n", material->rho0);
        fprintf(stderr,"cv: %g\n", material->cv);
    } else {
        fprintf(stderr,"a: %g\n", material->a);
        fprintf(stderr,"b: %g\n", material->b);
        fprintf(stderr,"A: %g\n", material->A);
        fprintf(stderr,"B: %g\n", material->B);
        fprintf(stderr,"rho0: %g\n", material->rho0);
        fprintf(stderr,"u0: %g\n", material->u0);
        fprintf(stderr,"us: %g\n", material->us);
        fprintf(stderr,"us2: %g\n", material->us2);
        fprintf(stderr,"alpha: %g\n", material->alpha);
        fprintf(stderr,"beta: %g\n", material->beta);
        fprintf(stderr,"cv: %g\n", material->cv);
    }
}

/*
 * These functions provide a more general interface for EOS calls
 * so the user can in principle add different EOS (e.g., an ideal
 * gas).
 */
double eosPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
    if (material->iMaterial == IDEALGAS)
    {
        /*
         * In this case an ideal gas EOS is used.
         */
#if 0
        if (pcSound != NULL) *pcSound = material->dConstGamma*(material->dConstGamma-1.0)*u;
        return ((material->dConstGamma-1.0)*rho*u);
#endif
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
    } else {
        return (tillPressureSound(material, rho, u, pcSound));
    }
}

double eosPressure(TILLMATERIAL *material, double rho, double u)
{
    return (eosPressureSound(material, rho, u, NULL));
}

double eosdPdrho(TILLMATERIAL *material, double rho, double u)
{
    if (material->iMaterial == IDEALGAS)
    {
        /*
         * In this case an ideal gas EOS is used.
         */
#if 0
        return ((material->dConstGamma-1.0)*u);
#endif
        /*
         * (CR) 25.12.17: Generalized the ideal gas EOS by introducing a volume
         * to each gas particle. In the limit b=0 an ideal gas is obtained.
         * (CR) 20.01.18: Found a bug in the expression for the ideal gas!
         */
        return ((material->dConstGamma-1.0)*u/(pow(1.0-material->b*rho, 2.0)));
    } else {
        return (tilldPdrho(material, rho, u));
    }
}

double eosdPdu(TILLMATERIAL *material, double rho, double u)
{
    if (material->iMaterial == IDEALGAS)
    {
        /*
         * In this case an ideal gas EOS is used.
         */
# if 0
        return ((material->dConstGamma-1.0)*rho);
#endif
        /*
         * (CR) 25.12.17: Generalized the ideal gas EOS by introducing a volume
         * to each gas particle. In the limit b=0 an ideal gas is obtained.
         */
        return ((material->dConstGamma-1.0)*rho/(1.0-material->b*rho));
    } else {
        return (tilldPdu(material, rho, u));
    }
}

double eosTempRhoU(TILLMATERIAL *material, double rho, double u)
{
    /*
     * Calculate T(rho,u) for a material. As an approximation we use
     *
     * u(rho,T) = uc(rho) + cv*T
     *
     * for condensed materials.
     */
    assert(material->cv > 0.0);

    /*
     * For an ideal gas the temperature is independent of the density.
     */
    if (material->iMaterial == IDEALGAS)
    {
        //        fprintf(stderr, "(CR) eosTempRhoU(): Ideal gas: rho= %g u= %g gamma= %g cv= %g T= %g\n", rho, u, material->dConstGamma, material->cv, u/material->cv);
        return(u/material->cv);
    } else {
        return(tillTempRhoU(material, rho, u));
    }
}

/*
 * Calculate rho(P, u) for a given EOS and material using bisection.
 */
double eosRhoPU(TILLMATERIAL *material, double P, double u)
{
    double a, b, c, Pa, Pb, Pc;

    a = 0.0;
    Pa = eosPressure(material, a, u);

    b = material->rhomax;
    Pb = eosPressure(material, b, u);

    /*
     * Make sure the root is bracketed.
     */
    while (Pb < P)
    {
        b = 2.0*b;
        Pb = eosPressure(material, b, u);
    }

    fprintf(stderr, "a= %g, Pa= %g, b= %g, Pb= %g\n", a, Pa, b, Pb);

    assert(Pa < P && Pb > P);

    while ((Pb-Pa) > 1e-10*Pc)
    {
        c = 0.5*(a + b);
        Pc = eosPressure(material, c, u);

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

/*
 * Calculate u(rho, P) for a given EOS and material using bisection.
 */
double eosURhoP(TILLMATERIAL *material, double rho, double P)
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

/*
 * Calculate rho(P, T) for a given EOS and material.
 */
double eosRhoPTemp(TILLMATERIAL *material, double P, double T)
{
    double rho;
    
    if (material->iMaterial == IDEALGAS)
    {
        // For the ideal gas there is an analytic expression.
        rho = P/((material->dConstGamma-1.0)*material->cv*T);
    } else {
        // If it is a Tillotson material call tillRhoPTemp().
        rho = tillRhoPTemp(material, P, T);
    }
    
    return rho;
}

/*
 * Calculate phi and gamma as given in Hu et al. (2009).
 */
double eosPhi(TILLMATERIAL *material, double rho, double u)
{
    return(eosdPdrho(material, rho, u));
}

double eosGamma(TILLMATERIAL *material, double rho, double u)
{
    return(1.0/rho*eosdPdu(material, rho, u));
}

double tilldPdrho_s(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate dP/drho at constant entropy.
     */

    return (1.0/(rho*rho)*(tillSoundSpeed2old(material,rho, u)-2.0*tillPressure(material,rho,u)/rho));
}

double tillSoundSpeed2old(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate the sound speed for the Tillotson EOS like we did in
     ** Gasoline. In the intermediate states its better however to do a
     ** linear interpolation in the sound speed because we have no
     ** discontinuity when we change from expanded hot to condensed states.
     */

    double eta, mu;
    double Pc, Pe;
    double c2c, c2e;
    double Gammac, Gammae, w0, y, z;

    eta = rho/material->rho0;
    mu = eta - 1.0;
    z = (1.0 - eta)/eta;
    w0 = u/(material->u0*eta*eta)+1.0;

    /*
     **  Here we evaluate, which part of the equation of state we need.
     */
    if (rho >= material->rho0) {
        /*
         **  condensed states (rho > rho0)
         */
        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;

        return ((Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho));
    } else if (u < material->us) {
        /* 
         ** cold expanded states (rho < rho0 and u < us)
         ** P is like for the condensed states
         */
        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;

        return ((Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho));
    } else if (u > material->us2) {
        /*
         ** expanded hot states (rho < rho0 and u > us2)
         */
        Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
        Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

        return ((Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho)));
    } else {
        /*
         **  intermediate states (rho < rho0 and us < u < us2)
         */
        y = (u - material->us)/(material->us2 - material->us);

        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
        Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
        Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

        return ((Gammac*(1.0-y)+Gammae*y+1.0)*(Pc*(1.0-y)+Pe*y)/rho+((material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho))*(1.0-y)+(material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho)))*y);
    }
}

double tillPressureSoundold(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
    /*
     ** Calculate the pressure and sound speed from the Tillotson EOS for a material.
     ** Set pcSound = NULL to only calculate the pressure forces. Here we used the old
     ** function that was originally implemented in Gasoline.
     */

    double eta, mu;
    double P,c2;
    double Gamma, epsilon0, y, x;

    eta = rho/material->rho0;
    mu = eta - 1.0;
    x  = 1.0 - eta;

    /*
     **  Here we evaluate, which part of the equation of state we need.
     */
    if (rho >= material->rho0) {
        /*
         **  condensed states (rho > rho0)
         */
        y = 0.0;
    } else if (u < material->us) {
        /* 
         ** cold expanded states (rho < rho0 and u < us)
         ** P is like for the condensed states
         */
        y = 0.0;
    } else if (u > material->us2) {
        /*
         ** expanded hot states (rho < rho0 and u > us2)
         */
        y = 1.0;
    } else {
        /*
         **  intermediate states (rho < rho0 and us < u < us2)
         */
        y = (u - material->us)/(material->us2 - material->us);
    }

    Gamma = material->a + material->b/(u/(material->u0*eta*eta)+1.0)*((1.0-y) + y*exp(-material->beta*x*x/eta/eta));
    epsilon0 = 1.0/(Gamma*rho)*((material->A*mu + material->B*mu*mu)*(1.0-y)+y*material->A*mu*exp(-(material->alpha*x/eta + material->beta*x*x/eta/eta)));
    P = Gamma*rho*(u + epsilon0);

    if (pcSound != NULL)
    {
        /* calculate the sound speed */
        c2 = (Gamma+1.0)*P/rho + 1.0/rho*((material->A+material->B*(eta*eta-1.0))*(1.0-y)+y*material->A*exp(-(material->alpha*x/eta+material->beta*x*x/eta/eta))*(1+mu/eta*(material->alpha+2*material->beta*x/eta)));
        /* make sure that c^2 > 0 for rho < rho0 */
        if (c2 < material->A/material->rho0){
            c2 = material->A/material->rho0;
        }

        *pcSound = sqrt(c2);
    }

    return (P);
}

double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
    /*
     ** Calculate the pressure and sound speed from the Tillotson EOS for a material.
     ** Set pcSound = NULL to only calculate the pressure forces.
     */

    double eta, mu;
    double Pc, Pe;
    double c2c, c2e;
    double Gammac, Gammae, w0, y, z;

    eta = rho/material->rho0;
    mu = eta - 1.0;
    z = (1.0 - eta)/eta;
    w0 = u/(material->u0*eta*eta)+1.0;

    /*
     **  Here we evaluate, which part of the equation of state we need.
     */
    if (rho >= material->rho0 || u < material->us) {
        /*
         **  Condensed states (rho > rho0) or expanded cold states.
         */
        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;

#ifdef TILL_PRESS_MELOSH
        /* Melosh 1989 suggests to cut the pressure in the expanded cold states for rho/rho < 0.8 */
        if (eta < 0.8)
        {
            //			fprintf(stderr,"Setting pressure to zero for eta=%g\n",eta);
            Pc = 0.0;
        }
#endif
#ifdef TILL_PRESS_NP
        /* The pressure is set to zero if it becomes negative. */
        if (Pc < 0.0) Pc = 0.0;
#endif
        if (pcSound != NULL)
        {
            /* Calculate the sound speed. */
            //			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
            // Altough for P > 0 the two expressions are the same, for P=0 this gives a bigger sound speed.
            c2c =material->a*u+material->b*u/(w0*w0)*(3.0*w0-2.0)+(material->A+2.0*material->B*mu)/material->rho0 + Pc/(rho*rho)*(material->a*rho+material->b*rho/(w0*w0));
            /*
             * Set the minimum sound speed to the uncompressed bulk sound speed to avoid issues when P<0.
             */
            c2c = MAX(c2c, material->A/material->rho0);
            //*pcSound = sqrt(c2c);
            *pcSound = c2c;
        }
        return (Pc);
    } else if (u > material->us2) {
        /*
         ** Expanded hot states (rho < rho0 and u > us2).
         */
        Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
        Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

        if (pcSound != NULL)
        {
            /* calculate the sound speed */
            c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));
            //*pcSound = sqrt(c2e);
            c2e = MAX(c2e, material->A/material->rho0);
            *pcSound = c2e;
        }

        return (Pe);
    } else {
        /*
         **  intermediate states (rho < rho0 and us < u < us2)
         */
        y = (u - material->us)/(material->us2 - material->us);

        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
        Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
        Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

#ifdef TILL_PRESS_MELOSH
        /* Melosh 1989 suggests to cut the pressure in the expanded cold states for rho/rho < 0.8 */
        if (eta < 0.8)
        {
            //			fprintf(stderr,"Setting pressure to zero for eta=%g\n",eta);
            Pc = 0.0;
        }
#endif
#ifdef TILL_PRESS_NP
        /* The pressure is set to zero if it becomes negative. */
        if (Pc < 0.0) Pc = 0.0;
#endif

        if (pcSound != NULL)
        {
            /* calculate the sound speed */
            //			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
            c2c =material->a*u+material->b*u/(w0*w0)*(3.0*w0-2.0)+(material->A+2.0*material->B*mu)/material->rho0 + Pc/(rho*rho)*(material->a*rho+material->b*rho/(w0*w0));
            c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));

            //*pcSound = sqrt(c2c*(1.0-y)+c2e*y);
            *pcSound = c2c*(1.0-y)+c2e*y;
            *pcSound = MAX(*pcSound, material->A/material->rho0);
        }

        return (Pc*(1.0-y)+Pe*y);
    }
}

double tillPressure(TILLMATERIAL *material, double rho, double u)
{
    /* Calculate the pressure from the Tillotson EOS for a material */
    double P = tillPressureSound(material, rho, u, NULL);

#ifdef TILL_PRESS_NP
    /* Make a pressure cut off, if P < 0. */
    if (P < 0.0 ) P = 0.0;
#endif

    return (P);
}

double tillPressureNP(TILLMATERIAL *material, double rho, double u)
{
    /* Calculate the pressure from the Tillotson EOS for a material */
    double P = tillPressureSound(material, rho, u, NULL);

    return (P);
}

double tilldPdrho(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate dP/drho at u=const.
     */
    double eta, mu;
    double dPcdrho, dPedrho;
    double w0, y, z;

    eta = rho/material->rho0;
    mu = eta - 1.0;
    z = (1.0 - eta)/eta;
    w0 = u/(material->u0*eta*eta)+1.0;

    /*
     ** Here we evaluate, which part of the equation of state we need.
     */
    if (rho >= material->rho0 || u < material->us) {
        /*
         ** Condensed (rho > rho0) or cold expanded states (rho < rho0 and u < us).
         */
        dPcdrho = material->a*u + material->b*u/(w0*w0)*(3.0*w0-2.0) + (material->A+2.0*material->B*mu)/material->rho0;

        return (dPcdrho);
    } else if (u > material->us2) {
        /*
         ** Expanded hot states (rho < rho0 and u > us2).
         */
        dPedrho = (material->a + material->b/w0*exp(-material->beta*z*z)*(2.0*material->beta*z/eta+3.0-2.0/w0))*u+material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta));		
        return (dPedrho);
    } else {
        /*
         **  Intermediate states (rho < rho0 and us < u < us2).
         */
        y = (u - material->us)/(material->us2 - material->us);

        dPcdrho = material->a*u + material->b*u/(w0*w0)*(3.0*w0-2.0) + (material->A+2.0*material->B*mu)/material->rho0;
        dPedrho = (material->a + material->b/w0*exp(-material->beta*z*z)*(2.0*material->beta*z/eta+3.0-2.0/w0))*u+material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta));		
        return (dPcdrho*(1.0-y)+dPedrho*y);
    }
}

double tilldPdu(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate dP/du at rho=const.
     */
    double eta, mu;
    double dPcdu, dPedu;
    double w0, y, z;

    eta = rho/material->rho0;
    mu = eta - 1.0;
    z = (1.0 - eta)/eta;
    w0 = u/(material->u0*eta*eta)+1.0;

    /*
     **  Here we evaluate, which part of the equation of state we need.
     */
    if (rho >= material->rho0 || u < material->us) {
        /*
         ** Condensed (rho > rho0) or cold expanded states (rho < rho0 and u < us).
         */
        dPcdu = material->a*rho + material->b*rho/(w0*w0);

        return (dPcdu);
    } else if (u > material->us2) {
        /*
         ** Expanded hot states (rho < rho0 and u > us2).
         */
        dPedu = (material->a + material->b/(w0*w0)*exp(-material->beta*z*z))*rho;		
        return (dPedu);
    } else {
        /*
         ** Intermediate states (rho < rho0 and us < u < us2).
         */
        y = (u - material->us)/(material->us2 - material->us);

        dPcdu = material->a*rho + material->b*rho/(w0*w0);
        dPedu = (material->a + material->b/(w0*w0)*exp(-material->beta*z*z))*rho;		

        return (dPcdu*(1.0-y)+dPedu*y);
    }
}

double tilldTdrho(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate dT/drho at u=const.
     */
    assert(material->cv > 0.0);
    return (-1.0/material->cv*tillPressure(material,rho,tillColdULookup(material,rho))*(rho*rho));
}

double tilldTdu(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate dT/du at rho=const.
     */

    assert(material->cv > 0.0);
    return (-1.0/material->cv);
}

double tillPressureRhoU(TILLMATERIAL material, double rho, double u)
{
    /* Calculate the pressure from the Tillotson EOS for a material */
    assert(0);
}

double tillTempRhoU(TILLMATERIAL *material, double rho, double u)
{
    /*
     ** Calculate T(rho,u) for a material. As an approximation
     ** we use u(rho,T) = uc(rho) + cv*T.
     */
    assert(material->cv > 0.0);
    return ((u-tillColdULookup(material,rho))/material->cv);
}

double tillTempRhoP(TILLMATERIAL *material, double rho, double P)
{
    /* Calculate T(rho,P) for a material */
    assert(0);
}

double tillURhoTemp(TILLMATERIAL *material, double rho, double T)
{
    /* Calculate u(rho,T) for a material */
    return(tillColdULookup(material, rho) + material->cv*T);
}

double tillRhoPTemp(TILLMATERIAL *material, double P, double T)
{
    /* 
     * Calculate rho(P,T) for a material. Because thermodynamical consistency requires
     * (dP/drho)_T > 0 this should always work except where we do the pressure cutoff.
     */
    double a, ua, Pa, b, ub, Pb, c, uc, Pc;

    Pc = 0.0;

    /*
     * Check, if P <= 0.0. In this case the algorithm will not work and the
     * function returns a negative value for the density to indicate that it
     * failed.
     */
    if (P <= 0.0)
    {
#if TILL_VERBOSE
        fprintf(stderr, "tillRhoPTemp: Called for negative pressure (P= %g, Pa= %g, Pb= %g). No correction is done.\n", P);
#endif
        return -1.0;
    }

    /*
     * We use rhoa=0 and rhob=rhomax as limits.
     */
    a = material->rhomin;

    // For rho=0.0 the pressure diverges so set a minimum density.
    if (a < 1e-5)
        a = 1e-5;

    // We should find a way to more precisely bracket P and avoid Pa <= 0.
    ua = tillURhoTemp(material, a, T);
    Pa = tillPressure(material, a, ua);

    b = material->rhomax*0.99;
    ub = tillURhoTemp(material, b, T);
    Pb = tillPressure(material, b, ub);

    // Make sure that Pb > P. If nescessary the lookup table has to be expanded.
    while (Pb < P)
    {
        /*
         * Careful, here a new lookup table has to be generated. Currently the new parameters are
         * set by hand but it would make more sense to move the initialization of all lookup table
         * related values in tillInitLookup().
         */
        material->rhomax *= 2.0;
        
        /* Set drho so that rho0 lies on the grid. */
        material->n = floor((material->rho0-material->rhomin)/(material->rhomax-material->rhomin)*material->nTableRho);
        material->drho =  (material->rho0-material->rhomin)/material->n;

        /* Set the actual rhomax. */ 
        material->rhomax = material->drho*(material->nTableRho-1);
        
        tillInitLookup(material);

#ifdef TILL_VERBOSE
        fprintf(stderr, "tillRhoPTemp: P > Pb, expanding lookup table (P= %g, Pb= %g, b= %g).\n", P, Pb, b);
#endif
        b = material->rhomax*0.99;

        ub = tillURhoTemp(material, b, T);
        Pb = tillPressure(material, b, ub);
    }

    /* What do we do for P=0 in the expanded cold states?*/
    //	fprintf(stderr,"tillRhoPTemp: starting with a= %g ua= %g Pa= %g b=%g ub= %g Pb= %g\n", a, ua, Pa, b, ub, Pb);
    /* Check if the root is bracketed. */
    if ((Pa >= P) || (Pb <=P))
    {
#ifdef TILL_VERBOSE
        fprintf(stderr, "tillRhoPTemp: Root can not be bracketed (P= %g, Pa= %g, Pb= %g).\n", P, Pa, Pb);
#endif
        return -1.0;
    }

    assert (Pa < P && P < Pb);	

    /*
     * Root bracketed by (a,b).
     */
    while (Pb-Pa > 1e-10) {
        c = 0.5*(a + b);
        uc = tillURhoTemp(material,c,T);
        Pc = tillPressure(material,c, uc);

        if (Pc < P) {
            a = c;
            Pa = Pc;
        }
        else {
            b = c;
            Pb = Pc;
        }
    }

    //fprintf(stderr,"tillRhoPTemp: rhoc= %g uc= %g Pc= %g P= %g T= %g\n", c, uc, Pc, P, T);
    /*
     * Return values.
     */
    return c;
}

double tillSoundSpeed(TILLMATERIAL *material, double rho, double u)
{
    /* Calculate sound speed for a material */
    double c;

    tillPressureSound(material, rho, u, &c);
    return (c);
}

double tillDensRatio(TILLMATERIAL material1, TILLMATERIAL material2, double P, double T)
{
    /* From Woolfson 2007 */
    assert(0);
}

/* 
 * Calculate rho(P,u) by doing a root finding using bisection.
 */
double tillRhoPU(TILLMATERIAL *material, double P, double u)
{
    double a, b, c, Pa, Pb, Pc;

    /*
     * Try to bracket the root with (a, b).
     */
    if (u == 0.0)
    {
        // For rho=0 AND u=0 Pc diverges!
        a = 1e-30;
    } else {
        // If u > 0 there is no problem
        a = 0.0;
    }

    Pa = tillPressure(material, a, u);
    b = 10.0*material->rho0;
    Pb = tillPressure(material, b, u);

    while (Pb < P)
    {
        b *= 2.0;
        Pb = tillPressure(material, b, u);
    }

    //    fprintf(stderr, "tillRhoPU: Root bracketed by: a=%15.7E Pa=%15.7E b=%15.7E Pb=%15.7E\n", a, Pa, b, Pb);
    //    assert(Pa < P && Pb > P);

    if (Pa >= P) return(a);
    if (Pb <= P) return(b);

    while (Pb-Pa > 1e-10*Pc)
    {
        c = 0.5*(a+b);
        Pc = tillPressure(material, c, u);
        //        fprintf(stderr, "a= %15.7E ua= %15.7E b= %15.7E ub= %15.7E c= %15.7E uc= %15.7E\n",a,ua,b,ub,c,uc);
        if (Pc < P)
        {
            a = c;
            Pa = Pc;
        } else {
            b = c;
            Pb = Pc;
        }
    }

    return (c);
}

double tilldudrho(TILLMATERIAL *material, double rho, double u)
{
    return (tillPressure(material,rho,u)/(rho*rho));
}

/*
 * Calculate the logarithmic derivative of u with respect to rho.
 *
 * du/dlogrho = P/rho
 *
 */
double tilldudlogrho(TILLMATERIAL *material, double logrho, double u)
{
    double rho = exp(logrho);
    return (tillPressure(material,rho,u)/(rho));
}

int tillSolveBC(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double rho1, double u1, double *prho2, double *pu2)
{
    /*
     * Given rho1, u1 (material 1) solve for rho2, u2 (material 2) at the interface between two material.
     *
     * The b.c. are:
     *
     * P1(rho1,u1)=P2(rho2,u2) and T1(rho1,u1)=T2(rho2,u2)
     *
     * Since P = P(rho, T) and T1=T2 we solve for P2(rho2)-P1=0.
     *
     * Returns 0 if successful or -1 if not.
     */
    double P, T;
    double a, ua, Pa, b, ub, Pb, c, uc, Pc;
    int iRet;

    iRet = -1;
    Pc = 0.0;

    /* Calculate P and T in material 1. */
    P = tillPressure(mat1, rho1, u1);
    T = tillTempRhoU(mat1, rho1, u1);

    fprintf(stderr,"modelSolveBC: P= %g, T= %g\n", P, T);

    /*
     * We use rho1 as an upper limit for rho2 assuming that the denser component is in the inner shell.
     */
    a = rho1;
    ua = tillURhoTemp(mat2, a, T);
    Pa = tillPressure(mat2, a, ua);

    b = mat2->rhomin;
    ub = tillURhoTemp(mat2, b, T);
    Pb = tillPressure(mat2, b, ub);

    fprintf(stderr,"modelSolveBC: starting with a=%g ua=%g Pa=%g b=%g ub=%g Pb=%g\n",a ,ua, Pa, b, ub, Pb);

    /*
     * Assert that the root is bracketed by (a, b).
     */
    if (Pa < P || Pb > P)
    {

        return(iRet);

    }

    //	assert (Pa > P && Pb < P);	

    /*
     * Root bracketed by (a,b).
     */
    while (Pa-Pb > 1e-10) {
        c = 0.5*(a + b);
        uc = tillURhoTemp(mat2,c,T);
        Pc = tillPressure(mat2,c, uc);

        if (Pc < P) {
            b = c;
            Pb = Pc;
        }
        else {
            a = c;
            Pa = Pc;
        }
        //		fprintf(stderr,"c:%.10g Pc:%.10g\n",c,Pc);
    }

    /*
     * Return values.
     */
    *prho2 = c;
    *pu2 = uc; 

    iRet = 0;

    return(iRet);
}

