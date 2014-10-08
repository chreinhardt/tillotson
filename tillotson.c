/*
 ** This file provides all the functions for the Tillotson EOS library.
 ** The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 ** and convenient to use equation of state that can describe matter over
 ** a large range of pressures, densities and internal energies.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "tillotson.h"

/* We need:
 *
 * tillInit: initialise the library
 *
 * tillInitTable: make the look up table for a) the cold curve and b) the isentropes
 *
 * tillPress: calculate the pressure for a given rho and u from the Tillotson EOS
 *  (also needed for the look up table)
 */

TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, double rhomax)
{
	/*
	 * Initialise a material from the Tillotson library
	 *
	 * We do:
	 * Initialize variables
	 * Convert quantities to code units
	 * Initialize lookup table for the cold curve
	 */

    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */

    TILLMATERIAL *material;
    
    material = malloc(sizeof(TILLMATERIAL));
    assert(material != NULL);

	material->iMaterial = iMaterial;
/*
    material->dKpcUnit = 2.06701e-13;
    material->dMsolUnit = 4.80438e-08;
*/
	/* This two parameters define the unit system we use */
    material->dKpcUnit = dKpcUnit;
    material->dMsolUnit = dMsolUnit;
	material->rhomax = rhomax;

	material->nTableMax = 10000;
	material->nTable = 0;
	material->delta =  material->rhomax/(material->nTableMax-2.0);

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

	/* Lookup table */
    material->cold = malloc(material->nTableMax*sizeof(struct lookup));
    assert(material->cold != NULL);

	/*
	** Set the Tillotson parameters for the material.
	*/

	switch(iMaterial)
	{
		case GRANITE:
			/*
			** Just set granite at the moment.
			*/
			material->a = 0.5;
			material->b = 1.3;
			material->u0 = 1.6e11; /* in ergs/g */
			material->rho0 = 2.7; /* g/cc */
			material->A = 1.8e11; /* ergs/cc */
			material->B = 1.8e11; /* ergs/cc */
			material->us = 3.5e10; /* ergs/g */
			material->us2 = 1.8e11; /* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.79e7; /* ergs/g K (or 790 J/kg K) */ 
			break;
		default:
			/* Unknown material */
			assert(0);
	}

    /*
    ** Convert energies and densities to code units!
    */
    material->u0 /= material->dErgPerGmUnit;
    material->us /= material->dErgPerGmUnit;
    material->us2 /= material->dErgPerGmUnit;
    material->rho0 /= material->dGmPerCcUnit;
    material->A /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
    material->B /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
 
 	material->cv /= material->dErgPerGmUnit;
    
    return(material);
}

void tillFinalizeMaterial(TILLMATERIAL *material)
{
	/* Free the memory */
	free(material);
}

/*
 * This is not implemented yet. Maybe I will not implement it at all.

double tillGamma(TILLMATERIAL *material,double rho,double u) {
    double eta = rho/material->rho0;
    double w0 = u/(material->u0*eta*eta) + 1.0;

    return(material->a + material->b/w0);
    }

*/

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
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		
		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
			*pcSound = sqrt(c2c);
		}
		return (Pc);
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		
		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
			*pcSound = sqrt(c2c);
		}
		return (Pc);
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
		Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));
			*pcSound = sqrt(c2e);
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

		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
			c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));

			*pcSound = sqrt(c2c*(1.0-y)+c2e*y);
		}
	
		return (Pc*(1.0-y)+Pe*y);
	}
}

double tillPressure(TILLMATERIAL *material, double rho, double u)
{
	/* Calculate the pressure from the Tillotson EOS for a material */
	return (tillPressureSound(material, rho, u, NULL));
}

double tillPressureRhoU(TILLMATERIAL material, double rho, double u)
{
	/* Calculate the pressure from the Tillotson EOS for a material */

}

double tillTempRhoU(TILLMATERIAL material, double rho, double u)
{
	/* Calculate T(rho,u) for a material */
}

double tillTempRhoP(TILLMATERIAL *material, double rho, double P)
{
	/* Calculate T(rho,P) for a material */
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
}

double dudrho(TILLMATERIAL *material, double rho, double u)
{
	return(tillPressure(material,rho,u)/(rho*rho));
}

int comparerho(const void* a, const void* b)
/*
** This function compares two entries in the look up table
** and returns -1 if a1.rho < a2.rho, 1 if a1.rho > a2.rho or 0 if
** they are equal (needed to sort the particles with qsort).
*/
{
	struct lookup a1 = *(const struct lookup*)(a);
    struct lookup a2 = *(const struct lookup*)(b);
    
    if (a1.rho < a2.rho) return -1;
    if (a1.rho > a2.rho) return 1;
    return 0;
}

void tillInitColdCurve(TILLMATERIAL *material)
{
	/* Generate the look up table for the cold curve */
    double rho;
    double u;
    double k1u,k2u;
	double h;
    int i;

	rho = material->rho0;
	u = 0.0;
	h = material->delta;
	i = 0;

	material->cold[i].rho = rho;
	material->cold[i].u = u;
	material->cold[i].dudrho = dudrho(material, rho, u);

	++i;

	/*
	** Integrate the condensed and expanded states separately.
	*/
    while (rho <= material->rhomax) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1u = h*dudrho(material,rho,u);
		k2u = h*dudrho(material,rho+0.5*h,u+0.5*k1u);

		u += k2u;
		rho += h;

	    material->cold[i].u = u;
	    material->cold[i].rho = rho;
		material->cold[i].dudrho = dudrho(material, rho, u);
	    ++i;
	}
	
	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	rho = material->rho0;
	u = 0.0;

    while (rho > 0.0) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1u = h*-dudrho(material,rho,u);
		k2u = h*-dudrho(material,rho+0.5*h,u+0.5*k1u);

		u += k2u;
		rho -= h;

	    material->cold[i].u = u;
	    material->cold[i].rho = rho;
		material->cold[i].dudrho = dudrho(material, rho, u);

	    ++i;
	}
	
	--i;	
   	material->nTable = i;

	/* Now sort the look up table. */
    qsort(&material->cold,material->nTable,sizeof(struct lookup),comparerho);
}

