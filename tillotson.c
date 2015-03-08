/*
 ** This file provides all the functions for the Tillotson EOS library.
 ** The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 ** and convenient to use equation of state that can describe matter over
 ** a large range of pressures, densities and internal energies.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"

#include "interpol/coeff.h"
#include "interpol/interpol.h"

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
	 * Allocate memory for the lookup tables
	 * Initialize lookup table and the cold curve
	 */

    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */

    TILLMATERIAL *material;
	int i;
	 
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
	
	/* Just as a first step we have equal steps in rho and v. */
	material->vmax = material->rhomax;

	material->nTableMax = 10000;
	/* Needs about 800M memory. */
	material->nTableMax = 10000;

	/* For debugging purpose. */
#ifdef TILL_USE_RK4
	material->nTableMax = 100;
#else
	material->nTableMax = 1000;
#endif
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
	
	/* We arrange the look up table as a 1D array with Lookup[i][j] = Lookup[i*Ntable+j] */
    material->Lookup = malloc(material->nTableMax*material->nTableMax*sizeof(double));
    assert(material->Lookup != NULL);

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
		case IRON:
			/*
			** Material parameters from Benz 1987.
			*/
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 9.5e10; /* in ergs/g */
			material->rho0 = 7.86; /* g/cc */
			material->A = 1.28e12; /* ergs/cc */
			material->B = 1.05e12; /* ergs/cc */
			material->us = 1.425e10; /* ergs/g */
			material->us2 = 8.45e10; /* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.449e7; /* ergs/g K */ 
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

	/* Set delta so that rho0 lies on the grid. */
	material->n = floor(material->rho0/material->rhomax*material->nTableMax);
 	material->delta =  material->rho0/material->n;
	
	/* Set the actual rhomax. */ 
  	material->rhomax = material->delta*(material->nTableMax-1);
  
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
	**  Here we evaluate, which part of the equation of state we need.
	*/
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		dPcdrho = (material->a+material->b/w0*(3.0-2.0/w0))*u + (material->A+2.0*material->B*mu)/material->rho0;
		
		return (dPcdrho);
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		dPcdrho = (material->a+material->b/w0*(3.0-2.0/w0))*u + (material->A+2.0*material->B*mu)/material->rho0;

		return (dPcdrho);
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		dPedrho = (material->a + material->b/w0*exp(-material->beta*z*z)*(2.0*material->beta*z/eta+3.0-2.0/w0))*u+material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta));		
		return (dPedrho);
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);

		dPcdrho = (material->a+material->b/(w0*w0)*(3.0-2.0/w0))*u + (material->A+2.0*material->B*mu)/material->rho0;
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
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		dPcdu = (material->a+material->b/(w0*w0))*rho;
		
		return (dPcdu);
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		dPcdu = (material->a+material->b/(w0*w0))*rho;
		return (dPcdu);
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		dPedu = (material->a + material->b/(w0*w0)*exp(-material->beta*z*z))*rho;		
		return (dPedu);
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);

		dPcdu = (material->a+material->b/(w0*w0))*rho;
		dPedu = (material->a + material->b/(w0*w0)*exp(-material->beta*z*z))*rho;		

		return (dPcdu*(1.0-y)+dPedu*y);
	}
}

double tilldTdrho(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dT/drho at u=const.
	*/
	return (-1.0/material->cv*tillPressure(material,rho,tillColdULookup(material,rho))*(rho*rho));
}

double tilldTdu(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dT/du at rho=const.
	*/
	return (-1.0/material->cv);
}

double tillPressureRhoU(TILLMATERIAL material, double rho, double u)
{
	/* Calculate the pressure from the Tillotson EOS for a material */

}

double tillTempRhoU(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate T(rho,u) for a material. As an approximation
	** we use u(rho,T) = uc(rho) + cv*T.
	*/
	return ((u-tillColdULookup(material,rho))/material->cv);
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

double tilldudrho(TILLMATERIAL *material, double rho, double u)
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
	struct lookup *isentrope;

	/* v = u(rho=rho0) */
	isentrope = tillSolveIsentrope(material,0.0);

	material->cold = isentrope;
	/* Now sort the look up table. */
	// qsort(material->cold,material->nTable,sizeof(struct lookup),comparerho);
}

void tillInitLookup(TILLMATERIAL *material)
{
	/*
	** Generate the look up table for the isentropic evolution of
	** the internal energy.
	*/

	struct lookup *isentrope;
	double v;
    int i,j;

	v = 0.0;
	//fprintf(stderr, "Starting integration...\n");

	/*
	** Integrate the isentropes for different v.
	*/
	for (i=0; i<material->nTableMax; i++)
	{
		isentrope = tillSolveIsentrope(material,v);
		
		/* Copy one row to the look up table. This is of course not efficitent at all. */
		for (j=0; j<material->nTableMax; j++)
		{
			/* Careful with the indices! */
			material->Lookup[material->nTableMax*j+i] = isentrope[j].u;		
		}
		
		free(isentrope);
		v += material->delta;
	}

	/* Initialize the coefficients for the interpolation. */
	SamplesToCoefficients(material->Lookup,material->nTableMax,material->nTableMax, TILL_SPLINE_DEGREE);
}
#ifdef TILL_USE_RK4
struct lookup *tillSolveIsentrope(TILLMATERIAL *material, double v)
{
	/*
	** Integrate one isentrope for the look up table. The parameter v corresponds to
	** u(rho=rho0) so v=0 gives the cold curve.
	*/
    double rho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	struct lookup *isentrope;
    isentrope = malloc(material->nTableMax*sizeof(struct lookup));

	rho = material->rho0;
	u = v;
	h = material->delta;

	i = material->n;

	isentrope[i].rho = rho;
	isentrope[i].u = u;
//	isentrope[i].dudrho = tilldudrho(material, rho, u);

	/*
	** Integrate the condensed and expanded states separately.
	*/
	for (i=material->n+1;i<material->nTableMax;i++)
	{
		float hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*tilldudrho(material,rho,u);
			k2u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k1u);
			k3u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k2u);
			k4u = hs*tilldudrho(material,rho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			rho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].rho = rho;
//		isentrope[i].dudrho = tilldudrho(material, rho, u);
	}
	
	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	rho = material->rho0;
	u = v;

	for (i=material->n-1;i>=0;i--)
	{
		float hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*-tilldudrho(material,rho,u);
			k2u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k1u);
			k3u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k2u);
			k4u = hs*-tilldudrho(material,rho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			rho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].rho = rho;
//		isentrope[i].dudrho = tilldudrho(material, rho, u);
	}

	return isentrope;
}
#else
struct lookup *tillSolveIsentrope(TILLMATERIAL *material, double v)
{
	/*
	** Integrate one isentrope for the look up table. The parameter v corresponds to
	** u(rho=rho0) so v=0 gives the cold curve.
	*/
    double rho;
    double u;
    double k1u,k2u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	struct lookup *isentrope;
    isentrope = malloc(material->nTableMax*sizeof(struct lookup));

	rho = material->rho0;
	u = v;
	h = material->delta;

	i = material->n;

	isentrope[i].rho = rho;
	isentrope[i].u = u;
//	isentrope[i].dudrho = tilldudrho(material, rho, u);

	/*
	** Integrate the condensed and expanded states separately.
	*/
	for (i=material->n+1;i<material->nTableMax;i++)
	{
		float hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
		{
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = hs*tilldudrho(material,rho,u);
			k2u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k1u);

			u += k2u;
			rho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].rho = rho;
//		isentrope[i].dudrho = tilldudrho(material, rho, u);
	}
	
	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	rho = material->rho0;
	u = v;

	for (i=material->n-1;i>=0;i--)
	{
		float hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
		{
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = hs*-tilldudrho(material,rho,u);
			k2u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k1u);

			u += k2u;
			rho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].rho = rho;
//		isentrope[i].dudrho = tilldudrho(material, rho, u);
	}
	
	return isentrope;
}
#endif

/* Prasenjits root finder */
float brent(float (*func)(TILLMATERIAL *,float,float,float),TILLMATERIAL *material,float a,float b,float rho,float u,float tol,int iOrder);

float tillFindUonIsentrope(TILLMATERIAL *material,float v,float rho)
{
	float iv,irho,u;
	/* Needed for the interpolation function. */
	iv = (material->nTableMax-1)*v/material->vmax;
	irho = (material->nTableMax-1)*rho/material->rhomax;

	return InterpolatedValue(material->Lookup,material->nTableMax,material->nTableMax,iv,irho,TILL_SPLINE_DEGREE);
}

float denergy(TILLMATERIAL *material,float v,float rho,float u)
{
	return (tillFindUonIsentrope(material,v,rho)-u);
}

/* Find isentrope for a given rho and u */
float tillFindEntropyCurve(TILLMATERIAL *material,float rho,float u,int iOrder)
{
	float tol=1e-6;
	return brent(denergy,material,0,material->vmax,rho,u,tol,iOrder);
}

double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder)
{
	/* Calculates u2 for a given rho1,u2,rho2. */
	double v;

	v = tillFindEntropyCurve(material,rho1,u1,iOrder);

	return tillFindUonIsentrope(material,v,rho2);
}
	 
double tillColdULookup(TILLMATERIAL *material,double rho)
{
    double x,xi;
	double drho;
    int i;

    i = material->nTableMax-1;
	
	/* What do we do if rho > rhomax */
	/* if (r >= material->r[i]) return(material->rho[i]*exp(-(r-material->r[i]))); */
	
	x = rho/material->delta;
	xi = floor(x);
	assert(xi >= 0.0);
	x -= xi;

	i = (int)xi;
	if (i < 0)
	{
		fprintf(stderr,"ERROR rho:%.14g x:%.14g xi:%.14g i:%d\n",rho,x,xi,i);
	}
    assert(i >= 0);

	if (i > material->nTableMax-2) fprintf(stderr,"ERROR: out of bounds rho:%.14g rhomax:%.14g i:%i nTableMax: %i\n",rho,material->rhomax,i,material->nTableMax);
	assert(i < material->nTableMax-1);

	if (i <= material->nTableMax-2)
	{
		/* linear interpolation for now. */
		return(material->cold[i].u*(1.0-x) + material->cold[i+1].u*x);
	}
	
	/* This would only be needed if we cut off the model between the last two steps.	
	if (i == material->nTable-2)
	{
		dr = material->r[i+1] - material->r[i];
		x = r/dr;
		xi = floor(x);
		x -= xi;
		return(material->rho[i]*(1.0-x) + material->rho[i+1]*x);
	*/
	/* What do we do if i >= nTable-1
	} else {
		i = material->nTable - 1;
		return(material->rho[i]*exp(-(r-material->r[i])));
	}*/
}
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2)
{
	/* Calculate u2 by solving the ODE */
    double rho;
    double u;
    double k1u,k2u;
	double h;

	rho = rho1;
	u = u1;
	/* Make smaller steps than we used for look up table. */
	h = material->delta/100.0;

	if (rho1 < rho2)
	{
		while (rho < rho2) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = h*tilldudrho(material,rho,u);
			k2u = h*tilldudrho(material,rho+0.5*h,u+0.5*k1u);
	
			u += k2u;
			rho += h;
		}
	} else if (rho1 > rho2) {
		while (rho > rho2) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = h*-tilldudrho(material,rho,u);
			k2u = h*-tilldudrho(material,rho+0.5*h,u+0.5*k1u);

			u += k2u;
			rho -= h;
		}
	}
	return u;
}


