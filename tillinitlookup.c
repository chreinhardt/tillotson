/*
 * Copyright (c) 2014-2018 Christian Reinhardt and Joachim Stadel.
 *
 * This file provides all the functions for the Tillotson EOS library.
 * The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 * and convenient to use equation of state that can describe matter over
 * a large range of pressures, densities and internal energies.
 * 
 * Basic functions:
 *
 * tillInitLookup: generate the look up table for the isentropes (allocate memory!)
 *
 * tillSolveIsentrope: solve the ODE du/drho=P/(rho*rho) for a given initial value v=u(rho0).
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"

int comparerho(const void* a, const void* b)
/*
** This function compares two entries in the look up table
** and returns -1 if a1.rho < a2.rho, 1 if a1.rho > a2.rho or 0 if
** they are equal (needed to sort the particles with qsort).
*/
{
	TILL_LOOKUP_ENTRY a1 = *(const TILL_LOOKUP_ENTRY*)(a);
    TILL_LOOKUP_ENTRY a2 = *(const TILL_LOOKUP_ENTRY*)(b);
    
    if (a1.rho < a2.rho) return -1;
    if (a1.rho > a2.rho) return 1;
    return 0;
}

/// CR: Needs to be adapted for the logrho table
void tillInitColdCurve(TILLMATERIAL *material)
{
    assert(0);
	/* Generate the look up table for the cold curve */
	TILL_LOOKUP_ENTRY *isentrope;

	/* Lookup table */
    material->cold = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));
    assert(material->cold != NULL);

	/* v = u(rho=rho0) */
	isentrope = tillSolveIsentrope(material, 0.0);

	material->cold = isentrope;
	/* Now sort the look up table. */
	// qsort(material->cold,material->nTable,sizeof(TILL_LOOKUP_ENTRY),comparerho);
}

/*
 * Generate the look up table for the isentropic evolution of the internal energy.
 */
void tillInitLookup(TILLMATERIAL *material)
{

	TILL_LOOKUP_ENTRY *isentrope;
	double v, dv;
    int i, j;

	/* We arrange the look up table as a 1D array with Lookup[i][j] = Lookup[i*Ntable+j] */
    material->Lookup = calloc(material->nTableRho*material->nTableV, sizeof(TILL_LOOKUP_ENTRY));
    assert(material->Lookup != NULL);

	v = 0.0;
	// (CR) There was a bug before. We need dv=vmax/(nTableMax-1).
	//dv = material->vmax/(material->nTableMax-1);
	dv = material->dv;

	fprintf(stderr,"tillInitLookup: Solving ODEs.\n");
    
	/*
     * Integrate the isentropes for different v.
     */
	for (j=0; j<material->nTableV; j++)
	{
//        fprintf(stderr,"tillInitLookup: v= %g.\n", v);

		isentrope = tillSolveIsentropeLogRho(material,v);
		
//        fprintf(stderr,"tillInitLookup: Copy data (v= %g).\n", v);

		/* Copy one row to the look up table. This is of course not efficient at all. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices! Lookup[i][j] = Lookup(rho,v). */
			material->Lookup[TILL_INDEX(i,j)] = isentrope[i];		
		}
		
		free(isentrope);
		v += dv;
	}

    /* Solve splines for both u and u1 in v storing the 2nd derivatives wrt v */
	fprintf(stderr,"tillInitLookup: Init splines.\n");
	tillInitSplines(material);
	fprintf(stderr,"tillInitLookup: Splines initialised.\n");	
}

TILL_LOOKUP_ENTRY *tillSolveIsentrope(TILLMATERIAL *material, double v)
{
	/*
	 * Integrate one isentrope for the look up table. The parameter v corresponds to u(rho=rho0)
     * so v=0 gives the cold curve.
     */
    double rho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	TILL_LOOKUP_ENTRY *isentrope;
    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));

	rho = material->rho0;
	u = v;
	h = material->drho;

	i = material->n;

	isentrope[i].rho = rho;
	isentrope[i].v = v;
	isentrope[i].u = u;
	isentrope[i].u1 = tilldudrho(material, rho, u); // du/drho
	
	/*
	** Integrate the condensed and expanded states separately.
	*/
	for (i=material->n+1;i<material->nTableRho;i++)
	{
		double hs = h/10.0;
		
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

			/* Assure that u >= 0.0. */
			if (u < 0.0) u = 0.0;

			rho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].v = v;
		isentrope[i].u1 = tilldudrho(material, rho, u);
	}

	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	rho = material->rho0;
	u = v;

	for (i=material->n-1;i>=0;i--)
	{
		double hs = h/10.0;
		
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

			/* Assure that u >= 0.0. */
			if (u < 0.0) u = 0.0;

			rho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].v = v;
		isentrope[i].u1 = tilldudrho(material, rho, u);
//#if 0
        /*
         * Avoid issues because P/rho^2 is diverging for rho=0.
         */
		if (i == 0)
		{
//			fprintf(stderr,"i=%i,rho(0)=%g,u(0)=%g,u1(0)=%g",i,isentrope[i].rho,isentrope[i].u,isentrope[i].u1);
//			fprintf(stderr,"  tilldudrho=%g P=%g\n",tilldudrho(material, rho, u),tillPressure(material, rho, u));
//			isentrope[i].u1 = 0.0;
			isentrope[i].u1 = isentrope[i+1].u1; // Set u1(0)=u1(drho)
		}
//#endif
	}
	
	return isentrope;
}

/*
 * Integrate one isentrope for the look up table using a log(rho) as a variable. The parameter v
 * corresponds to u(log(rho)=log(rho0)) so v=0 gives the cold curve.
 */
TILL_LOOKUP_ENTRY *tillSolveIsentropeLogRho(TILLMATERIAL *material, double v)
{
    double logrho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	TILL_LOOKUP_ENTRY *isentrope;
    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));

	logrho = log(material->rho0);
	u = v;
	h = material->dlogrho;

	i = material->n;

	isentrope[i].logrho = logrho;
	isentrope[i].v = v;
	isentrope[i].u = u;
	isentrope[i].u1 = tilldudlogrho(material, logrho, u); // du/drho
	
	/*
	** Integrate the condensed and expanded states separately.
	*/
	for (i=material->n+1;i<material->nTableRho;i++)
	{
		double hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*tilldudlogrho(material,logrho,u);
			k2u = hs*tilldudlogrho(material,logrho+0.5*hs,u+0.5*k1u);
			k3u = hs*tilldudlogrho(material,logrho+0.5*hs,u+0.5*k2u);
			k4u = hs*tilldudlogrho(material,logrho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;

			/* Assure that u >= 0.0. */
			if (u < 0.0) u = 0.0;

			logrho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].logrho = logrho;
		isentrope[i].v = v;
		isentrope[i].u1 = tilldudlogrho(material, logrho, u);
	}

	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	logrho = log(material->rho0);
	u = v;

	for (i=material->n-1;i>=0;i--)
	{
		double hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*-tilldudlogrho(material,logrho,u);
			k2u = hs*-tilldudlogrho(material,logrho+0.5*hs,u+0.5*k1u);
			k3u = hs*-tilldudlogrho(material,logrho+0.5*hs,u+0.5*k2u);
			k4u = hs*-tilldudlogrho(material,logrho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;

			/* Assure that u >= 0.0. */
			if (u < 0.0) u = 0.0;

			logrho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].logrho = logrho;
		isentrope[i].v = v;
		isentrope[i].u1 = tilldudlogrho(material, logrho, u);

        /*
         * Avoid issues because P/rho is diverging for rho=0.
         */
		if (i == 0)
		{
//			fprintf(stderr,"i=%i,rho(0)=%g,u(0)=%g,u1(0)=%g",i,isentrope[i].rho,isentrope[i].u,isentrope[i].u1);
//			fprintf(stderr,"  tilldudrho=%g P=%g\n",tilldudrho(material, rho, u),tillPressure(material, rho, u));
//			isentrope[i].u1 = 0.0;
			isentrope[i].u1 = isentrope[i+1].u1; // Set u1(0)=u1(dlogrho)
		}
	}
	
	return isentrope;
}

/*
 * Calculate u2(rho2) by solving the ODE
 *
 * du/dlogrho = P/rho
 *
 * with initial conitions (rho1, u1). Note that the integration is done in log(rho) because thhe
 * ODE is very stiff for small rho otherwise.
 */
double tillCalcU(TILLMATERIAL *material, double rho1, double u1, double rho2)
{
	/* Calculate u2 by solving the ODE */
    double logrho;
    double logrho2;
    double u;
    double k1u, k2u, k3u, k4u;
	double h;

	logrho = log(rho1);
	u = u1;
    logrho2 = log(rho2);

	/* Make smaller steps than we used for look up table. */
	h = material->dlogrho/100.0;
#ifdef TILL_VERBOSE
    fprintf(stderr, "tillCalcU: Solving ODE (rho1= %g u1= %g rho= %g).\n", rho1, u1, rho2);
#endif
	
    if (rho1 < rho2)
	{
		while (logrho < logrho2) {
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = h*tilldudlogrho(material, logrho, u);
			k2u = h*tilldudlogrho(material, logrho+0.5*h, u+0.5*k1u);
			k3u = h*tilldudlogrho(material, logrho+0.5*h, u+0.5*k2u);
			k4u = h*tilldudlogrho(material, logrho+h, u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			logrho += h;

            /* For the last step we set h so that rho == rho2. */
            if (fabs(logrho-logrho2) < h) {
                h = logrho2-logrho;
            }
		}
	} else if (rho1 > rho2) {
		while (logrho > logrho2) {
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = h*tilldudlogrho(material, logrho, u);
			k2u = h*tilldudlogrho(material, logrho+0.5*h, u+0.5*k1u);
			k3u = h*tilldudlogrho(material, logrho+0.5*h, u+0.5*k2u);
			k4u = h*tilldudlogrho(material, logrho+h, u+k3u);

			u -= k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			logrho -= h;

            /* For the last step we set h so that rho == rho2. */
            if (fabs(logrho-logrho2) < h) {
                h = logrho-logrho2;
            }

		}
	}
	return u;
}

/*
 * This function checks if a given (rho, u) is in the look up table or not.
 * 
 * Returns TILL_LOOKUP_SUCCESS if (rho, u) is in the table and an error code if not.
 */
int tillIsInTable(TILLMATERIAL *material, double rho, double u)
{
    double logrho;
    double v_eps = V_EPS;
    int i;

	/* Check if rho <= rhomin or rho >= rhomax */
	if (rho <= material->rhomin)
	{
		return TILL_LOOKUP_OUTSIDE_RHOMIN;
	}

	if (rho >= material->rhomax)
	{
		return TILL_LOOKUP_OUTSIDE_RHOMAX;
	}

    /*
     * Now find rho_i and rho_i+1 so that rho is bracketed (only works for equal spaced grid).
     * 
     * CR: Note that if rho lies on a grid point this can be problematic as, e.g., rho_i can either
     *     be in the interval [i-1, i] or [i, i+1] which is rare cases affects wether or not a
     *     a point close to vmax is in the table. Tests however show, that even if it fails (is in
     *     the table when it should not be) the brent rootfinder still works
     */
	i = tillLookupIndexLogRho(material, log(rho));
	assert(i >= 0 && i < material->nTableRho-1);

    /* Check if v(rho, u) < v_max. This requires interpolation. */
//#if 0
    if ((u < tillSplineIntU(material, material->vmax-v_eps, i)) &&
        (u < tillSplineIntU(material, material->vmax-v_eps, i+1)))
//#endif
#if 0    
    if ((u < material->Lookup[TILL_INDEX(i,material->nTableV-1)].u) &&
        (u < material->Lookup[TILL_INDEX(i+1,material->nTableV-1)].u))
#endif
    {
        /* Check if v(rho, u) > v_0 (so if u > u(rho, 0)). */
        if ((u > material->Lookup[TILL_INDEX(i,0)].u) && (u > material->Lookup[TILL_INDEX(i+1,0)].u))
        {
            /* u(rho, v) is definitely inside of the lookup table. */
#if 0
//        v_eps = 0.5*material->dv;
        fprintf(stderr, "rho= %15.7E u= %15.7E i= %i i+1= %i eps= %15.7E\n", rho, u, i, i+1, v_eps);
        fprintf(stderr, "u_i= %15.7E u_i+1= %15.7E\n", material->Lookup[TILL_INDEX(i,material->nTableV-1)].u, material->Lookup[TILL_INDEX(i+1,material->nTableV-1)].u);
        fprintf(stderr, "u_i= %15.7E u_i+1= %15.7E\n", tillSplineIntU(material, material->vmax-v_eps, i), tillSplineIntU(material, material->vmax-v_eps, i+1));

        fprintf(stderr, "du_i= %15.7E du_i+1= %15.7E\n", material->Lookup[TILL_INDEX(i,material->nTableV-1)].u-tillSplineIntU(material, material->vmax-v_eps, i), material->Lookup[TILL_INDEX(i+1,material->nTableV-1)].u-tillSplineIntU(material, material->vmax-v_eps, i+1));

        fprintf(stderr, "u_i-2= %15.7E u_i-1= %15.7E\n", material->Lookup[TILL_INDEX(i-2,material->nTableV-1)].u, material->Lookup[TILL_INDEX(i-1,material->nTableV-1)].u);
#endif
            return TILL_LOOKUP_SUCCESS;
        } else {
            /* u(rho, v) is below the cold curve. */
#ifdef TILL_VERBOSE
            fprintf(stderr, "tillIsInTable: Value (rho=%g, u=%g) below the cold curve!\n", rho, u);
#endif
            return TILL_LOOKUP_OUTSIDE_VMIN;
        }
    } else {
        /* u(rho, v) is larger than u(rho, v_max) so the lookup table has to be extended. */
#if 0
        fprintf(stderr, "rho= %15.7E u= %15.7E i= %i i+1= %i\n", rho, u, i, i+1);
        fprintf(stderr, "u_i= %15.7E u_i+1= %15.7E\n", material->Lookup[TILL_INDEX(i,material->nTableV-1)].u, material->Lookup[TILL_INDEX(i+1,material->nTableV-1)].u);
        fprintf(stderr, "u_i= %15.7E u_i+1= %15.7E\n", tillSplineIntU(material, material->vmax-v_eps, i), tillSplineIntU(material, material->vmax-v_eps, i+1));
#endif
        return TILL_LOOKUP_OUTSIDE_VMAX;
    }
}

int tillIsBelowColdCurve(TILLMATERIAL *material,double rho,double u)
{
	/*
     * This function checks if a given (rho,u) is in an unphysical state below the cold curve.
     *
     * Returns 1 (true) if (rho,u) is below the cold curve and 0 (false) if not.
	 */
	int iRet = 0;

    /// CR: This needs work!
    fprintf(stderr, "tillIsBelowColdCurve is currently not working!\n");
    exit(1);

	if (u < tillColdULookup(material, rho))
	{
//		printf("tillIsBelowColdCurve: value (%g,%g) below the cold curve (iMat=%i)!\n",rho,u,material->iMaterial);
		iRet = 1;
	}
	return(iRet);
}
