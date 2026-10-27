/*
 * Copyright (c) 2014-2018 Christian Reinhardt and Joachim Stadel.
 *
 * This file provides all the functions for the Tillotson EOS library.
 * The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 * and convenient to use equation of state that can describe matter over
 * a large range of pressures, densities and internal energies.
 * 
 * tillInitLookup: generate the look up table for the isentropes (allocates memory).
 *
 * tillSolveIsentropeLogRho: integrate du/dlogrho = P/rho for a given initial value v = u(rho0).
 *
 * tillCalcU: the same for an arbitrary start and end point (used outside the table).
 *
 * tillIsInTable: classify a state (rho, u) with respect to the table.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"

/*
 * Generate the look up table for the isentropic evolution of the internal energy.
 */
void tillInitLookup(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmin, double vmax)
{
	TILL_LOOKUP_ENTRY *isentrope;
	double v;
    int i, j;

    /* Check if the material is properly initialized. */
    assert(material != NULL);

    /* Min and max values for the lookup table (in code units). */
    material->rhomin = rhomin;
	material->rhomax = rhomax;
	material->vmin = vmin;
	material->vmax = vmax;
    assert(material->vmin > 0.0 && material->vmin < material->vmax);

    /* rhomin has to be larger than zero otherwise the logarithmic spacing does not work. */
    assert(material->rhomin > 0.0 && material->rhomin < material->rhomax);
        
    /* rhomin has to be smaller than rho0. */
    assert(material->rhomin < material->rho0);

    assert(material->rhomax > 0.0);
    assert(material->vmax > 0.0);

	/* Number of grid points for the look up table. */
	material->nTableRho = nTableRho;
	material->nTableV = nTableV;

    assert(material->nTableRho > 0);
    assert(material->nTableV > 0);

    /*
     * Set dlogrho so that log(rho0) lies on the grid. The index of rho0 is the integer part of
     * x = log(rho0/rhomin)/log(rhomax/rhomin)*nTableRho. For the standard table x is an integer
     * (2/3*1200 = 800) and the logarithms of the code-unit densities make it 800 - 1e-13 in some
     * unit systems: a plain floor then gives 799, a different grid (dlogrho, rhomax) and a
     * different table. The roundoff slack makes the grid the same in every unit system.
     */
    {
        const double x = (log(material->rho0)-log(material->rhomin))/(log(material->rhomax)-log(material->rhomin))*material->nTableRho;
        material->n = (int) floor(x + 1e-9*(1.0 + x));
    }
    assert(material->n >= 1 && material->n <= material->nTableRho-2);
    material->dlogrho = (log(material->rho0)-log(material->rhomin))/material->n;

    /* Set the actual rhomax (it can differ from the requested value after the correction). */
    material->rhomax = tillLookupRho(material, material->nTableRho-1);

    /*
     * The grid is uniform in s = log(1 + v/vmin) (see tillotson.h), so dv is the step in s and
     * the last row is at s = log(1 + vmax/vmin), i.e. v = vmax.
     */
    material->dv = log1p(material->vmax/material->vmin)/(material->nTableV-1);

	/* We arrange the look up table as a 1D array with Lookup[i][j] = Lookup[i*Ntable+j] */
    if (material->Lookup != NULL) free(material->Lookup);
    material->Lookup = calloc((size_t) material->nTableRho*material->nTableV, sizeof(TILL_LOOKUP_ENTRY));
    assert(material->Lookup != NULL);

	/*
     * Integrate the isentropes for different v = v(s_j).
     */
	for (j=0; j<material->nTableV; j++)
	{
		v = tillLookupV(material, j);
		isentrope = tillSolveIsentropeLogRho(material,v);
		
		/* Copy one row to the look up table. This is of course not efficient at all. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices! Lookup[i][j] = Lookup(rho,s). */
			material->Lookup[TILL_INDEX(material,i,j)] = isentrope[i];		
		}
		
		free(isentrope);
	}

    /* Solve splines for both u and u1 in v storing the 2nd derivatives wrt v */
	tillInitSplines(material);
}

/*
 * One classical RK4 step of du/dlogrho = P/rho from (logrho, u) with the signed
 * step h. For h < 0 the stages are evaluated at logrho + h/2 and logrho + h,
 * i.e. below the starting point, which is what a downward integration needs.
 */
static double tillRK4Step(TILLMATERIAL *material, double logrho, double u, double h)
{
    double k1u, k2u, k3u, k4u;

    k1u = h*tilldudlogrho(material, logrho, u);
    k2u = h*tilldudlogrho(material, logrho+0.5*h, u+0.5*k1u);
    k3u = h*tilldudlogrho(material, logrho+0.5*h, u+0.5*k2u);
    k4u = h*tilldudlogrho(material, logrho+h, u+k3u);

    return u + (k1u + 2.0*k2u + 2.0*k3u + k4u)/6.0;
}

/*
 * The standard table (TILL_TABLE_* in tillotson.h): temperatures at rho0 and densities in units
 * of rho0, so the table is the same physical table for every material and unit system.
 */
void tillInitLookupDefault(TILLMATERIAL *material)
{
    assert(material != NULL);
    tillInitLookup(material, TILL_TABLE_NRHO, TILL_TABLE_NV,
                   TILL_TABLE_RHOMIN_FAC*material->rho0, TILL_TABLE_RHOMAX_FAC*material->rho0,
                   material->cv*TILL_TABLE_TMIN, material->cv*TILL_TABLE_TMAX);
}

/*
 * Number of RK4 substeps per table interval dlogrho when integrating an isentrope. The lookup
 * error against a reference integration decreases only about linearly with the number of
 * substeps (the negative-pressure cut is a kink in du/dlogrho, where RK4 loses its order):
 * for granite on the standard table the mean error is 8e-5 / 1e-5 / 2e-6 and the maximum
 * (cold expanded states, T < 100 K at rho0) 3e-3 / 7e-4 / 2e-4 for 1 / 4 / 10 substeps;
 * above 100 K the maximum is 3e-4 / 2e-5 / 1e-5. 10 substeps cost 2.3 s per material.
 */
#ifndef TILL_ISENTROPE_NSUB
#define TILL_ISENTROPE_NSUB 10
#endif

/*
 * Integrate one isentrope for the look up table using log(rho) as the variable. The parameter v
 * corresponds to u(log(rho)=log(rho0)) so v=0 gives the cold curve.
 *
 * The condensed (logrho > log rho0) and expanded (logrho < log rho0) branches are integrated
 * separately, starting from the grid point i = n where rho = rho0, with TILL_ISENTROPE_NSUB RK4
 * substeps per table interval. The step is +dlogrho/NSUB upward and -dlogrho/NSUB downward.
 */
TILL_LOOKUP_ENTRY *tillSolveIsentropeLogRho(TILLMATERIAL *material, double v)
{
    const double logrho0 = log(material->rho0);
    const double hs = material->dlogrho/TILL_ISENTROPE_NSUB;
    TILL_LOOKUP_ENTRY *isentrope;
    double logrho, u;
    int i, s;

    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));
    assert(isentrope != NULL);

    i = material->n;
    isentrope[i].u = v;
    isentrope[i].u1 = tilldudlogrho(material, logrho0, v);

    /* Condensed states: integrate upward from rho0. */
    logrho = logrho0;
    u = v;
    for (i = material->n+1; i < material->nTableRho; i++) {
        for (s = 0; s < TILL_ISENTROPE_NSUB; s++) {
            u = tillRK4Step(material, logrho, u, hs);
            logrho += hs;
            /* Assure that u >= 0.0. */
            if (u < 0.0) u = 0.0;
        }
        /* Use the exact grid value, not the accumulated sum. */
        logrho = logrho0 + (i - material->n)*material->dlogrho;

        isentrope[i].u = u;
        isentrope[i].u1 = tilldudlogrho(material, logrho, u);
    }

    /* Expanded states: integrate downward from rho0 (negative step). */
    logrho = logrho0;
    u = v;
    for (i = material->n-1; i >= 0; i--) {
        for (s = 0; s < TILL_ISENTROPE_NSUB; s++) {
            u = tillRK4Step(material, logrho, u, -hs);
            logrho -= hs;
            /* Assure that u >= 0.0. */
            if (u < 0.0) u = 0.0;
        }
        logrho = logrho0 - (material->n - i)*material->dlogrho;

        isentrope[i].u = u;
        /*
         * Where u is clamped at 0 the tabulated curve is flat, so its slope is 0, not P(rho, 0)/rho
         * (which is positive if the cold term is positive at u = 0, e.g. for B > A). For the
         * built-in materials P(rho < rho0, 0) = 0 and both are the same.
         */
        isentrope[i].u1 = (u == 0.0) ? 0.0 : tilldudlogrho(material, logrho, u);
    }

    /* Avoid issues because P/rho is diverging for rho -> 0: set u1(0) = u1(dlogrho). */
    isentrope[0].u1 = isentrope[1].u1;

    return isentrope;
}

/*
 * Calculate u2(rho2) by solving the ODE
 *
 * du/dlogrho = P/rho
 *
 * with initial conditions (rho1, u1). The integration is done in log(rho) because the ODE is
 * very stiff for small rho otherwise. The interval is split into an integer number of equal
 * RK4 steps no larger than dlogrho/100 (or 1/100 of the interval if no table exists), so the
 * last step lands exactly on log(rho2) in either direction. As in tillSolveIsentropeLogRho, u is
 * clamped at zero after every step, so an expanding state ends on the cold curve (u = 0 below
 * rho0) instead of at a negative u.
 *
 * Returns TILL_FAIL_VALUE (counted in nLookupFail) if rho1 or rho2 is not positive and finite, if
 * u1 is not finite, or if the integration does not give a finite u. A successful integration is
 * not counted here: tillCalcU is also a plain ODE solver (e.g. for reference values), and the
 * caller that uses it as a fallback (tillLookupU) counts the fallback.
 */
double tillCalcU(TILLMATERIAL *material, double rho1, double u1, double rho2)
{
    double logrho1, logrho2, dlogrho;
    double hmax, h, u;
    long nsteps, s;

    /*
     * Guard the input: for a NaN or infinite rho, (long) ceil(|dlogrho|/hmax) below would be
     * undefined behaviour. For finite positive doubles |dlogrho| < 1500, so nsteps is bounded.
     */
    if (!(rho1 > 0.0) || !isfinite(rho1) || !(rho2 > 0.0) || !isfinite(rho2) || !isfinite(u1))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillCalcU: invalid input (rho1= %g u1= %g rho2= %g, iMat= %i).\n",
                         rho1, u1, rho2, material->iMaterial);

    logrho1 = log(rho1);
    logrho2 = log(rho2);
    dlogrho = logrho2 - logrho1;

    if (dlogrho == 0.0) return u1;

    if (material->dlogrho > 0.0) {
        /* Make smaller steps than we used for the look up table. */
        hmax = material->dlogrho/100.0;
    } else {
        hmax = fabs(dlogrho)/100.0;
    }
    assert(hmax > 0.0);

    nsteps = (long) ceil(fabs(dlogrho)/hmax);
    if (nsteps < 1) nsteps = 1;
    h = dlogrho/nsteps;     /* signed */

    TILL_WARN("tillCalcU: Solving ODE (rho1= %g u1= %g rho2= %g, %ld steps).\n", rho1, u1, rho2, nsteps);

    u = u1;
    for (s = 0; s < nsteps; s++) {
        u = tillRK4Step(material, logrho1 + s*h, u, h);
        /* Assure that u >= 0.0 (as in the table integration). */
        if (u < 0.0) u = 0.0;
    }

    if (!isfinite(u))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillCalcU: integration gave u= %g (rho1= %g u1= %g rho2= %g, iMat= %i).\n",
                         u, rho1, u1, rho2, material->iMaterial);

    return u;
}

/*
 * This function checks if a given (rho, u) is in the look up table or not.
 * 
 * Returns TILL_LOOKUP_SUCCESS if (rho, u) is in the table, one of the TILL_LOOKUP_OUTSIDE_*
 * codes if it lies outside a particular bound, and TILL_FAIL if u is not finite.
 */
int tillIsInTable(TILLMATERIAL *material, double rho, double u)
{
    double e, ulo, uhi, utol;
    int i;

    /* The table is closed at both ends with the tolerance TILL_EDGE_TOL (a NaN rho fails too). */
    if (!(rho >= material->rhomin*(1.0 - TILL_EDGE_TOL)))
        return TILL_LOOKUP_OUTSIDE_RHOMIN;
    if (rho > material->rhomax*(1.0 + TILL_EDGE_TOL))
        return TILL_LOOKUP_OUTSIDE_RHOMAX;

    /*
     * A non-finite u is not a state outside one particular bound but no state at all. NaN would
     * otherwise fail the test u <= uhi below and be reported as OUTSIDE_VMAX (and +inf would be
     * reported the same way), which EOSIsUsable tolerates for Tillotson because the pressure is
     * analytic. TILL_FAIL is mapped to EOS_FAIL by EOSIsInTable and rejected by EOSIsUsable.
     */
    if (!isfinite(u))
    {
        TILL_WARN("tillIsInTable: Value (rho=%15.7E, u=%15.7E, iMat= %i) has a non-finite u!\n", rho, u, material->iMaterial);
        return TILL_FAIL;
    }

    /* u on the cold curve (j = 0) and on the top row of the table at this rho. */
    tillLocateRho(material, rho, &i, &e);
    ulo = tillCubicIntGrid(material, i, e, 0);
    uhi = tillCubicIntGrid(material, i, e, material->nTableV-1);
    utol = tillUTol(material, u);

    if (u > uhi + utol)
        return TILL_LOOKUP_OUTSIDE_VMAX;
    if (u < ulo - utol)
    {
        TILL_WARN("tillIsInTable: Value (rho=%15.7E, u=%15.7E, iMat= %i) below the cold curve!\n", rho, u, material->iMaterial);
        return TILL_LOOKUP_OUTSIDE_VMIN;
    }

    return TILL_LOOKUP_SUCCESS;
}

/*
 * This function checks if a given (rho,u) is in an unphysical state below the cold curve.
 *
 * Returns 1 (true) if (rho,u) is below the cold curve, 0 (false) if not, and TILL_FAIL if
 * rho is outside the table or u is not finite.
 */
int tillIsBelowColdCurve(TILLMATERIAL *material, double rho, double u)
{
    int iRet = 1;

    if (!(rho >= material->rhomin*(1.0 - TILL_EDGE_TOL)) || (rho > material->rhomax*(1.0 + TILL_EDGE_TOL)))
        TILL_FAIL_RETURN(material, TILL_FAIL, "tillIsBelowColdCurve: rho= %g outside of the lookup table (iMat= %i, rhomin= %g, rhomax= %g).\n",
                         rho, material->iMaterial, material->rhomin, material->rhomax);

    if (!isfinite(u))
        TILL_FAIL_RETURN(material, TILL_FAIL, "tillIsBelowColdCurve: non-finite u= %g (rho= %g, iMat= %i).\n",
                         u, rho, material->iMaterial);

    if (u >= tillColdULookup(material, rho)) {
        iRet = 0;
    }

    return iRet;
}