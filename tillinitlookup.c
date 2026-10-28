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
 * tillInitColdCurve: only the grid and the cold curve and top isentrope (all that the
 * temperature and the table status need).
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
 * The grid of the table (tillInitLookup and tillInitColdCurve share it): rho and s axes and the
 * derived constants. Does not allocate.
 */
static void tillInitGrid(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmin, double vmax)
{
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

    /* Derived grid constants for the lookups (tillLookupIndexLogRho, tillLocateRho, the splines). */
    material->logrhomin = log(material->rhomin);
    material->invdlogrho = 1.0/material->dlogrho;
    material->invdv = 1.0/material->dv;
}

/*
 * The two edge isentropes, row 0 (v = 0, the cold curve) and row nTableV-1 (v = vmax), on the grid
 * of tillInitGrid. The same tillSolveIsentropeLogRho calls as for the table rows, with v from
 * tillLookupV, so the rows are bit-identical.
 */
static void tillInitEdgeCurves(TILLMATERIAL *material)
{
    if (material->ColdCurve != NULL) free(material->ColdCurve);
    if (material->TopCurve != NULL) free(material->TopCurve);
    material->ColdCurve = tillSolveIsentropeLogRho(material, tillLookupV(material, 0));
    material->TopCurve = tillSolveIsentropeLogRho(material, tillLookupV(material, material->nTableV-1));
}

/*
 * The grid and the edge curves only (see tillotson.h). A full table of a previous call is freed:
 * it may belong to another grid.
 */
void tillInitColdCurve(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmin, double vmax)
{
    assert(material != NULL);
    tillInitGrid(material, nTableRho, nTableV, rhomin, rhomax, vmin, vmax);
    if (material->Lookup != NULL) {
        free(material->Lookup);
        material->Lookup = NULL;
    }
    tillInitEdgeCurves(material);
}

void tillInitColdCurveDefault(TILLMATERIAL *material)
{
    assert(material != NULL);
    tillInitColdCurve(material, TILL_TABLE_NRHO, TILL_TABLE_NV,
                      TILL_TABLE_RHOMIN_FAC*material->rho0, TILL_TABLE_RHOMAX_FAC*material->rho0,
                      material->cv*TILL_TABLE_TMIN, material->cv*TILL_TABLE_TMAX);
}

/*
 * Generate the look up table for the isentropic evolution of the internal energy (and the edge
 * curves, which are its first and last row).
 */
void tillInitLookup(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmin, double vmax)
{
	TILL_LOOKUP_ENTRY *isentrope;
	double v;
    int i, j;

    assert(material != NULL);
    tillInitGrid(material, nTableRho, nTableV, rhomin, rhomax, vmin, vmax);
    tillInitEdgeCurves(material);

	/* We arrange the look up table as a 1D array with Lookup[i][j] = Lookup[i*Ntable+j] */
    if (material->Lookup != NULL) free(material->Lookup);
    material->Lookup = calloc((size_t) material->nTableRho*material->nTableV, sizeof(TILL_LOOKUP_ENTRY));
    assert(material->Lookup != NULL);

	/*
     * Integrate the isentropes for different v = v(s_j); the edge rows are the edge curves.
     */
	for (j=0; j<material->nTableV; j++)
	{
        if (j == 0) {
            isentrope = material->ColdCurve;
        } else if (j == material->nTableV-1) {
            isentrope = material->TopCurve;
        } else {
		    v = tillLookupV(material, j);
		    isentrope = tillSolveIsentropeLogRho(material,v);
        }
		
		/* Copy one row to the look up table. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices! Lookup[i][j] = Lookup(rho,s). */
			material->Lookup[TILL_INDEX(material,i,j)] = isentrope[i];		
		}
		
		if (j != 0 && j != material->nTableV-1) free(isentrope);
	}

    /* Solve splines for both u and u1 in v storing the 2nd derivatives wrt v */
	tillInitSplines(material);
}

/*
 * Fatal if the full isentrope table is not built (see tillotson.h): a programming error, the
 * caller must build it with tillInitLookup (in EOSlib: EOSinitIsentropicLookup) first.
 */
void tillRequireLookup(const TILLMATERIAL *material, const char *chFunc)
{
    if (material->Lookup == NULL) {
        fprintf(stderr, "%s: the isentrope table of Tillotson material %i is not built. Call tillInitLookup "
                "(EOSinitIsentropicLookup in EOSlib) before the isentropic lookups; tillInitColdCurve builds "
                "only the cold curve and the top isentrope.\n", chFunc, material->iMaterial);
        assert(0);
        exit(1);
    }
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
 * Local error tolerance of the adaptive isentrope integration (tillIntegrateIsentrope), relative
 * to |u| with the floor vmin (the energy scale of the table, as in tillUTol). It is used for the
 * table (tillSolveIsentropeLogRho) and for the direct integration outside it (tillCalcU). The
 * global error is about the number of steps times this.
 */
#ifndef TILL_ODE_RTOL
#define TILL_ODE_RTOL 1e-11
#endif

/* Upper bound on the number of attempted steps of one integration (a guard, not reached in practice). */
#define TILL_ODE_MAXSTEPS 10000000L

/*
 * Integrate du/dlogrho = P/rho from (x1, u1) to x2 (x = log(rho), either direction) with adaptive
 * RK4: each step of size h is compared with two steps of h/2, the difference estimates the local
 * error, the step is accepted if that is below TILL_ODE_RTOL*(|u| + vmin), and the next step is
 * scaled with the usual (tol/err)^(1/5), at most by 4. The steps are small only near the kinks of
 * P/rho (the negative-pressure cut, the region boundaries us and us2, the Melosh cut at 0.8 rho0),
 * where RK4 loses its order, and large where P/rho is smooth. The last step lands exactly on x2.
 * u is clamped at zero after every (half) step, so an expanding state ends on the cold curve
 * (u = 0 below rho0) instead of at a negative u.
 *
 * *ph is the trial step on input (its sign is ignored) and the last accepted step size before the
 * clamp to x2 on output, so a caller that integrates interval by interval keeps the step size.
 * Returns 0 on success, -1 if the integration gave a non-finite error estimate or did not reach
 * x2 within TILL_ODE_MAXSTEPS attempts.
 */
static int tillIntegrateIsentrope(TILLMATERIAL *material, double x1, double u1, double x2, double *ph, double *pu)
{
    const double dir = (x2 > x1) ? 1.0 : -1.0;
    double x = x1, u = u1, h = dir*fabs(*ph), hnext = h;
    long nattempt;

    if (x2 == x1) { *pu = u1; return 0; }

    for (nattempt = 0; nattempt < TILL_ODE_MAXSTEPS; nattempt++) {
        double uf, uh, err, tol, fac;
        int bLast = 0;

        /* Do not step past x2; the last step lands on it exactly. */
        if (fabs(h) >= fabs(x2 - x)) {
            h = x2 - x;
            bLast = 1;
        }

        uf = tillRK4Step(material, x, u, h);
        if (uf < 0.0) uf = 0.0;
        uh = tillRK4Step(material, x, u, 0.5*h);
        if (uh < 0.0) uh = 0.0;
        uh = tillRK4Step(material, x + 0.5*h, uh, 0.5*h);
        if (uh < 0.0) uh = 0.0;

        err = fabs(uh - uf)/15.0;
        tol = TILL_ODE_RTOL*(fabs(uh) + material->vmin);
        if (!isfinite(err) || !isfinite(uh)) return -1;

        if (err <= tol || fabs(h) <= 1e-12*(1.0 + fabs(x))) {
            /* Accepted (or refined to the roundoff limit of x). */
            u = uh;
            fac = (err > 0.0) ? fmin(0.9*pow(tol/err, 0.2), 4.0) : 4.0;
            if (bLast) {
                /* A clamped last step says little about the step size: keep the previous one if larger. */
                *ph = fmax(fabs(hnext), fabs(h)*fac);
                *pu = u;
                return 0;
            }
            x += h;
            h *= fac;
            hnext = h;
        } else {
            h *= fmax(0.9*pow(tol/err, 0.2), 0.1);
        }
    }
    return -1;
}

/*
 * Integrate one isentrope for the look up table using log(rho) as the variable. The parameter v
 * corresponds to u(log(rho)=log(rho0)) so v=0 gives the cold curve.
 *
 * The condensed (logrho > log rho0) and expanded (logrho < log rho0) branches are integrated
 * separately, starting from the grid point i = n where rho = rho0, interval by interval with the
 * adaptive integrator (tillIntegrateIsentrope), which ends every interval exactly on the grid
 * point and carries its step size over to the next interval. Up to tillotson 4.0.0 every interval
 * took TILL_ISENTROPE_NSUB = 10 fixed RK4 steps: the node values were good to ~1e-4 relative
 * (cold states in tension) and the build took 2.3 s per material.
 */
TILL_LOOKUP_ENTRY *tillSolveIsentropeLogRho(TILLMATERIAL *material, double v)
{
    const double logrho0 = log(material->rho0);
    TILL_LOOKUP_ENTRY *isentrope;
    double logrho, u, h;
    int i;

    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));
    assert(isentrope != NULL);

    i = material->n;
    isentrope[i].u = v;
    isentrope[i].u1 = tilldudlogrho(material, logrho0, v);

    /* Condensed states: integrate upward from rho0. */
    u = v;
    h = material->dlogrho;
    for (i = material->n+1; i < material->nTableRho; i++) {
        /* The exact grid values, not an accumulated sum. */
        const double xa = logrho0 + (i - 1 - material->n)*material->dlogrho;
        logrho = logrho0 + (i - material->n)*material->dlogrho;
        if (tillIntegrateIsentrope(material, xa, u, logrho, &h, &u) != 0) {
            fprintf(stderr, "tillSolveIsentropeLogRho: integration failed (iMat= %i, v= %g, i= %d).\n", material->iMaterial, v, i);
            assert(0);
            exit(1);
        }

        isentrope[i].u = u;
        isentrope[i].u1 = tilldudlogrho(material, logrho, u);
    }

    /* Expanded states: integrate downward from rho0. */
    u = v;
    h = material->dlogrho;
    for (i = material->n-1; i >= 0; i--) {
        const double xa = logrho0 - (material->n - i - 1)*material->dlogrho;
        logrho = logrho0 - (material->n - i)*material->dlogrho;
        if (tillIntegrateIsentrope(material, xa, u, logrho, &h, &u) != 0) {
            fprintf(stderr, "tillSolveIsentropeLogRho: integration failed (iMat= %i, v= %g, i= %d).\n", material->iMaterial, v, i);
            assert(0);
            exit(1);
        }

        isentrope[i].u = u;
        /*
         * Where u is clamped at 0 the tabulated curve is flat, so its slope is 0, not P(rho, 0)/rho
         * (which is positive if the cold term is positive at u = 0, e.g. for B > A). For the
         * built-in materials P(rho < rho0, 0) = 0 and both are the same.
         */
        isentrope[i].u1 = (u == 0.0) ? 0.0 : tilldudlogrho(material, logrho, u);
    }

    /*
     * Node 0 (rhomin) keeps its own slope like every other node, so the Hermite interval [0, 1] is
     * consistent with the data. P/rho does not diverge there: in the expanded hot region the
     * A mu exp(-(alpha z + beta z^2)) term underflows (z ~ 1e4 at 1e-4 rho0) and P/rho -> Gamma u,
     * and the cold term A mu + B mu^2 ~ B - A is negative and cut for B <= A (all built-in
     * materials). Only for B > A does P/rho grow like (B - A)/rho in the cold and intermediate
     * regions; the true slope is then still the right Hermite data. (Up to tillotson 4.0.0 u1(0)
     * was replaced by u1(1), which made the first interval 10-1000x less accurate.)
     */

    return isentrope;
}

/*
 * Calculate u2(rho2) by solving the ODE
 *
 * du/dlogrho = P/rho
 *
 * with initial conditions (rho1, u1), with the same adaptive integrator as the table
 * (tillIntegrateIsentrope; up to tillotson 4.0.0 fixed steps of dlogrho/100, ~1e5 steps for an
 * expansion from rho0 to 1e-7 rho0). The integration is done in log(rho) because the ODE is very
 * stiff for small rho otherwise. The first trial step is the table step if a table exists, else
 * 1/100 of the interval.
 *
 * Returns TILL_FAIL_VALUE (counted in nLookupFail) if rho1 or rho2 is not positive and finite, if
 * u1 is not finite, or if the integration does not give a finite u. A successful integration is
 * not counted here: tillCalcU is also a plain ODE solver (e.g. for reference values), and the
 * caller that uses it as a fallback (tillLookupU) counts the fallback.
 */
double tillCalcU(TILLMATERIAL *material, double rho1, double u1, double rho2)
{
    double logrho1, logrho2, h, u = 0.0;

    if (!(rho1 > 0.0) || !isfinite(rho1) || !(rho2 > 0.0) || !isfinite(rho2) || !isfinite(u1))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillCalcU: invalid input (rho1= %g u1= %g rho2= %g, iMat= %i).\n",
                         rho1, u1, rho2, material->iMaterial);

    logrho1 = log(rho1);
    logrho2 = log(rho2);

    if (logrho2 == logrho1) return u1;

    h = (material->dlogrho > 0.0) ? material->dlogrho : fabs(logrho2 - logrho1)/100.0;

    TILL_WARN("tillCalcU: Solving ODE (rho1= %g u1= %g rho2= %g).\n", rho1, u1, rho2);

    if (tillIntegrateIsentrope(material, logrho1, u1, logrho2, &h, &u) != 0 || !isfinite(u))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillCalcU: integration failed (rho1= %g u1= %g rho2= %g, iMat= %i).\n",
                         rho1, u1, rho2, material->iMaterial);

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
 * Returns 1 (true) if (rho,u) is below the cold curve, 0 (false) if not, and TILL_FAIL (counted
 * once) where the cold curve is not defined, i.e. exactly where tillColdULookup fails (rho > rhomax,
 * rho <= 0, NaN), or if u is not finite. Below rhomin the cold curve is zero (see tillColdULookup),
 * so u < 0 is below it there. Up to 4.0.0 this function failed for every rho outside [rhomin,
 * rhomax], while tillColdULookup, and with it EOSisbelowColdCurve, gave a valid answer below rhomin.
 */
int tillIsBelowColdCurve(TILLMATERIAL *material, double rho, double u)
{
    double ucold;

    if (!isfinite(u))
        TILL_FAIL_RETURN(material, TILL_FAIL, "tillIsBelowColdCurve: non-finite u= %g (rho= %g, iMat= %i).\n",
                         u, rho, material->iMaterial);

    /* Counted there if it fails (rho > rhomax, rho <= 0, NaN). */
    ucold = tillColdULookup(material, rho);
    if (tillFailed(ucold)) return TILL_FAIL;

    return (u < ucold) ? 1 : 0;
}