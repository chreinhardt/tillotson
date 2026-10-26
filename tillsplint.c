/*
 * Copyright (c) 2014-2018 Christian Reinhardt and Joachim Stadel.
 *
 * This file provides all the functions for the Tillotson EOS library.
 * The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 * and convenient to use equation of state that can describe matter over
 * a large range of pressures, densities and internal energies.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"

/*
 * Basic functions:
 *
 * tillInitSplines: initialise the cubic splines in u and u1.
 *
 * tillSplineIntU: do a cubic spline interpolation of u in v.
 *
 * tillSplineIntU1: do a cubic spline interpolation of u1 in v.
 *
 * tillCubicInt: interpolate a value (i,i+1) in logrho and (j,j+1) in s (uses Joachim's 2D
 *               interpolator).
 *
 * tillLookupU: u2 on the isentrope through (rho1, u1) at rho2.
 */

#include <stddef.h>   /* offsetof */

/*
 * Natural cubic spline in v along one row of the lookup table, for a field of
 * TILL_LOOKUP_ENTRY. The data is at byte offset offY in each entry, the second
 * derivatives d2y/dv2 are written to byte offset offY2. The grid in v is
 * uniform (dv), so the tridiagonal system is
 *
 *   y2[j-1] + 4 y2[j] + y2[j+1] = 6/dv^2 (y[j+1] - 2 y[j] + y[j-1]),   y2[0] = y2[n-1] = 0,
 *
 * solved by the usual forward sweep and back substitution (Thomas algorithm).
 * nTableV >= 3 is required.
 */
static void tillInitSplineField(TILLMATERIAL *material, size_t offY, size_t offY2)
{
    const int n = material->nTableV;
    const double h = material->dv;
    double *w;      /* work array of the forward sweep */
    int i, j;

#define TILL_FIELD(entry, off) (*(double *)((char *)&(entry) + (off)))

    assert(n >= 3);
    assert(h > 0.0);

    w = malloc(n*sizeof(double));
    assert(w != NULL);

    for (i = 0; i < material->nTableRho; i++) {
        TILL_LOOKUP_ENTRY *row = &material->Lookup[TILL_INDEX(material, i, 0)];

        /* Natural b.c. at v = 0. */
        TILL_FIELD(row[0], offY2) = 0.0;
        w[0] = 0.0;

        /* Forward sweep. */
        for (j = 1; j <= n-2; j++) {
            double p = 0.5*TILL_FIELD(row[j-1], offY2) + 2.0;
            double r = (TILL_FIELD(row[j+1], offY) - 2.0*TILL_FIELD(row[j], offY)
                        + TILL_FIELD(row[j-1], offY))/h;
            TILL_FIELD(row[j], offY2) = -0.5/p;
            w[j] = (3.0*r/h - 0.5*w[j-1])/p;
        }

        /* Natural b.c. at v = vmax, then back substitution. */
        TILL_FIELD(row[n-1], offY2) = 0.0;
        for (j = n-2; j >= 0; j--) {
            TILL_FIELD(row[j], offY2) = TILL_FIELD(row[j], offY2)*TILL_FIELD(row[j+1], offY2) + w[j];
        }
    }

    free(w);
#undef TILL_FIELD
}

/*
 * Calculate the second derivatives of u and u1 (du/dlogrho) with respect to v
 * for every row of the lookup table.
 */
void tillInitSplines(TILLMATERIAL *material)
{
    tillInitSplineU(material);
    tillInitSplineU1(material);
}

/* Spline of u in v: udv2 = d2u/dv2. */
void tillInitSplineU(TILLMATERIAL *material)
{
    tillInitSplineField(material, offsetof(TILL_LOOKUP_ENTRY, u), offsetof(TILL_LOOKUP_ENTRY, udv2));
}

/* Spline of u1 = du/dlogrho in v: u1dv2 = d2u1/dv2. */
void tillInitSplineU1(TILLMATERIAL *material)
{
    tillInitSplineField(material, offsetof(TILL_LOOKUP_ENTRY, u1), offsetof(TILL_LOOKUP_ENTRY, u1dv2));
}

/*
 * Cubic spline interpolation in s of the table field at byte offset offY (with second
 * derivatives at offY2) along row i. s is clamped to the table so that the last interval is
 * closed; callers check the range themselves.
 */
static double tillSplineIntField(TILLMATERIAL *material, double s, int i, size_t offY, size_t offY2)
{
    const double h = material->dv;
    const TILL_LOOKUP_ENTRY *row = &material->Lookup[TILL_INDEX(material, i, 0)];
    int klo = tillLookupIndexV(material, s);
    double a, b;

#define TILL_FIELD(entry, off) (*(const double *)((const char *)&(entry) + (off)))
    a = ((klo+1)*h - s)/h;
    b = 1.0 - a;

    return a*TILL_FIELD(row[klo], offY) + b*TILL_FIELD(row[klo+1], offY)
         + ((a*a*a - a)*TILL_FIELD(row[klo], offY2) + (b*b*b - b)*TILL_FIELD(row[klo+1], offY2))*(h*h)/6.0;
#undef TILL_FIELD
}

/* Cubic spline interpolation of u in s along row i. */
double tillSplineIntU(TILLMATERIAL *material, double s, int i)
{
    return tillSplineIntField(material, s, i, offsetof(TILL_LOOKUP_ENTRY, u), offsetof(TILL_LOOKUP_ENTRY, udv2));
}

/* Cubic spline interpolation of u1 = du/dlogrho in s along row i. */
double tillSplineIntU1(TILLMATERIAL *material, double s, int i)
{
    return tillSplineIntField(material, s, i, offsetof(TILL_LOOKUP_ENTRY, u1), offsetof(TILL_LOOKUP_ENTRY, u1dv2));
}

/*
 * Do a bicubic spline interpolation in logrho and v.
 *
 * The variables are:
 *
 * u[0] = u(v,0)
 * u[1] = u(v,1)
 * dudlogrho[0] = dudlogrho(v,0)
 * dudlogrho[1] = dudlogrho(v,1)
 * dudv[0] = dudv(v,0)
 * dudv[1] = dudv(v,1)
 * dudvdlogrho[0] = dudvdlogrho(v,0)
 * dudvdlogrho[1] = dudvdlogrho(v,1)
 *
 * where v is interpolated between vj and vj+1 using a cubic spline.
 */
void cubicint(double u[2],double dudlogrho[2], double dudv[2], double dudvdlogrho[2],
              double logrho[2], double logrhoint, double *intvalues) {
	double dx, e, e1;
	double ce[4];

	assert(intvalues != NULL);

	dx = logrho[1] - logrho[0];
	e = (logrhoint - logrho[0])/dx;
	e1 = e - 1;

	// these are the 4 Hermite functions
	ce[0] = (2*e + 1)*e1*e1;
	ce[1] = e*e1*e1;
	ce[2] = e*e*(3 - 2*e);
	ce[3] = e*e*e1;

	/*
	**    = ce[0]*u(v,0) + ce[1]*dudlogrho(v,0)*dx + ce[2]*u(v,1) + ce[3]*dudlogrho(v,1);
	** the above is written as 4 independent spline lookups in the table v lies between some j and j+1
	*/
	intvalues[0] = (ce[0]*u[0] + ce[1]*dudlogrho[0]*dx + ce[2]*u[1] + ce[3]*dudlogrho[1]*dx);
	intvalues[1] = (ce[0]*dudv[0] + ce[1]*dudvdlogrho[0]*dx + ce[2]*dudv[1] + ce[3]*dudvdlogrho[1]*dx);
}

/*
 * Use cubicint to interpolate u for a given rho and v.
 *
 * Note: We do a transformation of variables. The actual interpolation is done in logrho.
 */
double tillCubicInt(TILLMATERIAL *material, double rhoint, double vint) {
	double dv, A, B;
	int i, j;
	double *u, *dudlogrho, *dudv, *dudvdlogrho, *logrho, *intvalues;
    double logrhoint;
	double uint;

    /* Transform from rho to logrho. */
    logrhoint = log(rhoint);

	/* Allocate memory */
	u = malloc(2*sizeof(double));
	dudlogrho = malloc(2*sizeof(double));
	dudv = malloc(2*sizeof(double));
	dudvdlogrho = malloc(2*sizeof(double));
	logrho = malloc(2*sizeof(double));

	intvalues = malloc(4*sizeof(double));

    /*
     * Obtain the index i so that log(rho_i) and log(rho_i+1) bracket rho.
     */
	i = tillLookupIndexLogRho(material, logrhoint);

	// vint is between s[j] and s[j+1] (vint is the table variable s)
	j = tillLookupIndexV(material, vint);

	// dv = v[j+1]-v[j]
	// Uniform steps in v
	dv = material->dv;

	logrho[0] = tillLookupLogRho(material, i);
	logrho[1] = tillLookupLogRho(material, i+1);

	/* Do the spline look up in v. */
	u[0] = tillSplineIntU(material, vint, i);
	u[1] = tillSplineIntU(material, vint, i+1);

	dudlogrho[0] = tillSplineIntU1(material, vint, i);
	dudlogrho[1] = tillSplineIntU1(material, vint, i+1);

	/*
	** Calculate dudv from udv2
	** dudv = (u[j+1]-u[j])/dv - (3.0*A*A-1.0)/6.0*dv*udv2[j] + (3.0*B*B-1.0)/6.0*dv*udv2[j+1]
	**
	** Calculate dudvdrho from u1dv2
	** dudvdrho = (u1[j+1]-u1[j])/dv - (3.0*A*A-1.0)/6.0*dv*u1dv2[j] + (3.0*B*B-1.0)/6.0*dv*u1dv2[j+1]
	*/

	/* Calculate dudv from udv2 */
	A = ((j+1)*material->dv-vint)/material->dv;
	B = (vint-j*material->dv)/material->dv;
	
	dudv[0] = (material->Lookup[TILL_INDEX(material, i, j+1)].u-material->Lookup[TILL_INDEX(material, i, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i, j+1)].udv2;
	dudv[1] = (material->Lookup[TILL_INDEX(material, i+1, j+1)].u-material->Lookup[TILL_INDEX(material, i+1, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i+1, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i+1, j+1)].udv2;

	/* Calculate dudvdlogrho from u1dv2 */
	dudvdlogrho[0] = (material->Lookup[TILL_INDEX(material, i, j+1)].u1-material->Lookup[TILL_INDEX(material, i, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i, j+1)].u1dv2;
	dudvdlogrho[1] = (material->Lookup[TILL_INDEX(material, i+1, j+1)].u1-material->Lookup[TILL_INDEX(material, i+1, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i+1, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(material, i+1, j+1)].u1dv2;

	/*
	 * Do the interpolation for u(i,v), u(i+1,v), udlogrho(i,v), udlogrho(i+1,v)
	 */
	cubicint(u, dudlogrho, dudv, dudvdlogrho, logrho, logrhoint, intvalues);
	
	uint = intvalues[0];

	/* Free memory */
	free(u);
	free(dudlogrho);
	free(dudv);
	free(dudvdlogrho);
	free(logrho);
	free(intvalues);

	return uint;
}

double tillFindUonIsentrope(TILLMATERIAL *material, double s, double rho)
{
    return tillCubicInt(material, rho, s);
}

/* f(s) = u(rho, s) - u for the root finder in tillFindEntropyCurve. */
struct tillEntropyCurveCtx { double rho; double u; };

static double denergy(TILLMATERIAL *material, double s, const void *ctx)
{
    const struct tillEntropyCurveCtx *p = (const struct tillEntropyCurveCtx *) ctx;
    return tillFindUonIsentrope(material, s, p->rho) - p->u;
}

/*
 * Find the isentrope through (rho, u): returns its table variable s, or TILL_SENTINEL if (rho, u)
 * is below the cold curve (s < 0) or above the table (s > s_max). iOrder is only used in messages.
 */
double tillFindEntropyCurve(TILLMATERIAL *material, double rho, double u, int iOrder)
{
    struct tillEntropyCurveCtx ctx;
    const double s_max = tillLookupS(material, material->nTableV-1) - V_EPS;
    double fa, fb, s;

    ctx.rho = rho;
    ctx.u = u;

    fa = denergy(material, 0.0, &ctx);
    fb = denergy(material, s_max, &ctx);

    if (fa > 0.0)
        TILL_FAIL_RETURN(TILL_SENTINEL, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i) is below the cold curve.\n",
                         iOrder, rho, u, material->iMaterial);
    if (fb < 0.0)
        TILL_FAIL_RETURN(TILL_SENTINEL, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i) is above the table.\n",
                         iOrder, rho, u, material->iMaterial);

    /* s is O(1-10); 1e-10 in s is 1e-10 relative in v for v >> vmin. Brent converges superlinearly, so this costs little. */
    if (tillBrent(denergy, material, &ctx, 0.0, s_max, fa, fb, 1e-10, 1e-10, 200, &s) != TILL_SUCCESS)
        TILL_FAIL_RETURN(TILL_SENTINEL, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i): root finder failed.\n",
                         iOrder, rho, u, material->iMaterial);

    return s;
}

/*
 * Calculate u2 on the isentrope through (rho1, u1) at rho2. Inside the lookup table this is an
 * interpolation; for rho outside [rhomin, rhomax] or u1 above the table the ODE is integrated
 * directly. Returns TILL_SENTINEL if (rho1, u1) is below the cold curve.
 */
double tillLookupU(TILLMATERIAL *material, double rho1, double u1, double rho2, int iOrder)
{
    double s;
    int iRet;

    iRet = tillIsInTable(material, rho1, u1);

    if (iRet == TILL_LOOKUP_OUTSIDE_VMIN)
        TILL_FAIL_RETURN(TILL_SENTINEL, "tillLookupU: Particle %i (rho= %g, u= %g, iMat= %i) is below the cold curve.\n",
                         iOrder, rho1, u1, material->iMaterial);

    /* Check if the starting and end point are inside of the look up table. */
    if ((iRet == TILL_LOOKUP_SUCCESS) &&
        (rho2 >= material->rhomin*(1.0 - TILL_EDGE_TOL)) && (rho2 <= material->rhomax*(1.0 + TILL_EDGE_TOL)))
    {
        s = tillFindEntropyCurve(material, rho1, u1, iOrder);
        if (s <= TILL_SENTINEL) return TILL_SENTINEL;
        return tillFindUonIsentrope(material, s, rho2);
    }

#ifdef TILL_VERBOSE
    fprintf(stderr, "tillLookupU: value (%g, %g, %g) outside of look up table, doing direct integration (Error: %i).\n",
            rho1, u1, rho2, iRet);
#endif
    return tillCalcU(material, rho1, u1, rho2);
}

/*
 * Calculate u_cold(rho) from the lookup table. Returns TILL_SENTINEL if rho is outside the table.
 *
 * 25.09.2018: Adapted it for the logarithmic lookup table.
 */
double tillColdULookup(TILLMATERIAL *material, double rho)
{
	// Make sure that the look up table is initialized.
	assert(material->Lookup != NULL);

    // In case of a liquid there is not cold term in the expanded states.
    if (rho < material->rho0)
        return 0.0;

    if (!(rho >= material->rhomin*(1.0 - TILL_EDGE_TOL)) || (rho > material->rhomax*(1.0 + TILL_EDGE_TOL)))
        TILL_FAIL_RETURN(TILL_SENTINEL, "tillColdULookup: rho= %g outside of the lookup table (iMat= %i, rhomin= %g, rhomax= %g).\n",
                         rho, material->iMaterial, material->rhomin, material->rhomax);

    return tillCubicInt(material, rho, 0.0);
}

/*
 * Return the index i, so that log(rho_i) and log(rho_i+1) bracket log(rho) in the lookup table.
 */
int tillLookupIndexLogRho(TILLMATERIAL *material, double logrho)
{
    int i;

    // Uniform spacing in log(rho). The result is clamped to [0, nTableRho-2] so that (i, i+1) is
    // always a valid interval; the callers check the range of rho with TILL_EDGE_TOL.
    i = (int) floor((logrho-log(material->rhomin))/material->dlogrho);

    if (i < 0) return 0;
    if (i > material->nTableRho-2) return material->nTableRho-2;

    return i;
}

/*
 * Return the index j, so that s_j and s_j+1 bracket s in the lookup table (s, not v!).
 */
int tillLookupIndexV(TILLMATERIAL *material, double v)
{
    int j;

    // Uniform spacing in s; the argument is s. Clamped to [0, nTableV-2], see tillLookupIndexLogRho.
    j = (int) floor(v/material->dv);

    if (j < 0) return 0;
    if (j > material->nTableV-2) return material->nTableV-2;

    return j;
}

/*
 * Return rho(i).
 */
double tillLookupRho(TILLMATERIAL *material, int i)
{
    // Assume uniform spacing in log(rho).
    return (material->rhomin*exp(i*material->dlogrho));
}

/*
 * Return log(rho(i)).
 */
double tillLookupLogRho(TILLMATERIAL *material, int i)
{
    // Assume uniform spacing in log(rho).
    return (log(material->rhomin)+i*material->dlogrho);
}

/*
 * Return s(j). The grid is uniform in s.
 */
double tillLookupS(TILLMATERIAL *material, int j)
{
    return (j*material->dv);
}

/*
 * Return the physical v(j) = u(rho0) of isentrope j.
 */
double tillLookupV(TILLMATERIAL *material, int j)
{
    return tillVofS(material, tillLookupS(material, j));
}

/*
 * Mapping between the physical isentrope label v and the uniform table variable s.
 */
double tillVofS(TILLMATERIAL *material, double s)
{
    return material->vmin*expm1(s);
}

double tillSofV(TILLMATERIAL *material, double v)
{
    return log1p(v/material->vmin);
}