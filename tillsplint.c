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
 * Cubic Hermite interpolation between (y0, d0) at e = 0 and (y1, d1) at e = 1, where the
 * derivatives d are with respect to the full interval dx. Returns the value; if pdy is not NULL
 * the derivative with respect to e/dx is stored there.
 */
static inline double tillHermite(double y0, double y1, double d0, double d1, double dx, double e, double *pdy)
{
    const double e1 = e - 1.0;
    const double h00 = (2.0*e + 1.0)*e1*e1;     /* 1 - 3e^2 + 2e^3 */
    const double h10 = e*e1*e1;                 /* e - 2e^2 + e^3 */
    const double h01 = e*e*(3.0 - 2.0*e);
    const double h11 = e*e*e1;

    if (pdy != NULL) {
        *pdy = (6.0*e*e1*(y0 - y1)/dx) + (3.0*e*e - 4.0*e + 1.0)*d0 + (3.0*e*e - 2.0*e)*d1;
    }
    return h00*y0 + h10*d0*dx + h01*y1 + h11*d1*dx;
}

/* Position of logrho in the rho grid: interval i and fraction e in [0, 1] (clamped at the edges). */
void tillLocateRho(TILLMATERIAL *material, double rho, int *pi, double *pe)
{
    const double logrho = log(rho);
    int i = tillLookupIndexLogRho(material, logrho);
    *pi = i;
    *pe = (logrho - tillLookupLogRho(material, i))/material->dlogrho;
}

/*
 * u(rho, s_j) on a grid line of the table: Hermite interpolation in logrho between the rows i
 * and i+1 with the tabulated u1 = du/dlogrho. No spline in s is needed on a grid line. Used for
 * the cold curve (j = 0), the upper edge of the table and for bracketing in s.
 */
double tillCubicIntGrid(TILLMATERIAL *material, int i, double e, int j)
{
    const TILL_LOOKUP_ENTRY *p0 = &material->Lookup[TILL_INDEX(material, i, j)];
    const TILL_LOOKUP_ENTRY *p1 = &material->Lookup[TILL_INDEX(material, i+1, j)];

    return tillHermite(p0->u, p1->u, p0->u1, p1->u1, material->dlogrho, e, NULL);
}

/*
 * Cubic spline interpolation in s along row i of the field at offY (second derivatives at offY2),
 * with the derivative d/ds if pdy is not NULL. s must be inside the table (the interval index is
 * clamped).
 */
static inline double tillSplineIntFieldD(TILLMATERIAL *material, double s, int i, size_t offY, size_t offY2, double *pdy)
{
    const double h = material->dv;
    const TILL_LOOKUP_ENTRY *row = &material->Lookup[TILL_INDEX(material, i, 0)];
    const int klo = tillLookupIndexV(material, s);
    const double a = ((klo+1)*h - s)/h, b = 1.0 - a;
#define TILL_FIELD(entry, off) (*(const double *)((const char *)&(entry) + (off)))
    const double ylo = TILL_FIELD(row[klo], offY), yhi = TILL_FIELD(row[klo+1], offY);
    const double y2lo = TILL_FIELD(row[klo], offY2), y2hi = TILL_FIELD(row[klo+1], offY2);
#undef TILL_FIELD

    if (pdy != NULL) {
        *pdy = (yhi - ylo)/h - (3.0*a*a - 1.0)/6.0*h*y2lo + (3.0*b*b - 1.0)/6.0*h*y2hi;
    }
    return a*ylo + b*yhi + ((a*a*a - a)*y2lo + (b*b*b - b)*y2hi)*(h*h)/6.0;
}

/*
 * Interpolate u(rho, s) in the table: cubic splines in s along the two bracketing rows for u and
 * u1 = du/dlogrho, then cubic Hermite interpolation in logrho between the rows. If pduds is not
 * NULL, du/ds is returned as well (needed by the root finder).
 *
 * rho and s must be inside the table; the callers check this.
 */
double tillCubicIntD(TILLMATERIAL *material, double rho, double s, double *pduds)
{
    const size_t offU = offsetof(TILL_LOOKUP_ENTRY, u), offU2 = offsetof(TILL_LOOKUP_ENTRY, udv2);
    const size_t offU1 = offsetof(TILL_LOOKUP_ENTRY, u1), offU12 = offsetof(TILL_LOOKUP_ENTRY, u1dv2);
    double u0, u1, d0, d1, e;
    int i;

    tillLocateRho(material, rho, &i, &e);

    if (pduds == NULL) {
        u0 = tillSplineIntFieldD(material, s, i, offU, offU2, NULL);
        u1 = tillSplineIntFieldD(material, s, i+1, offU, offU2, NULL);
        d0 = tillSplineIntFieldD(material, s, i, offU1, offU12, NULL);
        d1 = tillSplineIntFieldD(material, s, i+1, offU1, offU12, NULL);
        return tillHermite(u0, u1, d0, d1, material->dlogrho, e, NULL);
    } else {
        /* d/ds of the Hermite is the Hermite of the d/ds of its inputs (the weights depend on e only). */
        double u0s, u1s, d0s, d1s;
        u0 = tillSplineIntFieldD(material, s, i, offU, offU2, &u0s);
        u1 = tillSplineIntFieldD(material, s, i+1, offU, offU2, &u1s);
        d0 = tillSplineIntFieldD(material, s, i, offU1, offU12, &d0s);
        d1 = tillSplineIntFieldD(material, s, i+1, offU1, offU12, &d1s);
        *pduds = tillHermite(u0s, u1s, d0s, d1s, material->dlogrho, e, NULL);
        return tillHermite(u0, u1, d0, d1, material->dlogrho, e, NULL);
    }
}

double tillCubicInt(TILLMATERIAL *material, double rho, double s)
{
    return tillCubicIntD(material, rho, s, NULL);
}

double tillFindUonIsentrope(TILLMATERIAL *material, double s, double rho)
{
    return tillCubicIntD(material, rho, s, NULL);
}

/*
 * Find the isentrope through (rho, u): returns its table variable s, or TILL_FAIL_VALUE if (rho, u)
 * is below the cold curve (s < 0) or above the table (s > s_max), or u is not finite. iOrder is
 * only used in messages.
 *
 * u(rho, s) is monotonically increasing in s, so the cell [s_j, s_j+1] containing the root is
 * found by bisection on the grid lines (cheap Hermite evaluations), and the root inside the cell
 * by a bisection-safeguarded Newton iteration on the full interpolant.
 */
double tillFindEntropyCurve(TILLMATERIAL *material, double rho, double u, int iOrder)
{
    const double h = material->dv;
    const double utol = tillUTol(material, u);
    double e, ulo, uhi, slo, shi, s, f, dfds;
    int i, jlo, jhi, iter;

    /* NaN would pass both range tests below and send the bisection astray. */
    if (!isfinite(u))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i) has a non-finite u.\n",
                         iOrder, rho, u, material->iMaterial);

    tillLocateRho(material, rho, &i, &e);

    jlo = 0;
    jhi = material->nTableV - 1;
    ulo = tillCubicIntGrid(material, i, e, jlo);
    uhi = tillCubicIntGrid(material, i, e, jhi);

    if (u < ulo - utol)
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i) is below the cold curve.\n",
                         iOrder, rho, u, material->iMaterial);
    if (u > uhi + utol)
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i) is above the table.\n",
                         iOrder, rho, u, material->iMaterial);

    /* Bisection on the grid lines: u(s_jlo) <= u < u(s_jhi). */
    while (jhi - jlo > 1) {
        int jm = (jlo + jhi)/2;
        double um = tillCubicIntGrid(material, i, e, jm);
        if (u < um) { jhi = jm; uhi = um; } else { jlo = jm; ulo = um; }
    }

    slo = jlo*h;
    shi = jhi*h;
    /* Start from the linear estimate inside the cell. */
    s = (uhi > ulo) ? slo + h*(u - ulo)/(uhi - ulo) : 0.5*(slo + shi);
    s = fmin(fmax(s, slo), shi);

    for (iter = 0; iter < 50; iter++) {
        f = tillCubicIntD(material, rho, s, &dfds) - u;
        if (fabs(f) <= utol) return s;
        /* Shrink the bracket (u is increasing in s). */
        if (f > 0.0) shi = s; else slo = s;
        if (shi - slo <= 1e-12*(1.0 + s)) return 0.5*(slo + shi);
        /* Newton step, bisection if it leaves the bracket or the slope is not positive. */
        if (dfds > 0.0) {
            double sn = s - f/dfds;
            s = (sn > slo && sn < shi) ? sn : 0.5*(slo + shi);
        } else {
            s = 0.5*(slo + shi);
        }
    }
    TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillFindEntropyCurve: Particle %i (rho= %g, u= %g, iMat= %i): root finder did not converge.\n",
                     iOrder, rho, u, material->iMaterial);
}

/*
 * u_cold(rho2), the end point of the cold isentrope (s = 0) at rho2, for the clamp in tillLookupU.
 * Inside the table (and below rho0, where it is zero) from the table; above rhomax the cold curve
 * is continued by integrating the ODE from (rhomax, u_cold(rhomax)), like any other isentrope
 * that leaves the table. Returns TILL_FAIL_VALUE for a NaN rho2 (from tillColdULookup).
 */
static double tillColdIsentropeU(TILLMATERIAL *material, double rho2)
{
    if (rho2 > material->rhomax*(1.0 + TILL_EDGE_TOL)) {
        double ucold = tillColdULookup(material, material->rhomax);
        if (tillFailed(ucold)) return TILL_FAIL_VALUE;
        return tillCalcU(material, material->rhomax, ucold, rho2);
    }
    return tillColdULookup(material, rho2);
}

/*
 * Calculate u2 on the isentrope through (rho1, u1) at rho2. Inside the lookup table this is an
 * interpolation; for rho outside [rhomin, rhomax] or u1 above the table the ODE is integrated
 * directly.
 *
 * If (rho1, u1) is below the cold curve there is no isentrope through it. The state is then
 * clamped to the cold curve, i.e. u2 = u_cold(rho2), the end point of the coldest isentrope
 * (s = 0) through (rho1, u_cold(rho1)). This is what ANEOSisentropicU and scvheosIsentropicU do
 * (T(rho1, u1) falls back to the lowest temperature and the isentrope of that state is followed);
 * the fallback is counted in nLookupFail and reported under EOSLIB_VERBOSE. Returns
 * TILL_FAIL_VALUE if u1 is not finite (or rho2 is invalid).
 */
double tillLookupU(TILLMATERIAL *material, double rho1, double u1, double rho2, int iOrder)
{
    double s;
    int iRet;

    iRet = tillIsInTable(material, rho1, u1);

    if (iRet == TILL_LOOKUP_OUTSIDE_VMIN) {
        material->nLookupFail++;
        TILL_WARN("tillLookupU: Particle %i (rho= %g, u= %g, iMat= %i) is below the cold curve, clamped to u_cold(rho2= %g).\n",
                  iOrder, rho1, u1, material->iMaterial, rho2);
        return tillColdIsentropeU(material, rho2);
    }

    /* Non-finite u1 (TILL_FAIL): there is no isentrope to follow, not even by direct integration. */
    if (iRet == TILL_FAIL)
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillLookupU: Particle %i (rho= %g, u= %g, iMat= %i) has a non-finite u.\n",
                         iOrder, rho1, u1, material->iMaterial);

    /* Check if the starting and end point are inside of the look up table. */
    if ((iRet == TILL_LOOKUP_SUCCESS) &&
        (rho2 >= material->rhomin*(1.0 - TILL_EDGE_TOL)) && (rho2 <= material->rhomax*(1.0 + TILL_EDGE_TOL)))
    {
        s = tillFindEntropyCurve(material, rho1, u1, iOrder);
        if (tillFailed(s)) return TILL_FAIL_VALUE;
        return tillFindUonIsentrope(material, s, rho2);
    }

    TILL_WARN("tillLookupU: (rho1= %g, u1= %g, rho2= %g, iMat= %i) outside of the lookup table (%s), doing direct integration.\n",
            rho1, u1, rho2, material->iMaterial, tillLookupStatusString(iRet));
    return tillCalcU(material, rho1, u1, rho2);
}

/*
 * Calculate u_cold(rho) = u(rho, s=0) from the lookup table. The cold curve is zero for all
 * rho < rho0 (the pressure cutoff keeps u constant on the cold isentrope in tension), including
 * densities below rhomin. Returns TILL_FAIL_VALUE for rho > rhomax.
 */
double tillColdULookup(TILLMATERIAL *material, double rho)
{
	// Make sure that the look up table is initialized.
	assert(material->Lookup != NULL);

    // In case of a liquid there is not cold term in the expanded states.
    if (rho < material->rho0)
        return 0.0;

    if (!(rho >= material->rhomin*(1.0 - TILL_EDGE_TOL)) || (rho > material->rhomax*(1.0 + TILL_EDGE_TOL)))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillColdULookup: rho= %g outside of the lookup table (iMat= %i, rhomin= %g, rhomax= %g).\n",
                         rho, material->iMaterial, material->rhomin, material->rhomax);

    {
        double e;
        int i;
        tillLocateRho(material, rho, &i, &e);
        return tillCubicIntGrid(material, i, e, 0);
    }
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