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
//#include "tillsplint.h"

/* Header file for the Numerical Recipes routines */
#include "nr/nrcubicspline.h"

/*
 * Basic functions:
 *
 * tillInitSplines: initialise the cubic splines in u and u1.
 *
 * tillSplineIntU: do a cubic spline interpolation of u in v.
 *
 * tillSplineIntU1: do a cubic spline interpolation of u1 in v.
 *
 * tillCubicInt: interpolate a value (i,i+1) in logrho and (j,j+1) in v (uses Joachim's 2D
 *               interpolator).
 */

void tillInitSplines(TILLMATERIAL *material)
{
	/*
     * Calculate the second derivatives for u and u1 in v. For this we use routines from the
     * Numerical Recipes.
	 */
	tillInitSplineU1(material);
	tillInitSplineU(material);
}

void tillInitSplineV(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u in v.
	*/
	int i,j,k,n;
	double yp1, ypn;	// dudv(v=0) and dudv(v=n-1)
	double p,qn,sig,un,*u;

	n = material->nTableV;
	/* Allocate memory for temporary array */
	u = malloc(n*sizeof(double));
	
	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	for (i=0; i<material->nTableRho; i++)
	{
		/*
		** Set up splines in u(v) for a given rho.
		*/

		/* Set b.c */	
		if (yp1 > 0.99e30)
			material->Lookup[TILL_INDEX(i,0)].udv2=u[0]=0.0;
		else {
			material->Lookup[TILL_INDEX(i,0)].udv2 = -0.5;
//			u[0]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
			/* (CR) This part has to be changed if dv is not constant! */
			u[0]=(3.0/(material->dv))*((material->Lookup[TILL_INDEX(i,1)].u-material->Lookup[TILL_INDEX(i,0)].u)/(material->dv)-yp1);
		}

		/* Calculate second derivatives of u in v */
		for (j=1;j<=n-2;j++) {
//			sig=(x[j]-x[j-1])/(x[j+1]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			sig=material->dv/(2.0*material->dv);
			p=sig*material->Lookup[TILL_INDEX(i,j-1)].udv2+2.0;
			material->Lookup[TILL_INDEX(i,j)].udv2=(sig-1.0)/p;
//			u[j]=(y[j+1]-y[j])/(x[j+1]-x[j]) - (y[j]-y[j-1])/(x[j]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(material->Lookup[TILL_INDEX(i,j+1)].u-material->Lookup[TILL_INDEX(i,j)].u)/(material->dv) - (material->Lookup[TILL_INDEX(i,j)].u-material->Lookup[TILL_INDEX(i,j-1)].u)/(material->dv);
//			u[j]=(6.0*u[j]/(x[j+1]-x[j-1])-sig*u[j-1])/p;
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(6.0*u[j]/(2.0*material->dv)-sig*u[j-1])/p;
		}

		/* Set b.c. */
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
//			Careful with index!! This array ends at n-1, so this code is already using the right indices.			
//			un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
			/* (CR) This part has to be changed if dv is not constant! */
			un=(3.0/(material->dv))*(ypn-(material->Lookup[TILL_INDEX(i,n-1)].u-material->Lookup[TILL_INDEX(i,n-2)].u)/(material->dv));
		}
		material->Lookup[TILL_INDEX(i,n-1)].udv2=(un-qn*u[n-2])/(qn*material->Lookup[TILL_INDEX(i,n-2)].udv2+1.0);
		for (k=n-2;k>=0;k--)
			material->Lookup[TILL_INDEX(i,k)].udv2=material->Lookup[TILL_INDEX(i,k)].udv2*material->Lookup[TILL_INDEX(i,k+1)].udv2+u[k];
	}
}

void tillInitSplineU1(TILLMATERIAL *material)
{
	/*
     * Calculate the second derivatives for u1 in v. For this we use routines from the
     * Numerical Recipes.
     */
	int i,j,k,n;
	double yp1, ypn;	// du1dv(v=0) and d1udv(v=n-1)
	double p,qn,sig,un,*u;

	n = material->nTableV;
	/* Allocate memory for temporary array */
	u = malloc(n*sizeof(double));
		
	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	for (i=0; i<material->nTableRho; i++)
	{
		/*
		** Set up splines in u(v) for a given rho.
		*/

		/* Set b.c */	
		if (yp1 > 0.99e30)
			material->Lookup[TILL_INDEX(i,0)].u1dv2=u[0]=0.0;
		else {
			material->Lookup[TILL_INDEX(i,0)].u1dv2 = -0.5;
//			u[0]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
			/* (CR) This part has to be changed if dv is not constant! */
			u[0]=(3.0/(material->dv))*((material->Lookup[TILL_INDEX(i,1)].u-material->Lookup[TILL_INDEX(i,0)].u)/(material->dv)-yp1);
		}

		/* Calculate second derivatives of u in v */
		for (j=1;j<=n-2;j++) {
//			sig=(x[j]-x[j-1])/(x[j+1]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			sig=material->dv/(2.0*material->dv);
			p=sig*material->Lookup[TILL_INDEX(i,j-1)].u1dv2+2.0;
			material->Lookup[TILL_INDEX(i,j)].u1dv2=(sig-1.0)/p;
//			u[j]=(y[j+1]-y[j])/(x[j+1]-x[j]) - (y[j]-y[j-1])/(x[j]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(material->Lookup[TILL_INDEX(i,j+1)].u-material->Lookup[TILL_INDEX(i,j)].u)/(material->dv) - (material->Lookup[TILL_INDEX(i,j)].u-material->Lookup[TILL_INDEX(i,j-1)].u)/(material->dv);
//			u[j]=(6.0*u[j]/(x[j+1]-x[j-1])-sig*u[j-1])/p;
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(6.0*u[j]/(2.0*material->dv)-sig*u[j-1])/p;
		}

		/* Set b.c. */
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
//			Careful with index!! This array ends at n-1, so this code is already using the right indices.			
//			un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
			/* (CR) This part has to be changed if dv is not constant! */
			un=(3.0/(material->dv))*(ypn-(material->Lookup[TILL_INDEX(i,n-1)].u-material->Lookup[TILL_INDEX(i,n-2)].u)/(material->dv));
		}
		material->Lookup[TILL_INDEX(i,n-1)].u1dv2=(un-qn*u[n-2])/(qn*material->Lookup[TILL_INDEX(i,n-2)].u1dv2+1.0);
		for (k=n-2;k>=0;k--)
			material->Lookup[TILL_INDEX(i,k)].u1dv2=material->Lookup[TILL_INDEX(i,k)].u1dv2*material->Lookup[TILL_INDEX(i,k+1)].u1dv2+u[k];
	}
}

void tillInitSplineU(TILLMATERIAL *material)
{
	tillInitSplineV(material);
}

/*
 * Do a cubic spline interpolation of u in v.
 *
 * (CR) 10.08.2016: Integrated splint() into the function, allowing to use the TILL_LOOKUP_ENTRY
 * structure for the interpolation.
 */
double tillSplineIntv(TILLMATERIAL *material, double v, int ilogrho)
{
	int klo,khi,k;
	double h,b,a;
	double u;

	/* Use v=k*material->dv */	
	klo = tillLookupIndexV(material, v);
	khi = klo+1;

	h = (khi-klo)*material->dv;
	assert(h != 0.0);
	assert(h == material->dv);

	a=(khi*material->dv-v)/h;
	b=(v-klo*material->dv)/h;

	u=a*material->Lookup[TILL_INDEX(ilogrho,klo)].u+b*material->Lookup[TILL_INDEX(ilogrho,khi)].u+((a*a*a-a)*material->Lookup[TILL_INDEX(ilogrho,klo)].udv2+(b*b*b-b)*material->Lookup[TILL_INDEX(ilogrho,khi)].udv2)*(h*h)/6.0;

	return u;
}

/*
 * Do a cubic spline interpolation of u in v.
 */
double tillSplineIntU(TILLMATERIAL *material, double v, int ilogrho)
{
	/*
	** Do a cubic spline interpolation of u in v.
	*/
	int klo,khi,k;
	double h,b,a;
	double u;
	
    /* Use v=k*material->dv */
	klo = tillLookupIndexV(material, v);
	khi = klo+1;

	h = (khi-klo)*material->dv;
	assert(h != 0.0);
	assert(h == material->dv);

	a=(khi*material->dv-v)/h;
	b=(v-klo*material->dv)/h;

	u=a*material->Lookup[TILL_INDEX(ilogrho,klo)].u+b*material->Lookup[TILL_INDEX(ilogrho,khi)].u+((a*a*a-a)*material->Lookup[TILL_INDEX(ilogrho,klo)].udv2+(b*b*b-b)*material->Lookup[TILL_INDEX(ilogrho,khi)].udv2)*(h*h)/6.0;

	return u;
}

/*
 * Do a cubic spline interpolation of u1 (du/dlogrho) in v.
 */
double tillSplineIntU1(TILLMATERIAL *material, double v, int ilogrho)
{
	int klo,khi,k;
	double h,b,a;
	double u1;
	
	/* Use v=k*material->dv */
	klo = tillLookupIndexV(material, v);
	khi = klo+1;

	h = (khi-klo)*material->dv;
	assert(h != 0.0);
	assert(h == material->dv);

	a=(khi*material->dv-v)/h;
	b=(v-klo*material->dv)/h;
	
	u1=a*material->Lookup[TILL_INDEX(ilogrho,klo)].u1+b*material->Lookup[TILL_INDEX(ilogrho,khi)].u1+((a*a*a-a)*material->Lookup[TILL_INDEX(ilogrho,klo)].u1dv2+(b*b*b-b)*material->Lookup[TILL_INDEX(ilogrho,khi)].u1dv2)*(h*h)/6.0;

	return u1;
}

#ifdef TILL_DEBUG_SPLINT
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv)
{
	/*
	** Do a cubic spline interpolation in rho.
	*/
	int klo,khi,k;
	double h,b,a;
	double u;

    /// CR: Not yet adapted for logrho.
    assert(0);

    /* Use rho = rhomin + i*drho. */
	klo = tillLookupIndexRho(material, rho);
	khi = klo+1;

	h = (khi-klo)*material->drho;
	assert(h != 0.0);
	assert(h == material->drho);

	a=(khi*material->drho-rho)/h;
	b=(rho-klo*material->drho)/h;
	
	return a*material->Lookup[TILL_INDEX(klo,iv)].u+b*material->Lookup[TILL_INDEX(khi,iv)].u+((a*a*a-a)*material->Lookup[TILL_INDEX(klo,iv)].udrho2+(b*b*b-b)*material->Lookup[TILL_INDEX(khi,iv)].udrho2)*(h*h)/6.0;
}
#endif

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
	double *ce;
	
	assert(intvalues != NULL);

	/* Allocate memory */
	ce = malloc(4*sizeof(double));
	assert(ce != NULL);

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
	
	// free memory
	free(ce);
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
	assert(i < material->nTableRho-1);

	// vint is between v[j] and v[j+1]
	j = tillLookupIndexV(material, vint);
	assert(j < material->nTableV-1);

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
	B = (vint-material->dv)/material->dv;
	
	dudv[0] = (material->Lookup[TILL_INDEX(i, j+1)].u-material->Lookup[TILL_INDEX(i, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j+1)].udv2;
	dudv[1] = (material->Lookup[TILL_INDEX(i+1, j+1)].u-material->Lookup[TILL_INDEX(i+1, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j+1)].udv2;

	/* Calculate dudvdlogrho from u1dv2 */
	dudvdlogrho[0] = (material->Lookup[TILL_INDEX(i, j+1)].u1-material->Lookup[TILL_INDEX(i, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j+1)].u1dv2;
	dudvdlogrho[1] = (material->Lookup[TILL_INDEX(i+1, j+1)].u1-material->Lookup[TILL_INDEX(i+1, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j+1)].u1dv2;

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

/* Prasenjits root finder */
double brent(double (*func)(TILLMATERIAL *,double,double,double),TILLMATERIAL *material,double a,double b,double rho,double u,double tol,int iOrder);

double tillFindUonIsentrope(TILLMATERIAL *material,double v,double rho)
{
	double iv,irho,u;
	
    return tillCubicInt(material, rho, v);
}

double denergy(TILLMATERIAL *material,double v,double rho,double u)
{
	return (tillFindUonIsentrope(material,v,rho)-u);
}

/* Find isentrope for a given rho and u */
double tillFindEntropyCurve(TILLMATERIAL *material,double rho,double u,int iOrder)
{
	double tol=1e-6;
    double v_eps = V_EPS;
    
    /* Note that using vmax rather than vmax-eps can cause an error in the brent root finder. */
	return brent(denergy,material,0.0,material->vmax-v_eps,rho,u,tol,iOrder);
}

double tillLookupU(TILLMATERIAL *material, double rho1, double u1, double rho2, int iOrder)
{
	/* Calculates u2 for a given rho1, u2 and rho2. */
	double v, u;
    int iRet;

    iRet = tillIsInTable(material, rho1, u1);

	/* Check if the starting and end point are inside of the look up table */
	if ((iRet == TILL_LOOKUP_SUCCESS) && (rho2 > material->rhomin) && (rho2 < material->rhomax))
	{
		/* Interpolate using the look up table */
		v = tillFindEntropyCurve(material, rho1, u1, iOrder);
		u = tillFindUonIsentrope(material, v, rho2);
	} else {
        /* Do direct integration unless the value is below the cold curve. */
#ifdef TILL_VERBOSE
		fprintf(stderr,"tillLookupU: value (%g, %g, %g) outside of look up table, doing direct integration (Error: %i).\n",
                rho1, u1, rho2, iRet);
#endif
        assert(iRet != TILL_LOOKUP_OUTSIDE_VMIN);
		u = tillCalcU(material, rho1, u1, rho2);
	}

	return u;
}

double eosLookupU(TILLMATERIAL *material, double rho1, double u1, double rho2, int iOrder)
{
	/* Calculates u2 for a given rho1, u2, rho2. */

    if (material->iMaterial == IDEALGAS)
    {
        /*
         * For the ideal gas EOS there is an analytic expression.
         */
//        return(u1*pow(rho2/rho1, material->dConstGamma-1.0));
        return(u1*pow(rho2/rho1*(1.0-material->b*rho1)/(1.0-material->b*rho2), material->dConstGamma-1.0));
    }  else {
        return tillLookupU(material, rho1, u1, rho2, iOrder);
    }
}

/*
 * Calculate u_cold(rho) from the lookup table.
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

    if ((rho < material->rhomin) || (rho >= material->rhomax))
    {
        fprintf(stderr, "tillColdULookup: rho= %g outside of the lookup table.\n", rho);
        assert(0);
    }

    return tillCubicInt(material, rho, 0.0);
}

/*
 * Return the index i, so that rho_i and rho_i+1 bracket rho in the lookup table.
 */
int tillLookupIndexRho(TILLMATERIAL *material, double rho)
{
    int i;

    // Assume uniform spacing in rho.
	i = floor((rho-material->rhomin)/material->drho);

    // Check if logrho is outside of the table.
    if (i < 0)
        return -1;

    if (i >= material->nTableRho-1)
        return (material->nTableRho-1);

    return i;
}

/*
 * Return the index i, so that log(rho_i) and log(rho_i+1) bracket log(rho) in the lookup table.
 *
 * Note: The index of rho(i) is i-1. Since the value of rho is still bracketed by i-1 and i this
 *       should be fine.
 */
int tillLookupIndexLogRho(TILLMATERIAL *material, double logrho)
{
    int i;

    // Assume uniform spacing in log(rho).
    i = (int) floor((logrho-log(material->rhomin))/material->dlogrho);

    // Check if logrho is outside of the table.
    if (i < 0)
        return -1;

    if (i >= material->nTableRho-1)
        return (material->nTableRho-1);

    return i;
}

/*
 * Return the index j, so that v_j and v_j+1 bracket v in the lookup table.
 */
int tillLookupIndexV(TILLMATERIAL *material, double v)
{
    int j;

    // Assume uniform spacing in v.
    j = (int) floor(v/material->dv);

    // Check if v is outside of the table.
    if (j < 0)
        return -1;

    if (j >= material->nTableV)
        return material->nTableV;

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
 * Return v(j).
 */
double tillLookupV(TILLMATERIAL *material, int j)
{
    // Assume uniform spacing in v.
    return (j*material->dv);
}
