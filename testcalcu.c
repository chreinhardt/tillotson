/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author: Christian Reinhardt
 * Date:   24.09.2018
 *
 * Test the function tillCalcU(). First a table is generated and printed to a file, then 
 * intermediate isentropes are interpolated a particle's evolution along each isentrope is
 * calculated by directly solving the ODE. The results can be plotted with testlookupulogrho.py.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) (((i)*tillMat->nTableV) + (j))

void main(int argc, char **argv) {
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
//	int nTableRho = 1000;
//	int nTableV = 1000;
//	double rhomax = 25.0;
//	double vmax = 25.0;
	int nTableRho = 100;
	int nTableV = 100;
	double rho, v, u;
    double rho1, rho2;
	FILE *fp = NULL;
	int i, j;

	double k = 0.0;
	double l = 0.0;

#ifdef TILL_PRESS_NP
	fprintf(stderr, "TILL_PRESS_NP.\n");
#endif
	fprintf(stderr, "Initializing material...\n");

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");

	/* Solve ODE and splines */
	tillInitLookup(tillMat);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", tillMat->rhomax, tillMat->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", tillMat->nTableRho, tillMat->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", tillMat->drho, tillMat->dv);
	fprintf(stderr,"\n");

	rho = 0.0;
	v = 0.0;
	u = 0.0;

	/*
     * Print the look up table to a file first.
     */	
	fp = fopen("lookup.txt", "w");
	assert(fp != NULL);

	for (i=0; i<tillMat->nTableRho; i+=1)
	{
		rho = tillLookupRho(tillMat, i);
		fprintf(fp, "%15.7E", rho);

		for (j=0; j<tillMat->nTableV; j+=1)
		{
			u = tillMat->Lookup[INDEX(i, j)].u;
			fprintf(fp, "%15.7E", u);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	/*
     * Interpolate values between the isentropes.
     */
	fp = fopen("testsplint.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Interpolating isentropes...\n");

	for (i=0; i<tillMat->nTableRho-1; i+=1)
	{
        // Logarithmic spacing
		rho = tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho);
		fprintf(fp, "%15.7E", rho);

		for (j=0;j<tillMat->nTableV-1;j+=1)
		{
            // Linear spacing
            v = tillMat->dv*(j+0.5);
            u = tillCubicInt(tillMat, rho, v);

			fprintf(fp, "%15.7E", u);
		}

		fprintf(fp,"\n");
	}
	fclose(fp);

	/*
     * Evolve particles along the isentropes doing direct intergration of the ODE.
     */
	fp = fopen("testcalcu.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Evolve a particle along an isentrope (direct integration)...\n");

#if 0
    /* First along the isentropes in the lookup table from rhomin to rhomax. */
	for (j=0; j<tillMat->nTableV-1; j+=1)
	{
        rho1 = tillLookupRho(tillMat, 1);

        // Linear spacing
        v = tillLookupV(tillMat, j);

        u = tillCubicInt(tillMat, rho1, v);
        
        fprintf(fp, "%15.7E%15.7E", rho1, u);
        
        rho2 = tillLookupRho(tillMat, tillMat->nTableRho-2);

	    fprintf(stderr, "i= %i rho1=%g u1= %g rho2= %g v= %g\n", j, rho1, u, rho2, v);

        u = tillCalcU(tillMat, rho1, u, rho2);

        fprintf(fp, "%15.7E%15.7E\n", rho2, u);
	}
#endif
    /* Then along the isentropes in the lookup table from rhomax to rhomin. */
	for (j=0; j<tillMat->nTableV-1; j+=1)
	{
        rho1 = tillLookupRho(tillMat, tillMat->nTableRho-2);

        // Linear spacing
        v = tillMat->dv*(j+0.5);

        u = tillCubicInt(tillMat, rho1, v);
        
        fprintf(fp, "%15.7E%15.7E", rho1, u);
        
        rho2 = tillLookupRho(tillMat, 1);

	    fprintf(stderr, "i= %i rho1=%g u1= %g rho2= %g v= %g\n", j, rho1, u, rho2, v);

        u = tillCalcU(tillMat, rho1, u, rho2);

        fprintf(fp, "%15.7E%15.7E\n", rho2, u);
	}
	fclose(fp);

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(tillMat);
}
