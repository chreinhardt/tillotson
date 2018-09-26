/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Test the bicubic interpolation tillCubicInt() in the lookup table. First a table is generated
 * and printed to a file, then values between the isentropes are interpolated. The results can be
 * plotted with testsplint.py.
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
	double rhomin = TILL_RHO_MIN;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	int nTableRho = 1000;
	int nTableV = 1000;
//	double rhomax = 25.0;
//	double vmax = 25.0;
//	int nTableRho = 100;
//	int nTableV = 100;
	double rho, v, u;
	struct lookup *isentrope;
	FILE *fp = NULL;
	int i, j;

	double k = 0.0;
	double l = 0.0;

#ifdef TILL_PRESS_NP
	fprintf(stderr, "TILL_PRESS_NP.\n");
#endif
	fprintf(stderr, "Initializing material...\n");

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
	
	fprintf(stderr, "Initializing the look up table...\n");

	/* Solve ODE and splines */
	tillInitLookup(tillMat, nTableRho, nTableV, rhomin, rhomax, vmax);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", tillMat->rhomax, tillMat->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", tillMat->nTableRho, tillMat->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", tillMat->drho, tillMat->dv);
	
	rho = 0.0;
	v = 0.0;
	u = 0.0;

	/*
     * Print the look up table to a file first.
     */	
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	for (i=0; i<tillMat->nTableRho; i+=1)
	{
		rho = tillLookupRho(tillMat, i);
		fprintf(fp, "%15.7E", rho);

		for (j=0;j<tillMat->nTableV;j+=1)
		{
			// v = j*tillMat->dv
			u = tillMat->Lookup[INDEX(i, j)].u;
			fprintf(fp, "%15.7E", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	/*
     * Interpolate values between the isentropes.
     */
	fp = fopen("testsplint.txt","w");
	assert(fp != NULL);

	for (i=0; i<tillMat->nTableRho-1; i+=1)
	{
        // Logarithmic spacing
		rho = tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho);
        fprintf(stderr, "rho= %15.7E, ", rho);
		fprintf(fp, "%15.7E", rho);

		for (j=0;j<tillMat->nTableV-1;j+=1)
		{
            v = tillMat->dv*(j+0.5);
            fprintf(stderr, "v= %15.7E\n", v);
            u = tillCubicInt(tillMat, rho, v);

			fprintf(fp, "%15.7E", u);
		}

		fprintf(fp,"\n");
	}

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(tillMat);
}
