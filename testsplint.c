/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     19.09.2018
 * Modified: 26.09.2018
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
//	double rhomax = 100.0;
//	double vmax = 1200.0;
	double rhomax = 25.0;
	double vmax = 25.0;
	int nTableRho = 100;
	int nTableV = 100;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
//	int nTableRho = 1000;
//	int nTableV = 1000;
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

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");

	/* Solve ODE and splines */
	tillInitLookup(tillMat);
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
        /*
         * Careful, its better to use Lookup[i, j].rho.
		 * rho = i*tillMat->drho;
         */
        rho = tillMat->Lookup[INDEX(i, j)].rho;
		rho = tillMat->rhomin+i*tillMat->drho;
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

	for (i=2; i<tillMat->nTableRho-1; i+=1)
	{
		rho = tillMat->rhomin+(i + 0.5)*tillMat->drho;
        fprintf(stderr, "rho= %15.7E, ", rho);
		rho = tillMat->Lookup[INDEX(i, 0)].rho+tillMat->drho*0.5;
        fprintf(stderr, "rho= %15.7E\n", rho);

		fprintf(fp, "%15.7E", rho);
		for (j=0;j<tillMat->nTableV-1;j+=1)
		{
            v = tillMat->dv*(j+0.5);
            u = tillCubicInt(tillMat, rho, v);

			fprintf(fp, "%15.7E", u);
		}

		fprintf(fp,"\n");
	}

#if 0
	/* Interpolate values between the isentropes */
	for (i=0;i<tillMat->nTableRho-2;i+=1)
	{
		// Middle of the interval (i,i+1)
		//rho = (i + 0.5)*tillMat->drho;
		l = 0.0;
		while (l < 0.9)
		{
			// Try
			rho = tillMat->rhomin+(i + l)*tillMat->drho;
			//rho = tillMat->Lookup[INDEX(i, 0)].rho;
			//rho += l*fabs((tillMat->Lookup[INDEX(i, 0)].rho-tillMat->Lookup[INDEX(i+1, 0)].rho));

			printf("%g", rho);
			for (j=0;j<tillMat->nTableV-1;j+=1)
			{
				// Middle of the interval (i,i+1)
				// v = (j + 0.5)*tillMat->dv;
				k = 0.5;
				while (k < 0.9)
				{
					// This does not work for non uniform steps in v
					//v = (j + k)*tillMat->dv;
					v = tillMat->vmax/pow(tillMat->nTableV-1,n)*pow(j+k,n);

					v = tillMat->vmax/tillMat->dv*(j+k);

					u = tillCubicInt(tillMat, rho, v);

					//fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);
					printf("  %.8g", u);
					k+=0.5;
				}
			}
		printf("\n");
		l+=0.5;
		}
	}
#endif

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(tillMat);
}
