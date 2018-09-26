/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     13.09.2018
 * Modified: 26.09.2018
 *
 * Generate the cold curve u_cold(rho) for a given material. Then interpolate using
 * tillColdULookup() and see if the results agree.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillMat->nTableV) + (j))

void main(int argc, char **argv) {
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomin = TILL_RHO_MIN;
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, u;
	FILE *fp = NULL;
    int i;

#ifdef TILL_PRESS_NP
	fprintf(stderr, "TILL_PRESS_NP.\n");
#endif
	fprintf(stderr, "Initializing material...\n");

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(tillMat, nTableRho, nTableV, rhomin, rhomax, vmax);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", tillMat->rhomax, tillMat->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", tillMat->nTableRho, tillMat->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", tillMat->dlogrho, tillMat->dv);
	fprintf(stderr,"\n");

	/*
	 * Print the cold curve to a file first.
	 */	
	fp = fopen("coldcurve.txt", "w");
	assert(fp != NULL);

	for (i=0; i<tillMat->nTableRho; i+=1)
	{
		rho = tillLookupRho(tillMat, i);
		fprintf(fp, "%15.7E",rho);

		u = tillMat->Lookup[INDEX(i, 0)].u;
		fprintf(fp, "%15.7E", u);
		
		fprintf(fp,"\n");
	}
	fclose(fp);

	/*
	 * Now print the interpolated data.
	 */	
    fp = fopen("testtillcoldu.txt", "w");
	assert(fp != NULL);
    
    fprintf(fp, "# rho             u\n");

	for (i=0; i<tillMat->nTableRho-1; i++)
	{
        // Choose a point between the grid points (logarithmic spacing)
		rho = tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho);
        
        u = tillColdULookup(tillMat, rho);

		fprintf(fp, "%15.7E%15.7E\n", rho, u);
	}
	fclose(fp);

    fprintf(stderr, "Done.\n");

    tillFinalizeMaterial(tillMat);
}
