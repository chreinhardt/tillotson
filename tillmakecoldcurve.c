/*
 * Generate the cold curve u_cold(rho) for a given material.
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
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;

//	double vmax = 28.0;
	double rho, u;
	int nTableRho = 1000;
	int nTableV = 1000;

    char MatName[256];
    
    int iMat;

	int i = 0;
	int j = 0;

	TILLMATERIAL *tillMat;
	TILL_LOOKUP_ENTRY *isentrope;

	fprintf(stderr, "Initializing material:");
	tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);

	fprintf(stderr, "%s\n", MatName);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(tillMat);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", tillMat->rhomax, tillMat->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", tillMat->nTableRho, tillMat->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", tillMat->drho, tillMat->dv);
	fprintf(stderr,"\n");

#if 0
    /*
     * Generate the cold curve from the lookup table.
     */
	FILE *fp = NULL;

	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	/* Print the lookup table to a file. */
	for (i=0;i<tillMat->nTableRho;i+=1)
	{
        /*
         * This is not the same as taking rho(i,0).
         */
		rho = i*tillMat->drho;

        u = tillMat->Lookup[INDEX(i, 0)].u;
        fprintf(fp, "%15.7E %15.7E\n", rho, u);
	}
	fclose(fp);
#endif

    /*
     * Print some header.
     */
    printf("# rho             u\n");

    /*
     * Solve for the cold curve.
     */
	isentrope = tillSolveIsentrope(tillMat, 0);

	for (j=0; j<tillMat->nTableRho; j++)
	{
		printf("%15.7E %15.7E\n", isentrope[j].rho, isentrope[j].u);
	}

    fprintf(stderr, "rhomax= %g umax=%g\n", isentrope[tillMat->nTableRho-1].rho, isentrope[tillMat->nTableRho-1].u);

    tillFinalizeMaterial(tillMat);
}
