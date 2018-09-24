/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     13.09.2018
 * Modified: 24.09.2018 
 *
 * Test tillIsInTable() by checking, if a grid of points is in the lookup table. Especially the
 * values close to vmax can be problematic.
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
	/*
	** Debug tillCubicIntRho(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	// vmax=25.0, rhomax=100.0, nTableV=10, nTableRho=4000
	//int nTableMax = 1000;
	int nTableRho = 100;
	int nTableV = 100;
	double rho, v, u;
	double umax;
    char ErrorString[256]; 
    int iRet;

	int i = 0;
	int j = 0;
	int n = 1;

	TILLMATERIAL *tillMat;

	fprintf(stderr, "Initializing material...\n");

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, n);
	
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
	u = 0.0;

	/* Create an output file for the look up table */
	FILE *fp = NULL;
	
	/*
	 * Print the look up table to a file first.
	 */	
	fp = fopen("lookup.txt", "w");
	assert(fp != NULL);

	for (i=0; i<tillMat->nTableRho; i+=1)
	{
		rho = i*tillMat->drho;
		fprintf(fp,"%g",rho);
		for (j=0; j<tillMat->nTableV; j+=1)
		{
			u = tillMat->Lookup[INDEX(i, j)].u;
			fprintf(fp,"  %g", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

#if 0
	/*
     * Test if different points (rho, u) are in the lookup table.
     */
	fp = fopen("testisintable.txt", "w");
	assert(fp != NULL);

	umax = tillMat->Lookup[INDEX(tillMat->nTableRho-1, tillMat->nTableV-1)].u;
	fprintf(stderr, "umax=%g\n", umax);
	rho = -tillMat->rhomax*0.05;
	u = -10;

	while (rho < tillMat->rhomax*1.05)
	{
		u = -0.05*umax;
		while (u < 1.05*umax)
		{
			if (tillIsInTable(tillMat, rho, u) != TILL_LOOKUP_SUCCESS)
			{
//				fprintf(stderr,"rho=%g, u=%g not in table!\n", rho,u);
				fprintf(fp, "%15.7E %15.7E\n", rho, u);
			}
			u += 0.01*umax;
		}
		rho += tillMat->drho;
		j++;
	}
	fclose(fp);
#endif

    /* 
     * Now check if points on last isentrope are treated correctly.
     */
    v = tillLookupV(tillMat, tillMat->nTableV-1);
    j = tillMat->nTableV-1;

	for (i=0; i<tillMat->nTableRho; i+=1)
	{
        fprintf(stderr, "i= %i\n", i);
		rho = tillLookupRho(tillMat, i);
        u = tillMat->Lookup[INDEX(i, j)].u;

        fprintf(stderr,"i= %i j= %i: Testing rho=%g, u=%g (v= %g)!\n", i, j, rho, u, v);

		iRet = tillIsInTable(tillMat, rho, u);

		if (iRet != TILL_LOOKUP_SUCCESS)
        {
            tillErrorString(iRet, ErrorString);
            fprintf(stderr,"i= %i j= %i: rho=%g, u=%g (v= %g) not in table (Error %s)!\n", i, j, rho, u, v, ErrorString);
//			fprintf(fp, "%15.7E %15.7E\n", rho, u);
        }
	}

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(tillMat);
}
