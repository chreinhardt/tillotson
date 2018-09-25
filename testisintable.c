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
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 100;
	int nTableV = 100;
	double rho, v, u;
	double umax;
    char ErrorString[256]; 
    int iRet;
	FILE *fp = NULL;
	int i = 0;
	int j = 0;
	int n = 1;

#ifdef TILL_PRESS_NP
	fprintf(stderr, "TILL_PRESS_NP.\n");
#endif
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

	
	/*
	 * Print the look up table to a file first.
	 */	
	fp = fopen("lookup.txt", "w");
	assert(fp != NULL);

	for (i=0; i<tillMat->nTableRho; i+=1)
	{
		rho = tillLookupRho(tillMat, i);
		fprintf(fp, "%15.7E",rho);

		for (j=0; j<tillMat->nTableV; j+=1)
		{
			u = tillMat->Lookup[INDEX(i, j)].u;
			fprintf(fp,"  %15.7E", u);
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
//    v -= 1e-8;
//    v -= 2.5*tillMat->dv;
    j = tillLookupIndexV(tillMat, v);
#if 0
	for (i=0; i<tillMat->nTableRho; i+=1)
	{
		rho = tillLookupRho(tillMat, i);
        u = tillMat->Lookup[INDEX(i, j)].u;

        fprintf(stderr, "\n");
        fprintf(stderr,"i= %i j= %i: Testing rho=%g, u=%g (v= %g)! index= %g = %i\n", i, j, rho, u, v, (log(rho)-log(tillMat->rhomin))/tillMat->dlogrho, (int) floor((log(rho)-log(tillMat->rhomin))/tillMat->dlogrho));


		iRet = tillIsInTable(tillMat, rho, u);

		if (iRet != TILL_LOOKUP_SUCCESS)
        {
            tillErrorString(iRet, ErrorString);
            fprintf(stderr,"i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) not in table (Error %s)!\n", i, j, rho, u, v, ErrorString);

            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);

//			fprintf(fp, "%15.7E %15.7E\n", rho, u);
        } else {
            fprintf(stderr, "i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) is in table.\n", i, j, rho, u, v);
            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);
            assert(0);
        }
	}
#endif
	for (i=0; i<tillMat->nTableRho; i+=1)
	{
        // Choose a point between the grid points (logarithmic spacing)
		rho = tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho);
        u = tillCubicInt(tillMat, rho, v);

        fprintf(stderr, "\n");
        fprintf(stderr,"i= %i j= %i: Testing rho=%g, u=%g (v= %g)! index= %g = %i\n", i, j, rho, u, v, (log(rho)-log(tillMat->rhomin))/tillMat->dlogrho, (int) floor((log(rho)-log(tillMat->rhomin))/tillMat->dlogrho));


		iRet = tillIsInTable(tillMat, rho, u);

		if (iRet != TILL_LOOKUP_SUCCESS)
        {
            tillErrorString(iRet, ErrorString);
            fprintf(stderr,"i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) not in table (Error %s)!\n", i, j, rho, u, v, ErrorString);

            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);

//			fprintf(fp, "%15.7E %15.7E\n", rho, u);
        } else {
            fprintf(stderr, "i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) is in table.\n", i, j, rho, u, v);
            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);
            assert(0);
        }
	}
	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(tillMat);
}
