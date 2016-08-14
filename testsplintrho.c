/*
 ** This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) (((i)*granite->nTableMax) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy. We store the second derivative
	** of u with respect to rho in udv2 just for debugging
	** splint(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;

	double rho, u;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	tillInitSplineRho(granite);

	fprintf(stderr, "Done.\n");

	rho = 0.0;
	u = 0.0;

	j = 51;
	 	
	for (j=0;j<granite->nTableMax;j++)
	{
		for (i=0;i<granite->nTableMax;i++)
		{
			//rho = granite->Lookup[INDEX(i,j)].rho;
			// Choose a point in the middle of the interval
			rho = granite->Lookup[INDEX(i,j)].rho + 0.5*granite->delta;
			u = tillSplineInterpolation(granite, rho, j);

			printf("%g %g\n", rho, u);
		}
		printf("\n");
	}

	tillFinalizeMaterial(granite);
}
