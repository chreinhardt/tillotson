/*
 * This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

//#define TILL_PRESS_NP

void main(int argc, char **argv)
{
	/*
     * Calculate the sound speed and check how the minimum sound speed affects
     * the result.
     */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 100.0;
	double rho, u, cs;
#if !defined(TILL_PRESS_NP) || !defined(TILL_PRESS_MELOSH)
    double u_max;
#endif
	int nTableRho = 100;
	int nTableV = 100;

	TILLMATERIAL *tillmat;
	FILE *fp1 = NULL;
	FILE *fp2 = NULL;

	int i = 0;
	int j = 0;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	fprintf(stderr, "Done.\n");

	/*
	 * Print cs on a rho x u grid.
	 */	
	fp1 = fopen("soundspeed_cutoff1.txt", "w");
	assert(fp1 != NULL);

    fp2 = fopen("soundspeed_cutoff2.txt", "w");
	assert(fp2 != NULL);

	fprintf(stderr, "Printing grid using tillPressureSound()...\n");

	/* Print a rho x u grid. */
	for (i=0; i<tillmat->nTableRho; i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0; j<tillmat->nTableV; j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
            tillPressureSound(tillmat, rho, u, &cs);
			fprintf(fp1, " %15.7E", cs);

            if (cs < sqrt(tillmat->A/tillmat->rho0))
                cs = sqrt(tillmat->A/tillmat->rho0);
            
			fprintf(fp2, " %15.7E", cs);
		}
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	fclose(fp1);
	fclose(fp2);

	tillFinalizeMaterial(tillmat);
}
