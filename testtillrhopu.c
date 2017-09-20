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
     * Calculate rho(P,u) for a given pressure and then see how well the result matches.
     */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
    double vmax = 100.0;
	double umax = 100.0;
	double rho, rho_int, u, P;
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
	 * Determine rho(P, u) on a rho x u grid.
	 */	
	fp1 = fopen("testtillrhopu1.txt","w");
	assert(fp1 != NULL);

	fp2 = fopen("testtillrhopu2.txt","w");
	assert(fp2 != NULL);
	
#if 0
    /* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;

        if (rho < tillmat->rho0) continue;

        u = 0.0;
        P = tillPressure(tillmat, rho, u);

        rho_int = tillRhoPU(tillmat, P, u);
    }
#endif
	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
        if (rho < tillmat->rho0) continue;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
            P = tillPressure(tillmat, rho, u);

			fprintf(fp1, " %15.7E", rho);
			fprintf(fp2, " %15.7E", tillRhoPU(tillmat, P, u));

            if (fabs(tillRhoPU(tillmat, P, u) - rho) > 1e-1) printf("%15.7E  %15.7E  %15.7E\n", rho, u, P);
		}
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	
    fclose(fp1);
	fclose(fp2);

	tillFinalizeMaterial(tillmat);
}
