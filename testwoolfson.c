/*
 * This is a simple program to test the functions provided in tillwoolfson.c.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"
#include "tillwoolfson.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

#define TILL_PRESS_NP

void main(int argc, char **argv)
{
	/*
     * Calculate a look up table for the coefficients f_ij for two materials
     * (e.g., iron and granite).
     */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
    double vmax = 100.0;
	double umax = 100.0;
	double P, T, rho1, rho2;
    double f_ij;
    double Pmin, Pmax, Tmin, Tmax;
    int nP, nT;
	int nTableRho = 100;
	int nTableV = 100;

	TILLMATERIAL **tillmat;
	FILE *fp = NULL;

	int i = 0;
	int j = 0;

    // Make a 100 x 100 grid.
    nP = 100;
    nT = 100;

    Pmin = 0.0;
    Pmax = 100.0;

    Tmin = 0.0;
    Tmax = 10000.0;

    tillmat = (TILLMATERIAL **) calloc(2, sizeof(TILLMATERIAL *));

	fprintf(stderr, "Initializing material...\n");

    tillmat[0] = tillInitMaterial(IRON, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	tillmat[1] = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);

    tillInitLookup(tillmat[0]);
    tillInitLookup(tillmat[1]);

    fprintf(stderr, "Done.\n");

	/*
	 * Calculate f_ij a P x T grid.
	 */	
	fp = fopen("testwoolfson.txt", "w");
	assert(fp != NULL);

	for (i=0; i<nP; i++)
    {
        P = Pmin + (Pmax - Pmin)/(nP-1)*i;

        fprintf(stderr, "P= %g\n", P);

		for (j=0; j<nT; j++)
		{
            T = Tmin + (Tmax - Tmin)/(nT-1)*j;
            f_ij = CalcWoolfsonCoeff(tillmat[0], tillmat[1], P, T);

			fprintf(fp, " %15.7E", f_ij);

#if 0
            u = j*tillmat->dv;
            P = eosPressure(tillmat, rho, u);

//            fprintf(stderr, "rho= %g, u= %g, P= %g\n", rho, u, P);

            P /= fabs(P);
            if (P > 0.0) P = 1e5;

			fprintf(fp, " %15.7E", P);
            if (P < 0.0)
                fprintf(stderr, "rho= %g, u= %g, P= %g\n", rho, u, P);
#endif
		}

		fprintf(fp, "\n");
	}

	
    fclose(fp);

	tillFinalizeMaterial(tillmat[0]);
	tillFinalizeMaterial(tillmat[1]);
}
