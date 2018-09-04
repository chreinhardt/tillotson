/*
 * Calculate the pressure along an isotherm.
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
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 10.0;
    double rhomin = 1e-2;
    double Tmax = 1e6;
    double Tmin = 1e-2;
    double drho;
        double dT;
    int nRho = 1000;
    int nT = 100;
    double vmax = 100.0;
	double umax = 100.0;
	double rho, u, P, T;
	TILLMATERIAL *tillmat;
    // To calculate u(rho, T) we need the lookup table.
	int nTableRho = 1000;
	int nTableV = 1000;
	int i = 0;
	int j = 0;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	fprintf(stderr, "Done.\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "Initializing lookup table...\n");
    tillInitLookup(tillmat);
	fprintf(stderr, "Done.\n");
	fprintf(stderr, "\n");
	/*
	 * Calculate P(rho, T=const) along different isotherms.
	 */	
	drho = (rhomax-rhomin)/(nRho-1);
    dT = (Tmax-Tmin)/(nT-1);

    /* Print a P(rho, T). */
	for (i=0; i<nRho; i++)
	{
        rho = rhomin + i*drho;
		printf("%15.7E", rho);

		for (j=0; j<nT; j++)
		{
            T = Tmin + j*dT;
			u = tillURhoTemp(tillmat, rho, T);

			printf("  %15.7E", eosPressure(tillmat, rho, u));
		}
		printf("\n");
	}

	tillFinalizeMaterial(tillmat);
}
