/*
 * Calculate the Hugoniot curve for a given material of the Tillotson EOS.
 *
 * Author:   Christian Reinhardt
 * Date:     05.02.2020
 * Modified: 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillMat->nTableV) + (j))

int main(int argc, char **argv)
{
    // Tillotson library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	int nTableRho = 1000;
	int nTableV = 1000;
    double rhomax = 100.0;
    double rhomin = 1e-4;
    double vmax = 100.0;
    double rho;
    double u;
    double P;
    int i, j;

    if (TILL_VERSION_MAJOR < 3)
    {
        fprintf(stderr, "Wrong version of the Tillotson library (%s).\n", TILL_VERSION_TEXT);
        exit(1);
    }

	fprintf(stderr, "Initializing material...\n");

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
    tillInitLookup(tillMat, nTableRho, nTableV, rhomin, rhomax, vmax);

    fprintf(stderr, "Done.\n");

    for (i=tillMat->n+1; i<tillMat->nTableRho; i++)
    {
        rho = tillLookupRho(tillMat, i);
        printf("%15.7E", rho);

        j = 0;
        u = tillMat->Lookup[INDEX(i, j)].u;

        printf("%15.7E", tillHugoniotPressure(tillMat, rho));
        printf("%15.7E", tillPressure(tillMat, rho, u));

        printf("\n");
    }

	tillFinalizeMaterial(tillMat);

    return 0;
}
