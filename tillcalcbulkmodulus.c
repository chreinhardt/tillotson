/*
 * Calculate the isentropic bulk modulus
 *
 * Ks=rho*dP/drho (s=const.)
 *
 * along different isentropes for a material.
 *
 * Author:   Christian Reinhardt
 * Date:     19.11.2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

void main(int argc, char **argv) {
	TILLMATERIAL *tillmat;
    /* Unit convertion factors */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    /* Lookup table */
	double rhomax = 100.0;
    double rhomin = 1e-4;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, u;
    double K, cs2, P;
    int nRho, nU;
	int iMat;
    int i, j;

    nRho = 100;
    nU = 10;

	if (argc != 2)
	{
        fprintf(stderr, "tillcalcbulkmodulus <iMat>\n");
		exit(1);
	}

	iMat = atoi(argv[1]);
	assert(iMat >= 0);

    /* Initalize the material and calculate the isentropes. */
	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    tillInitLookup(tillmat, nTableRho, nTableV, rhomin, rhomax, vmax);

    /*
     * Calculate the bulk modulus for different isentropes.
     */
    for (i=0; i<tillmat->nTableRho; i+=10)
    {
        rho = tillLookupRho(tillmat, i);
        printf("%15.7E", rho);

        for (j=0; j<tillmat->nTableV; j+= 10)
        {
            u = tillmat->Lookup[INDEX(i, j)].u;

            P = tillPressureSound(tillmat, rho, u, &cs2);
            K = rho*cs2;

            printf("%15.7E", K);
        }
        printf("\n");
    }

	tillFinalizeMaterial(tillmat);
}
