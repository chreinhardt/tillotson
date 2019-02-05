/*
 * Calculate u_cold(rho) for a given material.
 *
 * Author:   Christian Reinhardt
 * Date:     05.02.2019
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillMat->nTableV) + (j))

void main(int argc, char **argv)
{
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double rhomin = 1e-4;
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;
    char MatName[256];
    double rho, u;
    int iMat;

    if (argc != 3)
    {
        fprintf(stderr,"Usage: till_calc_u_cold <rho> <iMat>\n");
        exit(1);
	}

    rho = atof(argv[1]);
    assert(rho >= 0.0);

    iMat = atoi(argv[2]);
    assert(iMat > 0);

    // Initialize the Tillotson library
	tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    tillInitLookup(tillMat, nTableRho, nTableV, rhomin, rhomax, vmax);


    tilliMatString(tillMat, &MatName);

    u = tillCubicInt(tillMat, rho, 0.0);

    // Print output
    printf("iMat= %i (%s): rho= %g u= %g\n", iMat, MatName, rho, u);

    tillFinalizeMaterial(tillMat);
}
