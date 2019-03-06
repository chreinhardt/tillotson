/*
 * Calculate P(rho, u) for a material.
 *
 * Author:   Christian Reinhardt
 * Date:     01.10.2018
 * Modified: 06.03.2019
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

int main(int argc, char **argv) {
    // Tillotson EOS library
    TILLMATERIAL *tillmat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    int nTableRho = 1000;
    int nTableV = 1000;
    double rho_min = 1e-4;
    double rho_max = 100.0;
    double v_max = 1200.0;
	double rho, T, P;
	int iMat = 0;

	if (argc != 4)
	{
			fprintf(stderr, "eoscalcpressure_rhot <iMat> <rho> <T>\n");
			exit(1);
	}

	iMat = atoi(argv[1]);
	assert(iMat >= 0);

	rho = atof(argv[2]);
	assert(rho >= 0.0);

	T = atof(argv[3]);
	assert(T >= 0.0);

	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    tillInitLookup(tillmat, nTableRho, nTableV, rho_min, rho_max, v_max);

	P = eosPressureRhoT(tillmat, rho, T);
	printf("%i %.8g %.8g %.8g\n", iMat, rho, T, P);

	tillFinalizeMaterial(tillmat);

    return 0;
}
