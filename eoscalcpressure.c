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
	double rho, u, P;
	int iMat = 0;

	if (argc != 4)
	{
			fprintf(stderr, "eoscalcpressure <iMat> <rho> <u>\n");
			exit(1);
	}

	iMat = atoi(argv[1]);
	assert(iMat >= 0);

	rho = atof(argv[2]);
	assert(rho >= 0.0);

	u = atof(argv[3]);
	assert(u >= 0.0);

	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
	
	P = eosPressure(tillmat, rho, u);
	printf("%i %.8g %.8g %.8g\n", iMat, rho, u, P);

	tillFinalizeMaterial(tillmat);

    return 0;
}
