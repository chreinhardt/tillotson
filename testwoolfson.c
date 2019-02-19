/*
 * Calculate the correction factor f_ij for two given materials.
 *
 * Author:   Christian Reinhardt
 * Date:     21.08.2018
 * Modified: 13.02.2019
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"
#include "woolfson.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

#define TILL_PRESS_NP

int main(int argc, char **argv)
{
    // Tillotson library
	TILLMATERIAL **tillmat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	int nTableRho = 100;
	int nTableV = 100;
    double vmax = 100.0;
    double f_ij;
    double rhomin, rhomax, umin, umax;
    int i;

    if (TILL_VERSION_MAJOR < 3)
    {
        fprintf(stderr, "Wrong version of the Tillotson library (%s).\n", TILL_VERSION_TEXT);
        exit(1);
    }

    rhomin = 1e-4;
    rhomax = 10000.0;

    umin = 0.0;
    umax = 20000.0;

    tillmat = (TILLMATERIAL **) calloc(2, sizeof(TILLMATERIAL *));

	fprintf(stderr, "Initializing material...\n");

    tillmat[0] = tillInitMaterial(IRON, dKpcUnit, dMsolUnit);
	tillmat[1] = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);

    tillInitLookup(tillmat[0], nTableRho, nTableV, rhomin, rhomax, vmax);
    tillInitLookup(tillmat[1], nTableRho, nTableV, rhomin, rhomax, vmax);

    fprintf(stderr, "Done.\n");

    // Example for a 5ME model (30% iron, 70% granite)
    double rho1 = 50.3198;
    double u1 = 19.3885;
    double rho2 = 30.3324;
    double u2 = 38.6285;

    double P1 = tillPressure(tillmat[0], rho1, u1);
    double T1 = tillTempRhoU(tillmat[0], rho1, u1);

    double P2 = tillPressure(tillmat[1], rho2, u2);
    double T2 = tillTempRhoU(tillmat[1], rho2, u2);

    fprintf(stderr, "\n");
    fprintf(stderr, "material 1: imat= %i rho= %g u= %g P= %g T= %g\n", tillmat[0]->iMaterial, rho1, u1, P1, T1);
    fprintf(stderr, "material 2: imat= %i rho= %g u= %g P= %g T= %g\n", tillmat[1]->iMaterial, rho2, u2, P2, T2);
    fprintf(stderr, "\n");

    for (i=0; i<100000; i++) {
        f_ij = CalcWoolfsonCoeff(tillmat[0], tillmat[1], P1, T1);
    }

    fprintf(stderr, "Direct calculation: f_ij=rho1/rho2= %g (rho1/rho2= %g) f_ji=1/f_ij=rho2/rho1= %g (rho2/rho1= %g)\n", f_ij, rho1/rho2, 1.0/f_ij, rho2/rho1);
    fprintf(stderr, "\n");

	tillFinalizeMaterial(tillmat[0]);
	tillFinalizeMaterial(tillmat[1]);

    return 0;
}
