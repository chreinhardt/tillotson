/*
 * Calculate the generalized Gruneisen parameter
 *
 * Gamma(rho, u) = a + b/omega0
 *
 * along different isentropes for a material.
 *
 * Author:   Christian Reinhardt
 * Date:     01.02.2020
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillMat->nTableV) + (j))

double Gamma(TILLMATERIAL *tillMat, double rho, double u) {
    double eta;
    double omega0;

    assert(rho >= tillMat->rho0);
    eta = rho/tillMat->rho0;
    omega0 = u/(tillMat->u0*eta*eta)+1.0;

    return(tillMat->a + tillMat->b/omega0);
}

int main(int argc, char **argv) {
    TILLMATERIAL *tillMat;
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
    //int nSteps;
    int iMat;
    int i, j;

#if 0
    if (argc != 5)
    {
        fprintf(stderr, "tillcalcgamma <iMat> <rho_max> <u_start> <nSteps>\n");
        exit(1);
    }

    iMat = atoi(argv[1]);
    rho_max = atof(argv[2]);
    u_inital = atof(argv[3]);
    nSteps = atoi(argv[4]);

    assert(iMat >= 0);
#endif

    iMat = 1;

    /* Initalize the material and calculate the isentropes. */
    tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);

    rhomin = 0.99*tillMat->rho0;
    rhomax = 10000.0*tillMat->rho0;
    vmax = 25.0;
    //assert(rhomax > rhomin);


    tillInitLookup(tillMat, nTableRho, nTableV, rhomin, rhomax, vmax);

    rho = tillMat->rho0;
    u = 0.0;
    fprintf(stderr, "Gamma(rho=%15.7E, u=%15.7E)=%15.7E\n", rho, u, Gamma(tillMat, rho, u));

    i = 999;
    j = 999;
    rho = tillLookupRho(tillMat, i);
    u = tillMat->Lookup[INDEX(i, j)].u;
    fprintf(stderr, "Gamma(rho=%15.7E, u=%15.7E)=%15.7E\n", rho, u, Gamma(tillMat, rho, u));

    exit(1);
    fprintf(stderr, "i=%i\n", tillMat->n);
    /*
     * Calculate Gamma(rho, u) for the cold curve.
     */
    for (i=tillMat->n+1; i<tillMat->nTableRho; i+=10)
    {
        fprintf(stderr, "i=%i\n", i);
        rho = tillLookupRho(tillMat, i);
        printf("%15.7E", rho/tillMat->rho0);

        j = 0;
        u = tillMat->Lookup[INDEX(i, j)].u;

        printf("%15.7E", Gamma(tillMat, rho, u));

        printf("\n");
    }
#if 0
    /*
     * Calculate Gamma(rho, u) for different isentropes.
     */
    for (i=0; i<tillMat->nTableRho; i+=10)
    {
        rho = tillLookupRho(tillMat, i);
        printf("%15.7E", rho);

        for (j=0; j<tillMat->nTableV; j+= 10)
        {
            u = tillMat->Lookup[INDEX(i, j)].u;

            printf("%15.7E", Gamma(tillMat, rho, u));
        }
        printf("\n");
    }
#endif

    tillFinalizeMaterial(tillMat);

    return 0;
}
