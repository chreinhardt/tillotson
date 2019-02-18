/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     17.02.2019
 *
 * Calculate the pressure of a material with eosPressureRhoT() and check if it agrees with
 * eosPressure() for a given temperature.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

int main(int argc, char **argv) {
    TILLMATERIAL *tillmat;
    double dKpcUnit = 2.06701e-13;
    double dMsolUnit = 4.80438e-08;
    double rho, u;
    double drho, du;
    double rhomin = 1e-4;
    double rhomax = 100.0;
    double vmax = 1200.0;
    double umin = 0.0;
    double umax = 100;
    int nRho = 1000;
    int nU = 1000;
    double P1;
    double P2;
    double T;
    FILE *fp = NULL;
    int i = 0;
    int j = 0;

    if (TILL_VERSION_MAJOR < 3) {
        fprintf(stderr, "Wrong version of the Tilltson library (%s).\n", TILL_VERSION_TEXT);
        exit(1);
    }

    fprintf(stderr, "Initializing material...\n");
    tillmat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
    tillInitLookup(tillmat, nRho,  nU, rhomin, rhomax, vmax);

    fprintf(stderr, "Done.\n");

    /*
     * Print P on a rho x u grid.
     */	
    fp = fopen("testeospressurerhot.txt", "w");
    assert(fp != NULL);

    fprintf(stderr, "Printing grid (nRho=%i, nU= %i) using eosPressureRhoT()...\n", nRho, nU);

    drho = (rhomax-rhomin)/(nRho-1);
    du = (umax-umin)/(nU-1);

    /* Print a rho x u grid. */
    for (i=0; i<nRho; i+=1)
    {
        rho = rhomin + i*drho;

        for (j=0; j<nU; j+=1)
        {
            u = umin + j*du;

            // Check if (rho, u) are below the cold curve.
            if (eosIsBelowColdCurve(tillmat, rho, u)) {
                fprintf(fp," %15.7E", -1e10);
            } else {
                T = eosTempRhoU(tillmat, rho, u);
                P1 = eosPressure(tillmat, rho, u);
                P2 = eosPressureRhoT(tillmat, rho, T);

                fprintf(fp," %15.7E", P1-P2);
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    tillFinalizeMaterial(tillmat);

    return 0;
}
