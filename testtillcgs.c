/*
 * Verify that the Tillotson EOS is properly initialized if dKpcUnit and dMSolUnit are not set.
 *
 * Author:   Christian Reinhardt
 * Date:     16.03.2020
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillMat->nTableV) + (j))

int main(int argc, char **argv) {
    TILLMATERIAL *tillMat;
    /* Unit convertion factors */
    double dKpcUnit = 2.06701e-13;
    double dMsolUnit = 4.80438e-08;
    double rho, u;
    //int nSteps;
    int iMat;
    int i, j;

    iMat = 1;

    /* Initalize the material in cgs. */
    //tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    tillMat = tillInitMaterial(iMat, 0.0, 0.0);

    /* Print the material parameters. */
    tillPrintMat(tillMat);

    tillFinalizeMaterial(tillMat);

    return 0;
}
