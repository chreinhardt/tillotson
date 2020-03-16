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

int main(int argc, char **argv) {
    TILLMATERIAL *tillMat;
    /* Unit convertion factors */
    double rho, u;
    double P, cs;
    int iMat;

    iMat = 1;

    /* Initalize the material in cgs. */
    tillMat = tillInitMaterial(iMat, 0.0, 0.0);

    /* Print the material parameters. */
    tillPrintMat(tillMat);

    /* Print the pressure and sound speed at the reference state. */
    rho = tillMat->rho0;
    u = 0.0;

    /* Note that tillPressureSound returns cs2 not the sound speed. */
    P = tillPressureSound(tillMat, rho, u, &cs);

    printf("P= %15.7E cs= %15.7E K= %15.7E\n", P, sqrt(cs), rho*cs);

    tillFinalizeMaterial(tillMat);

    return 0;
}
