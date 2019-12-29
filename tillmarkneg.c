/*
 * Calculate for which values of (rho, u) the pressure is negative or the sound speed is imaginary.
 *
 * Author:   Christian Reinhardt
 * Date:     29.12.2019
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

int main(int argc, char **argv)
{
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double P, cs2;
    double rho, u;
    double rho_min, rho_max;
    double u_min, u_max;
    double drho;
    double du;
    int nRho = 100;
    int nU = 100;
    FILE *fp1;
    FILE *fp2;
    int i, j;

    // Initialize the Tillotson library
	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
    
    rho_min = 0.0;
    rho_max = 1.1*tillMat->rho0;
    u_min = 0.0;
    u_max = tillMat->us2;
    
    drho = (rho_max-rho_min)/(nRho-1);
    du = (u_max-u_min)/(nU-1);

    fp1 = fopen("press_neg.txt", "w");
    fp2 = fopen("cs_neg.txt", "w");

    for (i=0; i<nRho; i++) {
        for (j=0; j<nU; j++) {
            rho = rho_min + i*drho;
            u = u_min + j*du;

            P = tillPressureSoundNP(tillMat, rho, u, &cs2);

            if (P < 0.0)
                fprintf(fp1, "%15.7E %15.7E\n", rho, u);

            if (cs2 < 0.0)
                fprintf(fp2, "%15.7E %15.7E\n", rho, u);

        }
    }

    fclose(fp1);
    fclose(fp2);

    tillFinalizeMaterial(tillMat);

    return 0;
}
