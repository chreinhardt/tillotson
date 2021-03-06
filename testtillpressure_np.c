/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     04.03.2019
 * Modified: 05.03.2019
 *
 * Calculate the pressure of a material with tillPressureSoundNP() and check if it agrees with
 * tillPressure().
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

void main(int argc, char **argv) {
    // Tillotson EOS library
	TILLMATERIAL *tillmat;
    double dKpcUnit = 2.06701e-13;
    double dMsolUnit = 4.80438e-08;
    double vmax = 100.0;
    int iMat = GRANITE;
    // Define the size of the grid
    double rhomin = 1e-4;
	double rhomax = 10.0;
	double umax = 25.0;
	double umin = 0.0;
    double T_min = 0.0;
    double T_max = 1e4;
	int nRho = 1000;
	int nU = 1000;
	int nV = 1000;
	int nT = 100;
    double drho, dU, dT;
	double rho, u, T;
	FILE *fp = NULL;
	int i, j;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
	tillInitLookup(tillmat, nRho, nV, rhomin, rhomax, vmax);
    fprintf(stderr, "Done.\n");

    /*
     * Zoom in to the expanded states.
     */
    rhomax = 1.1*tillmat->rho0;
    umax = 1.1*tillmat->us2;

    fprintf(stderr, "rhomin= %g rhomax= %g umin= %g umax= %g\n", rhomin, rhomax, umin, umax);

    drho = (rhomax-rhomin)/(nRho-1);
    dU = (umax-umin)/(nU-1);
    dT = (T_max-T_min)/(nT-1);

	/*
	 * Print P on a rho x u grid.
	 */	
	fp = fopen("testtillpressure_np.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tillPressureSoundNP()...\n");

	/* Print a rho x u grid. */
	for (i=0; i<nU; i+=1)
	{
		u = umin + i*dU;
		for (j=0; j<nRho; j+=1)
		{
		    rho = rhomin + j*drho;
			fprintf(fp," %15.7E", tillPressureSoundNP(tillmat, rho, u, NULL));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testtillpressure.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tillPressure()...\n");

	/* Print a rho x u grid. */
	for (i=0; i<nU; i+=1)
	{
        u = umin + i*dU;

		for (j=0; j<nRho; j+=1)
		{
            rho = rhomin + j*drho;

			fprintf(fp," %15.7E", tillPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testtillpressure_np_diff.txt", "w");
	assert(fp != NULL);

	/* print a rho x u grid. */
	for (i=0; i<nU; i+=1)
	{
        u = umin + i*dU;

		for (j=0; j<nRho; j+=1)
		{
            rho = rhomin + j*drho;

            double P;
            P = fabs(tillPressureSoundNP(tillmat, rho, u, NULL)-tillPressure(tillmat, rho, u));

            if (P > 0.0) {
                fprintf(fp," %2i", 0);
            } else {
                fprintf(fp," %2i", 1);
            }
#if 0
            fprintf(fp," %15.7E", tillPressureSoundNP(tillmat, rho, u, NULL)-
                    tillPressure(tillmat, rho, u));
#endif
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

    /* Plot where P is negative. */
	fp = fopen("testtillpressure_np_region.txt", "w");
	assert(fp != NULL);

	/* print a rho x u grid. */
	for (i=0; i<nU; i+=1)
	{
        u = umin + i*dU;

		for (j=0; j<nRho; j+=1)
		{
            rho = rhomin + j*drho;
            if (tillPressureSoundNP(tillmat, rho, u, NULL) <= 0.0)
//            if (tillPressureSound(tillmat, rho, u, NULL) <= 0.0)
            {
                fprintf(fp," %2i", 0);
            } else {
                fprintf(fp," %2i", 1);
            }
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

    /*
     * Plot a few isotherms.
     */
	fp = fopen("testtillpressure_np_rhot.txt", "w");
	assert(fp != NULL);

	for (i=0; i<nRho; i+=1)
	{
        rho = rhomin + i*drho;

		fprintf(fp,"%15.7E", rho);

		for (j=0; j<nT; j+=1)
		{
            T = T_min + j*dT;
            u = eosURhoTemp(tillmat, rho, T);
			fprintf(fp," %15.7E", tillPressureSoundNP(tillmat, rho, u, NULL));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	tillFinalizeMaterial(tillmat);
}
