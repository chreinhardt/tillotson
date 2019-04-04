/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     04.03.2019
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
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    // Define the size of the grid
    double rhomin = 1e-4;
	double rhomax = 100.0;
    double vmax = 100.0;
	double umax = 100.0;
	double umin = 0.0;
	int nRho = 1000;
	int nU = 1000;
    double drho, dU;

	double rho;
    double u;

	TILLMATERIAL *tillmat;
	FILE *fp = NULL;

	int i = 0;
	int j = 0;

    drho = (rhomax-rhomin)/(nRho-1);
    dU = (umax-umin)/(nU-1);

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
	fprintf(stderr, "Done.\n");

	/*
	 * Print P on a rho x u grid.
	 */	
	fp = fopen("testtillpressure_np.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tillPressureSoundNP()...\n");

	/* Print a rho x u grid. */
	for (i=0; i<nRho; i+=1)
	{
		rho = rhomin + i*drho;
		for (j=0; j<nU; j+=1)
		{
			u = umin + j*dU;
			fprintf(fp," %15.7E", tillPressureSoundNP(tillmat, rho, u, NULL));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testtillpressure.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tillPressure()...\n");

	/* Print a rho x u grid. */
	for (i=0; i<nRho; i+=1)
	{
		rho = rhomin + i*drho;
		for (j=0; j<nU; j+=1)
		{
			u = umin + j*dU;
			fprintf(fp," %15.7E", tillPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testtillpressure_np_diff.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid for the ideal gas EOS...\n");

	/* Print a rho x u grid. */
	for (i=0; i<nRho; i+=1)
	{
		rho = rhomin + i*drho;

		for (j=0; j<nU; j+=1)
		{
			u = umin + j*dU;
			fprintf(fp," %15.7E", tillPressureSoundNP(tillmat, rho, u, NULL)-
                    tillPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	tillFinalizeMaterial(tillmat);
}
