/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     05.06.2018
 * Modified: 03.11.2018
 *
 * Calculate the pressure of an ideal gas using eosPressure() and compare to the analytic solution.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

//#define TILL_PRESS_NP

void main(int argc, char **argv)
{
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double rhomin = 1e-4;
	double rhomax = 100.0;
    double vmax = 100.0;
	double umax = 100.0;
    double drho, dv;
	double rho, rho_int, u, P, cs2, dPdrho, dPdu;
	int nRho = 100;
	int nV = 100;

	TILLMATERIAL *tillmat;
	FILE *fp1 = NULL;
	FILE *fp2 = NULL;

	int i = 0;
	int j = 0;

    drho = (rhomax-rhomin)/(nRho-1);
    dv = vmax/(nV-1);

    // Print version information
    fprintf(stderr,"Tillotson EOS library version: %s\n", TILL_VERSION_TEXT);

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(IDEALGAS, dKpcUnit, dMsolUnit);
	fprintf(stderr, "Done.\n");

	/*
	 * Determine P(rho, u) on a rho x u grid.
	 */	
	fp1 = fopen("testtillidealgas.pressure1.txt", "w");
	assert(fp1 != NULL);

	fp2 = fopen("testtillidealgas.pressure2.txt", "w");
	assert(fp2 != NULL);
	
	/* Print a rho x u grid. */
	for (i=0; i<nRho; i++)
	{
		rho = rhomin + i*drho;

        fprintf(stderr, "rho= %g\n", rho);
		for (j=0; j<nV; j++)
		{
			u = j*dv;

			fprintf(fp1, " %15.7E", eosPressure(tillmat, rho, u));
			fprintf(fp2, " %15.7E", (5.0/3.0-1.0)*rho*u);
		}
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	
    fclose(fp1);
	fclose(fp2);

	/*
	 * Determine cs2(rho, u) on a rho x u grid.
	 */	
	fp1 = fopen("testtillidealgas.soundspeed1.txt","w");
	assert(fp1 != NULL);

	fp2 = fopen("testtillidealgas.soundspeed2.txt","w");
	assert(fp2 != NULL);
	
	/* Print a rho x u grid. */
	for (i=0; i<nRho; i++)
	{
		rho = rhomin + i*drho;

		for (j=0; j<nV; j++)
		{
			u = j*dv;
            P = eosPressureSound(tillmat, rho, u, &cs2);
			fprintf(fp1, " %15.7E", cs2);
			fprintf(fp2, " %15.7E", (5.0/3.0)*(5.0/3.0-1.0)*u);
        }
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	
    fclose(fp1);
	fclose(fp2);

	/*
	 * Determine dPdrho(rho, u) on a rho x u grid.
	 */	
	fp1 = fopen("testtillidealgas.dPdrho1.txt","w");
	assert(fp1 != NULL);

	fp2 = fopen("testtillidealgas.dPdrho2.txt","w");
	assert(fp2 != NULL);
	
	/* Print a rho x u grid. */
	for (i=0; i<nRho; i++)
	{
		rho = rhomin + i*drho;

		for (j=0; j<nV; j++)
		{
			u = j*dv;

			fprintf(fp1, " %15.7E", eosdPdrho(tillmat, rho, u));
			fprintf(fp2, " %15.7E", (5.0/3.0-1.0)*u);
		}
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	
    fclose(fp1);
	fclose(fp2);

	/*
	 * Determine dPdu(rho, u) on a rho x u grid.
	 */	
	fp1 = fopen("testtillidealgas.dPdu1.txt","w");
	assert(fp1 != NULL);

	fp2 = fopen("testtillidealgas.dPdu2.txt","w");
	assert(fp2 != NULL);
	
	/* Print a rho x u grid. */
	for (i=0; i<nRho; i++)
	{
		rho = rhomin + i*drho;

		for (j=0; j<nV; j++)
		{
			u = j*dv;

			fprintf(fp1, " %15.7E", eosdPdu(tillmat, rho, u));
			fprintf(fp2, " %15.7E", (5.0/3.0-1.0)*rho);
		}
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	
    fclose(fp1);
	fclose(fp2);


	tillFinalizeMaterial(tillmat);
}
