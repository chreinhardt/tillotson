/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     01.10.2018
 *
 * Calculate the pressure of a material with eosPressure() and check if it agrees with
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
	double rho, u;

	TILLMATERIAL *tillmat;
	FILE *fp = NULL;

	int i = 0;
	int j = 0;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);
	fprintf(stderr, "Done.\n");

	/*
	 * Print P on a rho x u grid.
	 */	
	fp = fopen("testeospressure.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tillPressure()...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
//			fprintf(fp," %15.7E", eosPressure(tillmat, rho, u)-tillPressure(tillmat, rho, u));
			fprintf(fp," %15.7E", tillPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testeospressure2.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using eosPressure()...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
			fprintf(fp," %15.7E", eosPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testeospressure3.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid for the ideal gas EOS...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
			fprintf(fp," %15.7E", eosPressure(NULL, rho, u)-(5.0/3.0-1.0)*rho*u);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	tillFinalizeMaterial(tillmat);
}
