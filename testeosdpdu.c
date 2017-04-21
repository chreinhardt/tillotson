/*
 * This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	* Check, if dPdu agrees for eosdPdu and tilldPdu for
    * Tillotson materials and if eosdPdu agrees with
    * the analytic expression for an ideal gas.
    */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 100.0;
	double rho, u;
	int nTableRho = 100;
	int nTableV = 100;

	TILLMATERIAL *tillmat;
	FILE *fp = NULL;

	int i = 0;
	int j = 0;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	fprintf(stderr, "Done.\n");

	/*
	 * Print P on a rho x u grid.
	 */	
	fp = fopen("testeosdpdu1.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tilldPu()...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
			fprintf(fp," %15.7E", tilldPdu(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testeosdpdu2.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using eosdPu()...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			u = j*tillmat->dv;
			fprintf(fp," %15.7E", eosdPdu(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testeosdpdu3.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid for the ideal gas EOS...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			u = j*tillmat->dv;
			fprintf(fp," %15.7E", eosdPdu(NULL, rho, u)-(5.0/3.0-1.0)*rho);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	tillFinalizeMaterial(tillmat);
}
