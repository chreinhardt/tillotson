/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author: Christian Reinhardt
 * Date:   06.04.2019
 *
 * Test if eosRhoPTemp can bracket the root for all pressures and temperatures on a grid.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) (((i)*granite->nTableV) + (j))

int main(int argc, char **argv) {
    // Tillotson EOS
    TILLMATERIAL *tillmat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
    double rhomin = 1e-4;
	double vmax = 1200.0;
	int nRho = 1000;
	int nV = 1000;
    int iMat = GRANITE;
    // Grid
    double P_min = 0.0;
    double P_max = 1e3;
    double T_min = 0.0;
    double T_max = 1e3;
    int nP = 1000;
    int nT = 1000;
	double rho, u;
    double P, T;
    double dP, dT;
	double rhoPT;
	FILE *fp = NULL;
	int i, j;

    // Initialize the Tillotson EOS
	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
	tillInitLookup(tillmat, nRho, nV, rhomin, rhomax, vmax);
	fprintf(stderr, "Done.\n");
	
    dP = (P_max-P_nin)/(nP-1);
    dT = (T_max-T_min)/(nT-1);

	/*
	 * Print rho on a P x T grid.
	 */	
	fp = fopen("testtillrhoptemp.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid...\n");

	for (i=0; i<nT; i+=1)
	{
		T = T_min + i*dT;

		for (j=0; j<nP; j+=1)
		{
		    P = P_min + j*dP;
			fprintf(fp," %15.7E", tillRhoPTemp(tillmat, P, T));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}
