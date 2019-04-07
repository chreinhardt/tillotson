/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     06.04.2019
 * Modified: 07.04.2019
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
	double rho_max = 100.0;
    double rho_min = 1e-4;
	double v_max = 1200.0;
	int nRho = 1000;
	int nV = 1000;
    int iMat = GRANITE;
    // Grid
    double P_min = 1e-8;
    double P_max = 1e2;
    double T_min = 0.0;
    double T_max = 1e5;
    int nP = 1000;
    int nT = 1000;
	double rho, u;
    double P, T;
    double dP, dT, drho;
	FILE *fp = NULL;
	int i, j;

    // Initialize the Tillotson EOS
	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
	tillInitLookup(tillmat, nRho, nV, rho_min, rho_max, v_max);
	fprintf(stderr, "Done.\n");
	
    dP = (P_max-P_min)/(nP-1);
    dT = (T_max-T_min)/(nT-1);
    
    // Zoom in
	rho_max = 10.0;
    rho_min = 1e-4;
    drho = (rho_max-rho_min)/(nRho-1);

	/*
	 * Print rho on a P x T grid.
	 */	
	fp = fopen("testrhoptemp.txt", "w");
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

    /*
     * Compare the difference in pressure when using rho(P, T).
     */
	fp = fopen("testrhoptemp_diff.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid...\n");

	for (i=0; i<nT; i+=1)
	{
		T = T_min + i*dT;

		for (j=0; j<nP; j+=1)
		{
		    P = P_min + j*dP;

			rho = tillRhoPTemp(tillmat, P, T);

			fprintf(fp," %15.7E", eosPressureRhoT(tillmat, rho, T)-P);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

    /*
     * Compare the difference in rho due to root finding.
     */
	fp = fopen("testrhoptemp_rhot.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid...\n");

	for (i=0; i<nT; i+=1)
	{
		T = T_min + i*dT;

		for (j=0; j<nRho; j+=1)
		{
		    rho = rho_min + j*drho;
            P = eosPressureRhoT(tillmat, rho, T);

            // Avoid weird densities if the root finding fails.
            if (P <= 0.0) {
                printf("rho = %g T= %g P= %g\n", rho, T, P);
                fprintf(fp," %3i", -1);
            } else {
                if (fabs(tillRhoPTemp(tillmat, P, T)-rho) < 1e-3) {
                    fprintf(fp," %3i", 0);
                } else {
                    fprintf(fp," %3i", 1);
                }

//                fprintf(fp," %15.7E", tillRhoPTemp(tillmat, P, T)-rho);
            }
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(tillmat);

    return 0;
}
