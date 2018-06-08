/*
 * This is a simple program to test the Tillotson EOS library.
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
	/*
     * Calculate u(rho, P) for an ideal gas from the Tillotson EOS library and
     * compare the results to the analytic expression.
     */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
    double vmax = 100.0;
	double umax = 100.0;
	double rho, u, u1, u2, P;
	int nTableRho = 100;
	int nTableV = 100;

	TILLMATERIAL *tillmat;
	FILE *fp = NULL;

	int i = 0;
	int j = 0;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(IDEALGAS, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	fprintf(stderr, "Done.\n");

	/*
	 * Determine u(rho, P) on a rho x P grid.
	 */	
	fp = fopen("testeosurhop.txt", "w");
	assert(fp != NULL);


    for (i=1; i<tillmat->nTableRho; i++)
    {
        for (j=1; j<tillmat->nTableV; j++)
        {
            rho = i*tillmat->dv;
            u = j*tillmat->dv;
            P = eosPressure(tillmat, rho, u);

            fprintf(stderr, "rho= %g, u= %g, P= %g\n", rho, u, P);

            u1 = eosURhoP(tillmat, rho, P);
            u2 = P*(1.0-tillmat->b*rho)/((tillmat->dConstGamma-1.0)*rho);

            fprintf(fp, " %15.7E", u1-u2);

            if (fabs(u1-u2) > 1e-6)
                fprintf(stderr, "rho= %15.7E u1= %15.7E u2= %15.7E\n", rho, u1, u2);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");


    fclose(fp);

    tillFinalizeMaterial(tillmat);
}
