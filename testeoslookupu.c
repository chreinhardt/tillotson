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
     * Calculate rho(P, u) for an ideal gas from the Tillotson EOS library and
     * compare the results to the analytic expression.
     */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
    double vmax = 100.0;
	double umax = 100.0;
	double rho, rho1, rho2, u, P;
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
	 * Determine rho(P, u) on a P x u grid.
	 */	
	fp = fopen("testeosrhopu.txt","w");
	assert(fp != NULL);

	for (i=0; i<tillmat->nTableRho; i++)
    {
        rho = i*tillmat->drho;

		for (j=0; j<tillmat->nTableV; j++)
		{
			u = j*tillmat->dv;
            P = eosPressure(tillmat, rho, u);

//            fprintf(stderr, "rho= %g, u= %g, P= %g\n", rho, u, P);

            P /= fabs(P);
            if (P > 0.0) P = 1e5;

			fprintf(fp, " %15.7E", P);
            if (P < 0.0)
                fprintf(stderr, "rho= %g, u= %g, P= %g\n", rho, u, P);

		}

		fprintf(fp, "\n");
	}

    exit(1);
    /*
     * Starting from rho=0 or u=0 causes problems as P=0 for all rho.
     */
	for (i=1; i<tillmat->nTableRho; i++)
    {
        rho = i*tillmat->drho;

#if 0
        if (rho >= 1.0/tillmat->b)
        {
            fprintf(stderr, "i= %i: rho= %g, b= %g, 1/b= %g rho*b= %g\n", i, rho, tillmat->b, 1.0/tillmat->b, rho*tillmat->b);
            continue;
        }
#endif

		for (j=1; j<tillmat->nTableV; j++)
		{
			u = j*tillmat->dv;
            P = eosPressure(tillmat, rho, u);

            fprintf(stderr, "rho= %g, u= %g, P= %g\n", rho, u, P);

            rho1 = eosRhoPU(tillmat, P, u);
            rho2 = P/((tillmat->dConstGamma-1.0)*u+tillmat->b*P);

			fprintf(fp, " %15.7E", rho1-rho2);

            if (fabs(rho1-rho2) > 1e-10)
                fprintf(stderr, "rho1= %15.7E rho2= %15.7E\n", rho1, rho2);
		}
		fprintf(fp, "\n");
	}
	
    fclose(fp);

	tillFinalizeMaterial(tillmat);
}
