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

double dudrho(TILLMATERIAL *mat, double rho, double u)
{
    double P = (mat->dConstGamma-1.0)*rho*u/(1.0-mat->b*rho);
    return (P/(rho*rho));
//    return (eosPressure(mat, rho, u)/(rho*rho));
}

void main(int argc, char **argv)
{
	/*
     * Calculate u(rho) along an isentrope and compare to the result obtained
     * from eosLookupU().
     */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double rhomin = 1e-4;
	double rhomax = 10.0;
    double vmax = 100.0;
    double umin = 0.0;
	double umax = 100.0;
	double rho, rho1, rho2;
    double u, u1, u2;
    double drho, du;
    // Needed for RK2
    double k1u, k2u;
	int nTableRho = 100;
	int nTableV = 100;
    // The larger mean molecular weight the more points are needed for precission
    int nRho = 10000;
    int nU = 10;

	TILLMATERIAL *tillmat;
	FILE *fp = NULL;

	int i = 0;
	int j = 0;

	fprintf(stderr, "Initializing material...\n");

    // Do not convert to code units (does this really work?)
    dKpcUnit = 0.0;
	dMsolUnit = 0.0;

	tillmat = tillInitMaterial(IDEALGAS, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	fprintf(stderr, "Done.\n");

	/*
	 * Solve isentropes numerically for different inital u.
	 */
	fp = fopen("testeoslookupu.txt","w");
	assert(fp != NULL);

    rhomax = 1.0/tillmat->b;

    fprintf(stderr, "rhomax= %g b= %g\n", rhomax, tillmat->b);

    rhomax *= 0.9;

    drho = (rhomax - rhomin)/(nRho - 1);
    du = (umax - umin)/(nU - 1);

    for (i=0; i< nU; i++)
    {
        u = umin + i*du;

        rho = rhomin;
        
        for (j=0; j<nRho; j++)
        {
            /*
             * Midpoint Runga-Kutta (2nd order).
             */
            k1u = drho*dudrho(tillmat, rho, u);
		    k2u = drho*dudrho(tillmat, rho+0.5*drho, u+0.5*k1u);

		    u += k2u;
            rho += drho;
            
            /*
             * Note: All isentropes are written into two columns and then split
             * into different vectors in python.
             */
            fprintf(fp, "%15.7E  %15.7E\n", rho, u);
        }
    }

#if 0
    u = 0.1;
    rho = rhomin;

    for (j=0; j<nRho; j++)
    {
        /*
         * Midpoint Runga-Kutta (2nd order).
         */
        k1u = drho*dudrho(tillmat, rho, u);
		k2u = drho*dudrho(tillmat, rho+0.5*drho, u+0.5*k1u);

        u += k2u;
        rho += drho;
        
        fprintf(fp, "%15.7E  %15.7E\n", rho, u);
    }

    rho1 = rhomin;
    u1 = 0.1;

    rho2 = 0.9*rhomax;

    u2 = eosLookupU(tillmat, rho1, u1, rho2, i);

    printf("rho1= %15.7E u1= %15.7E rho2= %15.7E u2= %15.7E\n", rho1, u1, rho2, u2);
#endif
    /*
     * Use analytic solution.
     */

    for (i=0; i< nU; i++)
    {
        u1 = umin + i*du;
        rho1 = rhomin;

         printf("%15.7E  %15.7E\n", rho1, u1);
    
        for (j=0; j<nRho-1; j++)
        {
            rho2 = (j+1)*drho; 
            u2 = eosLookupU(tillmat, rho1, u1, rho2, i);
        
            printf("%15.7E  %15.7E\n", rho2, u2);

            rho1 = rho2;
            u1 = u2;
        }
    }

    fclose(fp);

	tillFinalizeMaterial(tillmat);
}
