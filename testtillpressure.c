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

double Pressure(TILLMATERIAL *mat, double rho, double u)
{
    double eta;
    double mu;
    
    double P1,P2,P;

    eta = rho/mat->rho0;
    mu = eta - 1.0;

    if (rho >= mat->rho0)
    {
        /*
         * Condensed states.
         */
        P = (mat->a+mat->b/(u/(mat->u0*eta*eta)+1))*rho*u + mat->A*mu + mat->B*mu*mu;
        return (P);
    } else if (u < mat->us) {
        /*
         * Expanded cold states (same as above).
         */
        P = (mat->a+mat->b/(u/(mat->u0*eta*eta)+1))*rho*u + mat->A*mu + mat->B*mu*mu;
        return (P);
    } else if (u > mat->us2) {
        /*
         * Expanded hot states.
         */
        P = mat->a*rho*u + (mat->b*rho*u/(u/(mat->u0*eta*eta)+1) +
            mat->A*mu*exp(-mat->beta*(mat->rho0/rho-1.0)))*exp(-mat->alpha*pow(mat->rho0/rho-1.0,2.0));
        return (P);
    } else {
        /*
         * Expanded intermediate states.
         */
        P1 = (mat->a+mat->b/(u/(mat->u0*eta*eta)+1))*rho*u + mat->A*mu + mat->B*mu*mu;
        P2 = mat->a*rho*u + (mat->b*rho*u/(u/(mat->u0*eta*eta)+1) +
             mat->A*mu*exp(-mat->beta*(mat->rho0/rho-1.0)))*exp(-mat->alpha*pow(mat->rho0/rho-1.0,2.0));

        P = ((u-mat->us)*P2 + (mat->us2-u)*P1)/(mat->us2-mat->us);
        return (P);
    }
}

void main(int argc, char **argv)
{
	/*
     * Calculate the pressure one from tillPressure() and once directly from the
     * Tillotson EOS to check that both agree.
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
	fp = fopen("testtillpressure.txt","w");
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
			fprintf(fp," %15.7E", tillPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("testtillpressure2.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using Pressure()...\n");

	/* Print a rho x u grid. */
	for (i=0;i<tillmat->nTableRho;i+=1)
	{
		rho = i*tillmat->drho;
		for (j=0;j<tillmat->nTableV;j+=1)
		{
			// v = j*tillmat->dv
			u = j*tillmat->dv;
			fprintf(fp," %15.7E", Pressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	tillFinalizeMaterial(tillmat);
}
