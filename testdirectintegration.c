/*
 ** This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*granite->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	 * Debug function tillCalcU() that directly integrates u from rho1 to rho2.
	 */
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
#if 0
	double rhomax = 25.0;
	double vmax = 26.32;
	int nTableRho = 1000;
	int nTableV = 1000;
#endif
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;

	double rho, u, rho1, u1, rho2, u2;
	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"\n");

//#if 0
	/* Create an output file for the look up table */
	FILE *fp = NULL;

	/*
	 * Print the look up table to a file first.
	 */	
	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	/* Print the lookup table to a file. */
	for (i=0;i<granite->nTableRho;i+=1)
	{
		rho = i*granite->drho;
		fprintf(fp,"%g",rho);
		for (j=0;j<granite->nTableV;j+=1)
		{
			// v = j*granite->dv
			u = granite->Lookup[INDEX(i, j)].u;
			fprintf(fp,"  %g", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
//#endif

	rho1 = TILL_RHO_MIN;

	fprintf(stderr,"rho= %g\n",rho1);

	i = granite->nTableRho-5;

	i = granite->nTableRho-800;
	j = granite->nTableV-10;

	// Start from an isentrope that is in the lookup table
	rho1 = granite->Lookup[INDEX(i, j)].rho;
	u1 = granite->Lookup[INDEX(i, j)].u;

#if 0
//	rho2 = TILL_RHO_MIN*0.1;
	rho2 = 0.5*rho1;
	u2 = 0.0;
	u2 = tillCalcU(granite, rho1, u1, rho2);

	fprintf(stderr,"rho1= %g u1= %g rho2= %g u2= %g\n",rho1, u1, rho2, u2);
	printf("%15.7E%15.7E%15.7E%15.7E\n",rho1, u1, rho2, u2);
#endif
//#if 0	

	//rho2 = TILL_RHO_MIN*0.1;
	rho2 = rho1 - 1e-3;
	printf("%15.7E%15.7E%15.7E\n",rho1, u1,tillLookupU(granite,rho1,u1,rho2, 0));


	while (rho2 > 1e-5)
	{
		u2 = tillCalcU(granite, rho1, u1, rho2);
//		printf("%15.7E%15.7E%15.7E%15.7E\n",rho1, u1, rho2, u2);

		printf("%15.7E%15.7E%15.7E\n",rho2, u2,tillLookupU(granite,rho1,u1,rho2, 0));
		u1 = u2;
		rho1 = rho2;
		rho2 -= 1e-3;
	}
//#endif
	tillFinalizeMaterial(granite);
}
