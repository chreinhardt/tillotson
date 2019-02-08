/*
 * Calculate the isentrope from (rho1, u1) to (rho2, u2) doing direct integration.
 *
 * Author:   Christian Reinhardt
 * Date:     07.02.2019
 * Modified: 08.02.2019
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*granite->nTableV) + (j))

void main(int argc, char **argv) {
    TILLMATERIAL *tillMat;
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double rhomin = 1e-4;
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, u;
    double drho;
    double rho1, u1, rho2, u2;
    int iMat;
    int nStep = 100;
	int i;

	if (argc != 5) {
		fprintf(stderr,"Usage: calcisentrope <rho1> <u1> <rho2> <iMat>\n");
		exit(1);
	}
	
    rho1 = atof(argv[1]);
    assert(rho1 > 0.0);

	u1 = atof(argv[2]);
    assert(u1 > 0.0);

    rho2 = atof(argv[3]);
    assert(rho2 > 0.0);

    iMat = atoi(argv[4]);
    assert(iMat >= 0);

	fprintf(stderr, "Initializing material...\n");
    
    tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);

#if 0
	fprintf(stderr, "Initializing the look up table...\n");
	granite = tillInitMaterial(iMat, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"\n");
#endif

	fprintf(stderr,"Starting from rho1=%g u1=%g\n", rho1, u1);

    drho = (rho1-rho2)/(nStep-1);
    u = u1;

    printf("#           rho              u\n");
    printf("%15.7E%15.7E\n", rho1, u);

    for (i=0; i<nStep-1; i++)
    {
        rho = rho1 - i*drho;

        u = tillCalcU(tillMat, rho, u, rho-drho);
        rho -= drho;

        printf("%15.7E%15.7E\n", rho, u);
    }

    u2 = tillCalcU(tillMat, rho1, u1, rho2);

    fprintf(stderr, "rho1= %g u1= %g rho2= %g u2= %g\n", rho1, u1, rho2, u2);

	tillFinalizeMaterial(tillMat);
}
