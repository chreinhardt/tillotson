/*
 * Plot the pressure on a rho x u and a rho x T grid.
 *
 * Author:   Christian Reinhardt
 * Date:     03.03.2019
 * Modified: 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

int main(int argc, char **argv)
{
    // Tillotson EOS library
	TILLMATERIAL *tillmat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rho_max = 100.0;
	double v_max = 100.0;
	int nRho = 1000;
	int nV = 1000;
    int nU = 1000;
    int iMat = GRANITE;
	double rho, u;
    double rho_min = 1e-4;
    double u_min = 0.0;
    double u_max = 100.0;
    double drho;
    double du;
    FILE *fp;
	int i, j;

	fprintf(stderr, "Initializing material...\n");
	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
	fprintf(stderr, "Done.\n");

	//tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);

	/*
	 * Print P on a rho x u grid.
	 */	
	fp = fopen("plottillpressure.txt", "w");
	assert(fp != NULL);

	fprintf(stderr, "Printing grid using tillPressure()...\n");

    drho = (rho_max-rho_min)/(nRho-1);
    du = (u_max-u_min)/(nU-1);

	/* Print a rho x u grid. */
	for (i=0; i<nRho; i++)
	{
		rho = rho_min + i*drho;

		for (j=0; j<nU; j++)
		{
			u = u_min + j*du;
			fprintf(fp," %15.7E", tillPressure(tillmat, rho, u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

#if !defined(TILL_PRESS_NP) || !defined(TILL_PRESS_MELOSH)
    /*
     * Print where the pressure is negative.
     */
    fp = fopen("till_granite_np.txt","w");
	assert(fp != NULL);

	fprintf(stderr, "Printing where P = 0.0 ...\n");

    rho = 0.0; 
	while (rho < tillmat->rho0)
	{
        u = 0.0;
        u_max = u;
        while (u < tillmat->us2)
        {
            if (tillPressure(tillmat, rho, u) < 0.0)
            {
			    u_max = u;
            }
            u += 0.01;
		}

		fprintf(fp,"%15.7E %15.7E\n", rho,u_max);
        rho += 0.01;
	}
#if 0
    rho = 0.0; 
	while (rho < tillmat->rho0)
	{
		fprintf(fp,"%15.7E", rho);
        u = 0.0;
        while (u < tillmat->us2)
        {
            if (tillPressure(tillmat, rho, u) < 0.0)
            {
			    fprintf(fp," %15.7E", u);
            }
            u += 0.01;
		}

		fprintf(fp,"\n");
        rho += 0.01;
	}
#endif
    fclose(fp);
#endif


	fprintf(stderr, "Printing, where the two functions differ...\n");
    printf("# rho          u\n");

	tillFinalizeMaterial(tillmat);

    return 0;
}
