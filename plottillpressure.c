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
    int nT = 1000;
    int iMat = GRANITE;
    double rho, u, T;
    double rho_min = 1e-4;
    double u_min = 0.0;
    double u_max = 100.0;
    double T_min = 0.0;
    double T_max = 1e4;
    double drho;
    double du;
    double dT;
    FILE *fp;
    int i, j;

    fprintf(stderr, "Initializing material...\n");
    tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Done.\n");

    tillInitLookup(tillmat, nRho, nV, rho_min, rho_max, v_max);

    /*
     * Zoom in to the expanded states.
     */
    rho_max = 1.1*tillmat->rho0;
    u_max = 1.1*tillmat->us2;

    /*
     * Print P on a rho x u grid.
     */	
    fp = fopen("plottillpressure.txt", "w");
    assert(fp != NULL);

    fprintf(stderr, "Printing grid using tillPressure()...\n");

    drho = (rho_max-rho_min)/(nRho-1);
    du = (u_max-u_min)/(nU-1);

    for (j=0; j<nU; j++)
    {
        u = u_min + j*du;

        for (i=0; i<nRho; i++)
        {
            rho = rho_min + i*drho;
            fprintf(fp," %15.7E", tillPressure(tillmat, rho, u));
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    /* Print where the pressure is negative. */
    fp = fopen("plottillpressure2.txt", "w");
    assert(fp != NULL);

    fprintf(stderr, "Printing where P = 0.0 ...\n");

    for (j=0; j<nU; j++)
    {
        u = u_min + j*du;

        for (i=0; i<nRho; i++)
        {
            rho = rho_min + i*drho;
            if (tillPressure(tillmat, rho, u) <= 0.0) {
                fprintf(fp, "%2i", 0);
            } else {
                fprintf(fp, "%2i", 1);
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    /*
     * Print P on a rho x T grid.
     */	
    fp = fopen("plottillpressure-rhot.txt", "w");
    assert(fp != NULL);

    fprintf(stderr, "Printing grid using tillPressure()...\n");

    drho = (rho_max-rho_min)/(nRho-1);
    dT = (T_max-T_min)/(nT-1);

    for (j=0; j<nT; j++)
    {
        T = T_min + j*dT;

        for (i=0; i<nRho; i++)
        {
            rho = rho_min + i*drho;
            fprintf(fp," %15.7E", eosPressureRhoT(tillmat, rho, T));
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    /* Print where the pressure is negative. */
    fp = fopen("plottillpressure2-rhot.txt", "w");
    assert(fp != NULL);

    fprintf(stderr, "Printing where P = 0.0 ...\n");

    for (j=0; j<nT; j++)
    {
        T = T_min + j*dT;

        for (i=0; i<nRho; i++)
        {
            rho = rho_min + i*drho;
            if (eosPressureRhoT(tillmat, rho, T) <= 0.0) {
                fprintf(fp, "%2i", 0);
            } else {
                fprintf(fp, "%2i", 1);
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);


#if 0
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
#endif
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

//    fprintf(stderr, "Printing, where the two functions differ...\n");
    
    printf("# rho            u\n");

    fp = fopen("plottillpressure-axis.txt", "w");
    assert(fp != NULL);
    
    fprintf(fp, "# rho            u              T\n");

    /* Print a rho x u grid. */
    for (i=0; i<nRho; i++)
    {
        rho = rho_min + i*drho;
        u = u_min + i*du;
        T = T_min + i*dT;

        fprintf(fp, "%15.7E%15.7E%15.7E\n", rho, u, T);
    }

    fclose(fp);

    tillFinalizeMaterial(tillmat);

    return 0;
}
