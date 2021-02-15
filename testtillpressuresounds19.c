/*
 * Calculate the pressure and sound speed in the cold and intermediate expanded states for
 * different assumptions.
 *
 * Author:  Christian Reinhardt
 * Date:    14.02.2021
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

//#define TILL_PRESS_NP

/*
 * Calculate the pressure and sound speed from the Tillotson EOS as described in Stewart (2019).
 * The pressure is set to 1e-10 if P < 0 but negative pressures are included in the intermediate
 * expanded states.
 */
double tillPressureSoundS19(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
    double Pmin = 1e-10;
    double eta, mu;
    double Pc, Pe, Pint;
    double c2c, c2e;
    double Gammac, Gammae, w0, y, z;

    eta = rho/material->rho0;
    mu = eta - 1.0;
    z = (1.0 - eta)/eta;
    w0 = u/(material->u0*eta*eta)+1.0;

    /*
     *  Here we evaluate, which part of the equation of state we need.
     */
    if ((rho >= material->rho0) || (u < material->us)) {
        /*
         *  Condensed states (rho > rho0) or expanded cold states.
         */
        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;

        /* Set Pc to the minimum pressure. */
        if (Pc < Pmin)
            Pc = Pmin;

        if (pcSound != NULL)
        {
            /* Calculate the sound speed. */
            c2c =material->a*u+material->b*u/(w0*w0)*(3.0*w0-2.0)+(material->A+2.0*material->B*mu)/material->rho0 + Pc/(rho*rho)*(material->a*rho+material->b*rho/(w0*w0));
            *pcSound = sqrt(c2c);
        }
        return Pc;
    } else if (u > material->us2) {
        /*
         * Expanded hot states (rho < rho0 and u > us2).
         */
        Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
        Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

        if (pcSound != NULL)
        {
            /* calculate the sound speed */
            c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));

            *pcSound = sqrt(c2e);
        }

        return Pe;
    } else {
        /*
         *  intermediate states (rho < rho0 and us < u < us2)
         */
        y = (u - material->us)/(material->us2 - material->us);

        Gammac = material->a + material->b/w0;
        Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
        Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
        Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

        Pint = (Pc*(1.0-y)+Pe*y);

        /* Set Pint to the minimum pressure. */
        if (Pint < Pmin)
            Pint = Pmin;
        
        if (pcSound != NULL)
        {
            /* calculate the sound speed */
            c2c =material->a*u+material->b*u/(w0*w0)*(3.0*w0-2.0)+(material->A+2.0*material->B*mu)/material->rho0 + Pc/(rho*rho)*(material->a*rho+material->b*rho/(w0*w0));
            c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));

            *pcSound = sqrt(c2c*(1.0-y)+c2e*y);
        }


        return Pint;
    }
}

int main(int argc, char **argv) {
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
    int iMat = GRANITE;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    // Zoom into the expanded states
    double rhomin = 1e-4;
	double rhomax = 10.0;
    double umin = 0.0;
    double umax = 25.0;
    int nRho = 1000;
    int nU = 100;
    double drho;
    double du;
	double rho, u;
    double P, PS19;
    double cs, csS19;
    FILE *fp1, *fp2;
    int i, j;

	fprintf(stderr, "Initializing material...\n");
	tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
	fprintf(stderr, "Done.\n");
	fprintf(stderr, "\n");

    /* Print the rho and u axis. */
	drho = (rhomax-rhomin)/(nRho-1);
    du = (umax-umin)/(nU-1);

    fp1 = fopen("testtillpressuresounds19_rhoaxis.txt", "w");
    fp2 = fopen("testtillpressuresounds19_uaxis.txt", "w");

    assert(fp1 != NULL);
    assert(fp2 != NULL);

    for (i=0; i<nRho; i++) {
        rho = rhomin + i*drho;
        fprintf(fp1, "%15.7E\n", rho);
    }
        
    for (j=0; j<nU; j++) {
        u = umin + j*du;
        fprintf(fp2, "%15.7E\n", u);
    }
    fclose(fp1);
    fclose(fp2);

	/* Calculate P(rho, u) and cs(rho, u) at the grid points and print the relative error. */	
    fp1 = fopen("testtillpressures19.txt", "w");
    fp2 = fopen("testtillsounds19.txt", "w");

    assert(fp1 != NULL);
    assert(fp2 != NULL);

    for (j=0; j<nU; j++) {
        u = umin + j*du;
        for (i=0; i<nRho; i++) {
            rho = rhomin + i*drho;
            P = tillPressureSound(tillMat, rho, u, &cs);
            PS19 = tillPressureSoundS19(tillMat, rho, u, &csS19);
            fprintf(fp1, "  %15.7E", (P-PS19)/P);
            fprintf(fp2, "  %15.7E", (cs-csS19)/cs);
        }
        fprintf(fp1, "\n");
        fprintf(fp2, "\n");
    }

    fclose(fp1);
    fclose(fp2);

	tillFinalizeMaterial(tillMat);

    return 0;
}
