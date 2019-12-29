/*
 * Calculate the region where cs2<0 for the Tillotson EOS.
 *
 * Author:   Christian Reinhardt
 * Date:     29.12.2019
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"


int main(int argc, char **argv) {
    // Tillotson EOS library
	TILLMATERIAL *tillMat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rho;
    double a, b, c;
    double Pa, Pb, Pc;
    double cs2a, cs2b, cs2c;

	tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);	
	
	/*
	** Find where P<0 for 0 < rho < rho0.
	*/
	rho=TILL_RHO_MIN;

#ifdef TILL_PRESS_MELOSH
	/* In this case the pressure is set to zero for eta<0.8. */
	fprintf(stderr,"TILL_PRESS_MELOSH defined!\n");
	exit(1);
#endif
	
	while (rho <= tillMat->rho0)
	{
		/* Do bisection to find where P<0. */
		a = 1e-10;
		b = tillMat->us2;
		c = 0.0;

		Pa = tillPressureSoundNP(tillMat, rho, a, &cs2a);
		Pb = tillPressureSoundNP(tillMat, rho, b, &cs2b);
		cs2c = 0.0;

		while (b-a > 1e-14)
		{
//			printf("rho=%g a=%g b=%g c=%g Pa=%g Pb=%g Pc=%g\n",rho,a,b,c,Pa,Pb,Pc);
			c = 0.5*(a + b);
			Pc = tillPressureSoundNP(tillMat, rho, c, &cs2c);

			if (cs2c < 0)
			{
				// Set a = c
				a = c;
				cs2a = cs2c;
			} else {
				// Set b = c
				b = c;
				cs2b = cs2c;
			}
		}

		printf("%15.7E %15.7E\n", rho, c);
		rho += 0.01;
	}

	fprintf(stderr,"Done.\n");

	tillFinalizeMaterial(tillMat);

    return 0;
}
