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

void main(int argc, char **argv) {
	double dKpcUnit;
	double dMsolUnit;
	double rho, u, P;
	double rhomax;
	
	TILLMATERIAL *granite;

/*
	if (argc != 5) {
	fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <Tcore> <gamma>  >myball.std\n");
	exit(1);
	}
    nDesired = atoi(argv[1]);
    mTot = atof(argv[2]);
    Tcore = atof(argv[3]);
    gamma = atof(argv[4]);
*/
    //gamma = 2.0; /* Choose a gamma between 0.2 and 0.6 */

}
