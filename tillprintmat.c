/*
 * Print material data for a given material.
 *
 * Author:   Christian Reinhardt
 * Date:     01.10.2018 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) ((i*tillMat->nTableV) + (j))

void main(int argc, char **argv)
{
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	double rho, u;
	int nTableRho = 1000;
	int nTableV = 1000;

    int iMat;

	TILLMATERIAL *tillMat;

    if (argc != 2)
    {
        fprintf(stderr,"Usage: tillprintmat <iMat>\n");
        exit(1);
	}

    iMat = atoi(argv[1]);

    assert(iMat >= 0);

	tillMat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);

    printf("\n");
    tillPrintMat(tillMat);
    printf("\n");
    
    tillFinalizeMaterial(tillMat);
}
