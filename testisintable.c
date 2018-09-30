/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     13.09.2018
 * Modified: 29.09.2018 
 *
 * Test tillIsInTable() by checking, if a grid of points is in the lookup table. Especially the
 * values close to v=0 and v=vmax can be problematic.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) (((i)*tillMat->nTableV) + (j))

void main(int argc, char **argv) {
    // Tillotson EOS library
    TILLMATERIAL *tillMat;
    double dKpcUnit = 2.06701e-13;
    double dMsolUnit = 4.80438e-08;
    double rhomin = TILL_RHO_MIN;
    double rhomax = 100.0;
    double vmax = 1200.0;
    int nTableRho = 100;
    int nTableV = 100;
    double rho, v, u;
    double umax;
    double rho2, u2;
    char ErrorString[256]; 
    int iRet;
    FILE *fp = NULL;
    int i = 0;
    int j = 0;

#ifdef TILL_PRESS_NP
    fprintf(stderr, "TILL_PRESS_NP.\n");
#endif
    fprintf(stderr, "Initializing material...\n");

    tillMat = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit);

    fprintf(stderr, "Initializing the look up table...\n");

    /* Solve ODE and splines */
    tillInitLookup(tillMat, nTableRho, nTableV, rhomin, rhomax, vmax);
    fprintf(stderr, "Done.\n");

    fprintf(stderr,"\n");
    fprintf(stderr,"rhomax: %g, vmax: %g \n", tillMat->rhomax, tillMat->vmax);
    fprintf(stderr,"nTableRho: %i, nTableV: %i \n", tillMat->nTableRho, tillMat->nTableV);
    fprintf(stderr,"drho: %g, dv: %g \n", tillMat->drho, tillMat->dv);
    fprintf(stderr,"\n");

    rho = 0.0;
    u = 0.0;


    /*
     * Print the look up table to a file first.
     */	
    fp = fopen("lookup.txt", "w");
    assert(fp != NULL);

    for (i=0; i<tillMat->nTableRho; i+=1)
    {
        rho = tillLookupRho(tillMat, i);
        fprintf(fp, "%15.7E",rho);

        for (j=0; j<tillMat->nTableV; j+=1)
        {
            u = tillMat->Lookup[INDEX(i, j)].u;
            fprintf(fp,"  %15.7E", u);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    /* 
     * Now check if points on the cold curve and the last isentrope are treated correctly.
     *
     * Here we also use points that are on the grid to check, if rounding errors in
     * tillLookupIndexLogRho() affect the test.
     */
    fp = fopen("testisintable2.txt", "w");
    assert(fp != NULL);
#if 0
    /* First do the cold curve. */
    v = 0.0;

    while (v < 1e-2)
    {
        fprintf(stderr, "v= %15.7E\n", v);
        for (i=0; i<tillMat->nTableRho-1; i+=1)
        {
            for (j=0; j<100; j++)
            {
                rho = tillMat->rhomin*exp((i+0.01*j)*tillMat->dlogrho);
                u = tillCubicInt(tillMat, rho, v);
/*
                fprintf(stderr, "\n");
                fprintf(stderr,"i= %i: Testing rho=%g, u=%g (v= %g)! Rho: index= %g = %i\n", i, rho, u, v, (log(rho)-log(tillMat->rhomin))/tillMat->dlogrho, (int) floor((log(rho)-log(tillMat->rhomin))/tillMat->dlogrho));
*/
                iRet = tillIsInTable(tillMat, rho, u);

                if (iRet != TILL_LOOKUP_SUCCESS)
                {
                    tillErrorString(iRet, ErrorString);
                    fprintf(stderr,"i= %i: rho=%15.7E, u=%15.7E (v= %15.7E) not in table (Error %s)!\n", i, rho, u, v, ErrorString);
#if 0
                    fprintf(stderr, "Calling tillLookupU.\n");
                    fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
                    u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
                    fprintf(stderr, "u2= %g\n", u);
#endif
                    //			fprintf(fp, "%15.7E %15.7E\n", rho, u);
                } else {
                    /* Test if tillLookupU() really works. */
                    rho2 = tillMat->rhomin*exp((i + 0.02*j)*tillMat->dlogrho);
                    fprintf(stderr, "i= %i: rho=%15.7E, u=%15.7E (v= %15.7E) is in table.\n", i, rho, u, v);
                    fprintf(stderr, "Calling tillLookupU.\n");
                    fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, rho2);
                    u = tillLookupU(tillMat, rho, u, rho2, 0);
                    fprintf(stderr, "u2= %g\n", u);
//                    assert(0);
                }
                
                fprintf(fp, "%15.7E%15.7E%2i\n", rho, u, iRet);
            }
        }
        v += 1e-2/(100-1);
    }
#endif
    /* Now vmax (note that tillCubicInt fails if v==vmax). */
    v = tillMat->vmax-1e-9;
    //v = tillMat->vmax-0.99999*tillMat->dv;

#if 0
    fprintf(stderr, "Checking v==vmax (v= %15.7E, vmax= %15.7E, dv= %15.7E).\n", v, tillMat->vmax,
            tillMat->dv);
    fprintf(stderr, "v_N-1= %15.7E v_N-2 %15.7E dv= %g\n", tillLookupV(tillMat, tillMat->nTableV-1), 
            tillLookupV(tillMat, tillMat->nTableV-2), tillLookupV(tillMat, tillMat->nTableV-1)-tillLookupV(tillMat, tillMat->nTableV-2));
    fprintf(stderr, "delta_v = %15.7E (dv= %g)\n", tillMat->vmax-v, 0.99999*tillMat->dv);
#endif
    while (v > tillMat->vmax-tillMat->dv)
    {
        fprintf(stderr, "v= %15.7E\n", v);
        for (i=0; i<tillMat->nTableRho-1; i+=1)
        {
            for (j=0; j<100; j++)
            {
                rho = tillMat->rhomin*exp((i+0.01*j)*tillMat->dlogrho);
                u = tillCubicInt(tillMat, rho, v);
/*
                fprintf(stderr, "\n");
                fprintf(stderr,"i= %i: Testing rho=%g, u=%g (v= %g)! Rho: index= %g = %i\n", i, rho, u, v, (log(rho)-log(tillMat->rhomin))/tillMat->dlogrho, (int) floor((log(rho)-log(tillMat->rhomin))/tillMat->dlogrho));
*/
                iRet = tillIsInTable(tillMat, rho, u);

                if (iRet != TILL_LOOKUP_SUCCESS)
                {
                    tillErrorString(iRet, ErrorString);
                    fprintf(stderr,"i= %i: rho=%15.7E, u=%15.7E (v= %15.7E) not in table (Error %s)!\n", i, rho, u, v, ErrorString);
#if 0
                    fprintf(stderr, "Calling tillLookupU.\n");
                    fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
                    u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
                    fprintf(stderr, "u2= %g\n", u);
#endif
                    //			fprintf(fp, "%15.7E %15.7E\n", rho, u);
                } else {
                    /* Test if tillLookupU() really works. */
//#if 0
                    rho2 = tillMat->rhomin*exp((i + 0.02*j)*tillMat->dlogrho);
                    fprintf(stderr, "i= %i: rho=%15.7E, u=%15.7E (v= %15.7E) is in table.\n", i, rho, u, v);
                    fprintf(stderr, "Calling tillLookupU.\n");
                    fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, rho2);
                    u2 = tillLookupU(tillMat, rho, u, rho2, 0);
                    fprintf(stderr, "u2= %g\n", u2);
//#endif
                }
                
                fprintf(fp, "%15.7E%15.7E%2i\n", rho, u, iRet);
            }
        }
        v  -= tillMat->dv/(100-1);
    }


    fclose(fp);

#if 0
    v = tillLookupV(tillMat, tillMat->nTableV-1);
    //    v -= 1e-8;
    //    v -= 2.5*tillMat->dv;
    j = tillLookupIndexV(tillMat, v);

    for (i=0; i<tillMat->nTableRho; i+=1)
    {
        rho = tillLookupRho(tillMat, i);
        u = tillMat->Lookup[INDEX(i, j)].u;

        fprintf(stderr, "\n");
        fprintf(stderr,"i= %i j= %i: Testing rho=%g, u=%g (v= %g)! index= %g = %i\n", i, j, rho, u, v, (log(rho)-log(tillMat->rhomin))/tillMat->dlogrho, (int) floor((log(rho)-log(tillMat->rhomin))/tillMat->dlogrho));


        iRet = tillIsInTable(tillMat, rho, u);

        if (iRet != TILL_LOOKUP_SUCCESS)
        {
            tillErrorString(iRet, ErrorString);
            fprintf(stderr,"i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) not in table (Error %s)!\n", i, j, rho, u, v, ErrorString);

            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);

            //			fprintf(fp, "%15.7E %15.7E\n", rho, u);
        } else {
            fprintf(stderr, "i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) is in table.\n", i, j, rho, u, v);
            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);
            assert(0);
        }
    }
#endif

#if 0
    v = tillLookupV(tillMat, tillMat->nTableV-1);
    v -= 1e-10;
    //    v -= 2.5*tillMat->dv;
    j = tillLookupIndexV(tillMat, v);

    for (i=0; i<tillMat->nTableRho-1; i+=1)
    {
        // Choose a point between the grid points (logarithmic spacing)
        rho = tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho);

        // Note that this requires v < vmax.
        u = tillCubicInt(tillMat, rho, v);

        fprintf(stderr, "\n");
        fprintf(stderr,"i= %i j= %i: Testing rho=%g, u=%g (v= %g)! index= %g = %i\n", i, j, rho, u, v, (log(rho)-log(tillMat->rhomin))/tillMat->dlogrho, (int) floor((log(rho)-log(tillMat->rhomin))/tillMat->dlogrho));


        iRet = tillIsInTable(tillMat, rho, u);

        if (iRet != TILL_LOOKUP_SUCCESS)
        {
            tillErrorString(iRet, ErrorString);
            fprintf(stderr,"i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) not in table (Error %s)!\n", i, j, rho, u, v, ErrorString);

            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);

            //			fprintf(fp, "%15.7E %15.7E\n", rho, u);
        } else {
            fprintf(stderr, "i= %i j= %i: rho=%15.7E, u=%15.7E (v= %15.7E) is in table.\n", i, j, rho, u, v);
            fprintf(stderr, "Calling tillLookupU.\n");
            fprintf(stderr, "rho1= %g u1= %g rho2= %g ", rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho));
            u = tillLookupU(tillMat, rho, u, tillMat->rhomin*exp((i + 0.5)*tillMat->dlogrho), 0);
            fprintf(stderr, "u2= %g\n", u);
            assert(0);
        }
    }
#endif
    fprintf(stderr,"Done.\n");
    tillFinalizeMaterial(tillMat);
}
