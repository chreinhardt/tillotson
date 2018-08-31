/*
 * Copyright (c) 2018 Christian Reinhardt.
 *
 * This file provides all the functions to do the density correction at the material interface
 * proposed in Woolfson (2007).
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"
#include "woolfson.h"

/*
 * Functions:
 *
 * Initialize/Finalize:
 *
 * InitWoolfsonCoeffTable: Allocate memory and generate the lookup table.
 *
 * WoolfsonCoeffInterpol: Linear interpolate in the lookup table to find f_ij(P, T).
 *
 * CalcWoolfsonCoeff: Calculate the correction coefficients f_ij as in Woolfson (2007).
 */
double CalcWoolfsonCoeff(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double P, double T)
{
    /*
     * Calculate the coefficient
     *
     * f_ij := rho_mat1(P, T)/rho_mat2(P, T)
     *
     * needed to correct the density at a material interface.
     */
    double rho1;
    double rho2;

    // Check if there is indeed a material interface.
    if (mat1->iMaterial == mat2->iMaterial)
        return 1.0;

    /*
     * In the low density region a density correction can be problematic (dPdrho not monotonic,
     * interpretation of mixed phases in the Tillotson EOS difficult), so the correction factor
     * is one in this case.
     */
    if ((rho1 < mat1->rho0) || (rho2 < mat2->rho0))
        return 1.0;

    rho1 = eosRhoPTemp(mat1, P, T);
    rho2 = eosRhoPTemp(mat2, P, T);

//    fprintf(stderr, "rho1= %g rho2= %g\n", rho1, rho2);

    // If the density is unphysical return 1 so the density is not corrected.
    if ((rho1 <= 0.0) || (rho2 <= 0.0))
        return 1.0;

    return (rho1/rho2);
}

/*
 * Allocate memory for the look up table as described in the Numerical recipes.
 */
WOOLFSON_COEFF_TABLE_ENTRY **CoeffMatrixAlloc(int nRow, int nCol)
{
    WOOLFSON_COEFF_TABLE_ENTRY **Lookup;
    int i;


	Lookup = (WOOLFSON_COEFF_TABLE_ENTRY **) calloc(nRow, sizeof(WOOLFSON_COEFF_TABLE_ENTRY*));
	Lookup[0] = (WOOLFSON_COEFF_TABLE_ENTRY *) calloc(nRow*nCol, sizeof(WOOLFSON_COEFF_TABLE_ENTRY));

	assert(Lookup != NULL);

	/* Set a pointer to each row. */
	for (i=1; i<nRow; i++)
	{
		Lookup[i] = Lookup[i-1]+nCol;
	}

    return Lookup;
}

/*
 * This function does allocate memory and generate the lookup table for the coefficients.
 */
WOOLFSON_COEFF_TABLE* InitWoolfsonCoeffTable(TILLMATERIAL *Mat1, TILLMATERIAL *Mat2, int nP,
                                             int nT, double Pmin, double Pmax, double Tmin,
                                             double Tmax)
{
    WOOLFSON_COEFF_TABLE *table;
    double P, dP;
    double T, dT;
    int i, j;

    table = (WOOLFSON_COEFF_TABLE *) calloc(1, sizeof(WOOLFSON_COEFF_TABLE));
    assert(table != NULL);

    table->nTableP = nP;
    table->nTableT = nT;
    
    table->Pmin = Pmin;
    table->Pmax = Pmax;
    table->Tmin = Tmin;
    table->Tmax = Tmax;

    // Make sure the EOS is initialized.
    assert(Mat1 != NULL);
    if (Mat1->Lookup == NULL)
        tillInitLookup(Mat1);

    assert(Mat2 != NULL);
    if (Mat2->Lookup == NULL)
        tillInitLookup(Mat2);

    // Allocate nP x nT table
    table->Lookup = CoeffMatrixAlloc(nP, nT);

    // Generate a linear grid in P and T.
    dP = (table->Pmax-table->Pmin)/(table->nTableP-1);
    dT = (table->Tmax-table->Tmin)/(table->nTableT-1);

    for (i=0; i< table->nTableP; i++)
    {
        P = table->Pmin + dP*i;
        for (j=0; j<table->nTableT; j++)
        {
            T = table->Tmin + dT*j;

            table->Lookup[i][j].P = P;
            table->Lookup[i][j].T = T;

            table->Lookup[i][j].f = CalcWoolfsonCoeff(Mat1, Mat2, P, T);
        }
    }
    
    return table;
}

/*
 * This function frees all allocated memory.
 */
WOOLFSON_COEFF_TABLE* FinalizeWoolfsonCoeffTable(TILLMATERIAL *Mat1, TILLMATERIAL *Mat2, int nP,
                                                 int nT, double Pmin, double Pmax, double Tmin,
                                                 double Tmax)
{
    
}

/*
 * Do bisection to find the values P_i and P_i+1 that bracket P in one column of the the lookup
 * table assuming it is sorted. The function returns the index i.
 */
int WoolfsonLookupPIndex(WOOLFSON_COEFF_TABLE_ENTRY **Lookup, double P, unsigned int nRow, unsigned int iCol)
{
	unsigned int iLower, iUpper, i;
	
	iLower = 0;
	iUpper = nRow-1;

	/*
	 * Make sure that P is in the lookup table. If this is not the case return -1 or nRow.
     */
	if (P < Lookup[iLower][iCol].P)
	{
		return -1;
	} else if (P > Lookup[iUpper][iCol].P)
	{
		return nRow;
	}

	assert(Lookup[iLower][iCol].P < Lookup[iUpper][iCol].P);
	
	/* Do bisection. */
	while (iUpper-iLower > 1)
	{
		/* Compute the midpoint. */
		i = (iUpper + iLower) >> 1;

		if (P >= Lookup[i][iCol].P)
		{
			iLower = i;
		} else {
			iUpper = i;
		}
	}
	
	if (P == Lookup[0][iCol].P) return 0;

	/* 
	 * Return Pmax-1 as the lower index, so the desired value is always bracketed by (i,i+1).
	 */
	if (P == Lookup[nRow-1][iCol].P)
	{
			return (nRow-2);
	}

	return iLower;
}

/*
 * Do bisection to find the values T_j and T_j+1 that bracket u in one row of the the lookup table
 * assuming it is sorted. The function returns the index j.
 */
int WoolfsonLookupTIndex(WOOLFSON_COEFF_TABLE_ENTRY **Lookup, double T, unsigned int iRow, unsigned int nCol)
{
	unsigned int iLower,iUpper,i;
	
	iLower = 0;
	iUpper = nCol-1;

	/*
	 * Make sure that T is in the lookup table. If this is not the case return -1 or nCol.
     */
	if (T < Lookup[iRow][iLower].T)
	{
		return -1;
	} else if (T > Lookup[iRow][iUpper].T)
    {
		return nCol;
	}

	assert(Lookup[iRow][iLower].T < Lookup[iRow][iUpper].T);
	
	/* Do bisection. */
	while (iUpper-iLower > 1)
	{
		/* Compute the midpoint. */
		i = (iUpper + iLower) >> 1;

		if (T >= Lookup[iRow][i].T)
		{
			iLower = i;
		} else {
			iUpper = i;
		}
	}
	
	if (T == Lookup[iRow][0].T) return 0;
	/* 
	 * Return umax-1 as the lower index, so the desired value is always bracketed by (i,i+1).
	 */
	if (T == Lookup[iRow][nCol-1].T)
	{
			return (nCol-2);
	}

	printf("T= %g T(i,jmax)= %g iLower= %i T(%i, %i)= %g\n", T, Lookup[iRow][nCol-1].T, iLower, iLower, nCol-1, Lookup[iLower][nCol-1].T);
	return iLower;
}

/*
 * Determine the coefficient f_ij(P, T) from the look up table using linear interpolation.
 */
double WoolfsonCoeffInterpol(WOOLFSON_COEFF_TABLE *table, double P, double T)
{
    double x, y;
    int i, j;

    /*
     * Since P(i, j) = P(i) and T(i, j) = T(j) we can do two 1D bisection to determine i and j.
     */
    i = WoolfsonLookupPIndex(table->Lookup, P, table->nTableP, 0);
    printf("i=%i ",i);
    
    j = WoolfsonLookupTIndex(table->Lookup, T, 0, table->nTableT);
    printf("j=%i ",j);

    /*
     * Check if the data is in the lookup table.
     */
    if (i < 0 || i >= table->nTableP)
    {
        printf("P= %15.7E i= %i is outside of the table.\n", P, i);
        assert(0);
    }

    if (j < 0 || j >= table->nTableT)
    {
        printf("T= %15.7E j= %i is outside of the table.\n", T, j);
        assert(0);
    }

    /*
     * We assume, that the grid is evenly spaced in P and T.
     */
    x = (table->Lookup[i+1][j].P - P)/(table->Lookup[i+1][j].P - table->Lookup[i][j].P);
    y = (table->Lookup[i][j+1].T - T)/(table->Lookup[i][j+1].T - table->Lookup[i][j].T);

    return (x*y*table->Lookup[i][j].f + x*(1.0-y)*table->Lookup[i][j+1].f +
        (1.0-x)*y*table->Lookup[i+1][j].f + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].f);
}

int PrintWoolfsonCoeffTable(WOOLFSON_COEFF_TABLE *table)
{
    int i, j;

    if (table == NULL)
        return 0;

    for (i=0; i< table->nTableP; i++)
    {
        for (j=0; j<table->nTableT; j++)
        {
//            printf(" %15.7E", table->Lookup[i][j].f);
            printf(" %15.7E", table->Lookup[i][j].T);
        }
        printf("\n");
    }

    return 1;
}

