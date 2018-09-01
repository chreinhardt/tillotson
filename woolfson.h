/*
 * The header file for woolfson.c.
 */
#ifndef WOOLFSON_HINCLUDED
#define WOOLFSON_HINCLUDED

#include "tillotson.h"

// Minimum pressure to avoid problems in root finding if P(rho) is not monotonic
#define WOOLFSON_MIN_PRESSURE 1e-3

typedef struct woolfson_coeff_table_entry
{
//	double u;
//	double rho;
	double P;
	double T;
	double f;
} WOOLFSON_COEFF_TABLE_ENTRY;

typedef struct woolfson_coeff_table
{
	TILLMATERIAL **tillMat;		/* Array that contains pointers to the materials */
	int nMat;					/* How many materials are there? */
	int nTableP;				/* Number of entries in the look up table in P */
	int nTableT;				/* and in T */
	double Pmin;				/* Min value for the lookup table */
	double Pmax;				/* Max value for the lookup table */
	double Tmin;
	double Tmax;

	double dP;
	double dT;
	
	/* An array of lookup tables for f_ij */
	WOOLFSON_COEFF_TABLE_ENTRY **Lookup;
} WOOLFSON_COEFF_TABLE;

double CalcWoolfsonCoeff(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double P, double T);

// Basic functions for the lookup table
WOOLFSON_COEFF_TABLE_ENTRY **CoeffMatrixAlloc(int nRow, int nCol);
WOOLFSON_COEFF_TABLE* InitWoolfsonCoeffTable(TILLMATERIAL *Mat1, TILLMATERIAL *Mat2, int nP,
                                             int nT, double Pmin, double Pmax, double Tmin,
                                             double Tmax);

int WoolfsonLookupPIndex(WOOLFSON_COEFF_TABLE_ENTRY **Lookup, double P, unsigned int nRow,
                         unsigned int iCol);
int WoolfsonLookupTIndex(WOOLFSON_COEFF_TABLE_ENTRY **Lookup, double T, unsigned int iRow,
                         unsigned int nCol);
double WoolfsonCoeffInterpol(WOOLFSON_COEFF_TABLE *table, double P, double T);

int PrintWoolfsonCoeffTable(WOOLFSON_COEFF_TABLE *table);

#endif

