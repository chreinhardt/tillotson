/*
 * The header file for tillwoolfson.c.
 */
#ifndef TILLWOOLFSON_HINCLUDED
#define TILLWOOLFSON_HINCLUDED

#include "tillotson.h"

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

//double tillInitInterfaceLookup(TILLMATERIAL *material1 TILLMATERIAL *material2, double rho, double u);
double CalcWoolfsonCoeff(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double P, double T);

#endif

