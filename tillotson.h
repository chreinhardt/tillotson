/*
 ** The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_HINCLUDED
#define TILLOTSON_HINCLUDED
/*
#include <sys/time.h>
*/

#define GRANITE 1
#define IRON 2

struct lookup
{
	double u;
	double rho;
	double d2udrho2; /* used for cubic spline interpolation */
};

typedef struct tillMaterial
{
	int iMaterial; /* What material is it? */
	int nTableMax; /* Max. number of entries in the look up table */
	int nTable; /* number of entries in the look up table */
	double rhomax; /* Max value for the lookup table */

	/* Unit convertion factors */
	double dKpcUnit;
	double dMsolUnit;
	
	double dGasConst;
	double dErgPerGmUnit;
	double dGmPerCcUnit;
	double dSecUnit;

	/* A material is defined by 10 Tillotson parameters */
	double a;
	double b;
	double A;
	double B;
	double rho0;
	double u0;
	double us;
	double us2;
	double alpha;
	double beta;

	 /* the specific heat capacity */
	double cv;

	/* The cold curve */	
	struct lookup *cold;
	double delta;
} TILLMATERIAL;

TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, double rhomax);
void tillFinalizeMaterial(TILLMATERIAL *material);
double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tillSoundSpeed2old(TILLMATERIAL *material, double rho, double u);
double tillPressureSoundold(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressure(TILLMATERIAL *material, double rho, double u);
double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tilldPdu(TILLMATERIAL *material, double rho, double u);
double tilldTdrho(TILLMATERIAL *material, double rho, double u);
double tilldTdu(TILLMATERIAL *material, double rho, double u);
double tillTempRhoU(TILLMATERIAL *material, double rho, double u);
double tillSoundSpeed(TILLMATERIAL *material, double rho, double u);
void tillInitColdCurve(TILLMATERIAL *material);
double tillColdULookup(TILLMATERIAL *material,double rho);
#endif

