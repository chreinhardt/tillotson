/*
 * The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_HINCLUDED
#define TILLOTSON_HINCLUDED

//#include "tillwoolfson.h"

/*
 * Version.
 */
#define TILL_VERSION_TEXT    "3.4.1"
#define TILL_VERSION_MAJOR   3
#define TILL_VERSION_MINOR   4
#define TILL_VERSION_PATCH   1

/*
 * Error codes.
 */
#define TILL_SUCCESS 0
#define TILL_FAIL   -1

enum till_error_lookup
{
    TILL_LOOKUP_SUCCESS = 0,
    TILL_LOOKUP_OUTSIDE_RHOMIN,
    TILL_LOOKUP_OUTSIDE_RHOMAX,
    TILL_LOOKUP_OUTSIDE_VMIN,
    TILL_LOOKUP_OUTSIDE_VMAX
};

/*
 * Phase.
 */
#define TILL_PHASE_FAIL      -1
#define TILL_PHASE_CONDENSED  1
#define TILL_PHASE_EXP_COLD   2
#define TILL_PHASE_EXP_HOT    3
#define TILL_PHASE_EXP_INT    4

/*
 * Warning: redefined the material constants on 10.10.2017!
 */
#define IDEALGAS    0
#define GRANITE     1
#define IRON        2
#define BASALT      3
#define ICE         4
#define WATER       5
#define DUNITE      6

/* 
 * ANEOS and M-ANEOS.
 */
#define ANEOS_IRON   100
#define ANEOS_DUNITE 101



#define TILL_N_MATERIAL_MAX 102

/* Compile the code with a pressure cutoff at negative pressures in the expanded cold states. */
#ifndef TILL_PRESS_NP
#define TILL_PRESS_NP
#endif

/* We write the look up table as a 1D array where a(i,j)=a(i*N+j) */
#define TILL_INDEX(i, j) (((i)*material->nTableV) + (j))

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

/* Define a minimum density for the look up table */
#ifndef TILL_RHO_MIN
#define TILL_RHO_MIN 1e-4
#endif

/* Define eps so that v <= v_max-eps. */
#define V_EPS 1e-8

/* Define FALSE and TRUE. */
//const int FALSE = 0;
//const int TRUE = 1;

typedef struct till_lookup_entry
{
	double u;
	double rho;
    double logrho;
	double v;
    double u1;		// du/drho
    /*
     * The following 2 variables are the second derivatives of the above variables (u and u1) with
     * respect to v (which is the value of the constant entropy curve (adiabat) at rho_0). Both of
     * these are obtained by fitting splines to u and u1 runs in the v axis.
     */
    double udv2; 	// d2u/dv2
    double u1dv2;	// d2/dv2(du/drho)
#ifdef TILL_DEBUG_SPLINT
	// Only for debugging
	double udrho2;  // d2u/drho2
#endif
} TILL_LOOKUP_ENTRY;

typedef TILL_LOOKUP_ENTRY* TILL_LOOKUP;

typedef struct tillMaterial
{
	int iMaterial;	/* What material is it? */
//	int nTableMax;	/* Max. number of entries in the look up table */
	int nTableRho;	/* Number of entries in the look up table in rho */
	int nTableV;	/* and in v */
	int n;			/* number of steps from rho to zero */
	double rhomin;	/* Min value for the lookup table */
	double rhomax;	/* Max value for the lookup table */
	double vmax;	/* Max value for the lookup table */
	/* Unit convertion factors */
	double dKpcUnit;
	double dMsolUnit;
	
	double dGasConst;
	double dErgPerGmUnit;
	double dGmPerCcUnit;
	double dSecUnit;

	/*
     * A material is defined by 10 Tillotson parameters.
     */
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

    /*
     * The constants needed for the ideal gas EOS.
     */
    double dConstGamma;
    double dMeanMolMass;        // in units of m_hydrogen

	 /* the specific heat capacity (assumed constant) */
	double cv;

	/* The cold curve */	
	TILL_LOOKUP_ENTRY *cold;
	double drho;
	double dlogrho;
	double dv;

	/* A look up table for u(rho) along an isentrope */
//	TILL_LOOKUP *Lookup;	// this is an array of pointers
	TILL_LOOKUP_ENTRY *Lookup;	// this is an array of pointers
} TILLMATERIAL;

/* Data structure needed for the GSL root finder. */
struct PressureRhoT_GSL_Params {
    TILLMATERIAL *material;
    double P;
    double T;
};

/*
 * Tillotson.c
 */
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit);
void tillFinalizeMaterial(TILLMATERIAL *material);

void tilliMatString(TILLMATERIAL *material, char *MatName);
void tillErrorString(int iError, char *ErrorString);
void tillPrintMat(TILLMATERIAL *material, FILE *fp);

// Some functions that provide a general interface for EOS calls where other EOS, e.g., an ideal gas EOS can be implemented
double eosPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound);
double eosPressureSoundRhoT(TILLMATERIAL *material, double rho, double T, double *pcSound);
double eosPressure(TILLMATERIAL *material, double rho, double u);
double eosPressureRhoT(TILLMATERIAL *material, double rho, double T);
double eosdPdrho(TILLMATERIAL *material, double rho, double u);
double eosdPdu(TILLMATERIAL *material, double rho, double u);
double eosTempRhoU(TILLMATERIAL *material, double rho, double u);
double eosRhoPU(TILLMATERIAL *material, double P, double u);
double eosURhoP(TILLMATERIAL *material, double rho, double P);
double eosURhoTemp(TILLMATERIAL *material, double rho, double T);
double eosRhoPTemp(TILLMATERIAL *material, double P, double T);
int eosSolveBC(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double rho1, double u1, double *prho2, double *pu2);
double eosPhi(TILLMATERIAL *material, double rho, double u);
double eosGamma(TILLMATERIAL *material, double rho, double u);

double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tillPressureSoundNP(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressure(TILLMATERIAL *material, double rho, double u);
double tillPressureNP(TILLMATERIAL *material, double rho, double u);
double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tilldPdu(TILLMATERIAL *material, double rho, double u);
double tilldTdrho(TILLMATERIAL *material, double rho, double u);
double tilldTdu(TILLMATERIAL *material, double rho, double u);
double tillTempRhoU(TILLMATERIAL *material, double rho, double u);

// Not implemented yet
double tillTempRhoP(TILLMATERIAL *material, double rho, double P);
double tillURhoTemp(TILLMATERIAL *material, double rho, double T);
double PressureRhoT_GSL(double rho, void *params);
double tillRhoPTemp(TILLMATERIAL *material, double P, double T);
double tillSoundSpeed(TILLMATERIAL *material, double rho, double u);
double tillRhoPU(TILLMATERIAL *material, double P, double u);
double tilldudrho(TILLMATERIAL *material, double rho, double u);
double tilldudlogrho(TILLMATERIAL *material, double logrho, double u);
int tillSolveBC(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double rho1, double u1, double *prho2, double *pu2);

/*
 * tillinitlookup.c
 */
void tillInitColdCurve(TILLMATERIAL *material);
void tillInitLookup(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmax);
TILL_LOOKUP_ENTRY *tillSolveIsentrope(TILLMATERIAL *material, double v);
TILL_LOOKUP_ENTRY *tillSolveIsentropeLogRho(TILLMATERIAL *material, double v);
/* Use bsstep.c from the Numerical Recipes */
TILL_LOOKUP_ENTRY *tillSolveIsentropeBS(TILLMATERIAL *material, double v);
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2);
int tillIsInTable(TILLMATERIAL *material,double rho,double u);
int eosIsBelowColdCurve(TILLMATERIAL *material,double rho,double u);
int tillIsBelowColdCurve(TILLMATERIAL *material,double rho,double u);

/*
 * tillsplint.c
 */
void tillInitSplines(TILLMATERIAL *material);
void tillInitSplineRho(TILLMATERIAL *material);
void tillInitSplineV(TILLMATERIAL *material);
void tillInitSplineU(TILLMATERIAL *material);
void tillInitSplineU1(TILLMATERIAL *material);
// Just for debugging
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv);
double tillSplineIntv(TILLMATERIAL *material, double v, int irho);
// These two are needed for the interpolator
double tillSplineIntU(TILLMATERIAL *material, double v, int irho);
double tillSplineIntU1(TILLMATERIAL *material, double v, int irho);

double tillCubicInt(TILLMATERIAL *material, double rho, double v);
void cubicint(double u[2],double dudrho[2], double dudv[2], double dudvdrho[2], double rho[2], double rhoint, double *intvalues);

double tillFindUonIsentrope(TILLMATERIAL *material,double v,double rho);
/* Used for the root finder */
double denergy(TILLMATERIAL *material,double v,double rho,double u);
double tillFindEntropyCurve(TILLMATERIAL *material,double rho,double u,int iOrder);
double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);
double tillCubicIntRho(TILLMATERIAL *material, double rhoint, int iv);
double tillColdULookup(TILLMATERIAL *material,double rho);

int tillLookupIndexRho(TILLMATERIAL *material, double rho);
int tillLookupIndexLogRho(TILLMATERIAL *material, double logrho);
int tillLookupIndexV(TILLMATERIAL *material, double v);

double tillLookupRho(TILLMATERIAL *material, int i);
double tillLookupLogRho(TILLMATERIAL *material, int i);
double tillLookupV(TILLMATERIAL *material, int j);

// A general version of tillLookupU() that can be used as an interface for different EOS
double eosLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);

#endif

