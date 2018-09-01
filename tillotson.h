/*
 ** The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_HINCLUDED
#define TILLOTSON_HINCLUDED

//#include "tillwoolfson.h"

#include "interpol/coeff.h"
#include "interpol/interpol.h"
#include "nr/nrcubicspline.h"

/*
 * Version.
 */
#define TILL_VERSION_TEXT    "2.1.0"
#define TILL_VERSION_MAJOR   2
#define TILL_VERSION_MINOR   1
#define TILL_VERSION_PATCH   0


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

#define ANEOS_IRON   31
#define ANEOS_DUNITE 32


#define TILL_N_MATERIAL_MAX 32

/* We write the look up table as a 1D array where a(i,j)=a(i*N+j) */
#define TILL_INDEX(i, j) (((i)*material->nTableV) + (j))

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

/* Define a minimum density for the look up table */
//#define TILL_RHO_MIN 1e-2
//#define TILL_RHO_MIN 5e-1
#define TILL_RHO_MIN 0.0

/* Define FALSE and TRUE. */
//const int FALSE = 0;
//const int TRUE = 1;

typedef struct till_lookup_entry
{
	double u;
	double rho;
	double v;
    double u1;		// du/drho
    /*
    ** The following 2 variables are the second derivatives of the above
    ** variables (u and u1) with respect to v (which is the value of the
    ** constant entropy curve (adiabat) at rho_0). Both of these are obtained
    ** by fitting splines to u and u1 runs in the v axis.
    */
    double udv2; 	// d2u/dv2
    double u1dv2;	// d2/dv2(du/drho)
#ifdef TILL_DEBUG_SPLINT
	// Only for debugging
	double udrho2;
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
	double iExpV;	/* Set to 1 for uniform steps in v */
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
//	double delta;
	double drho;
	double dv;

	/* A look up table for u(rho) along an isentrope */
//	TILL_LOOKUP *Lookup;	// this is an array of pointers
	TILL_LOOKUP_ENTRY *Lookup;	// this is an array of pointers
} TILLMATERIAL;

//TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, int nTableMax, double rhomax, double vmax);
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, int nTableRho, int nTableV, double rhomax, double vmax, int iExpV);
void tillFinalizeMaterial(TILLMATERIAL *material);

// Some functions that provide a general interface for EOS calls where other EOS, e.g., an ideal gas EOS can be implemented
double eosPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound);
double eosPressure(TILLMATERIAL *material, double rho, double u);
double eosdPdrho(TILLMATERIAL *material, double rho, double u);
double eosdPdu(TILLMATERIAL *material, double rho, double u);
double eosTempRhoU(TILLMATERIAL *material, double rho, double u);
double eosRhoPU(TILLMATERIAL *material, double P, double u);
double eosURhoP(TILLMATERIAL *material, double rho, double P);
double eosRhoPTemp(TILLMATERIAL *material, double P, double T);
double eosPhi(TILLMATERIAL *material, double rho, double u);
double eosGamma(TILLMATERIAL *material, double rho, double u);

double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tillSoundSpeed2old(TILLMATERIAL *material, double rho, double u);
double tillPressureSoundold(TILLMATERIAL *material, double rho, double u, double *pcSound);
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
double tillRhoPTemp(TILLMATERIAL *material, double P, double T);
double tillSoundSpeed(TILLMATERIAL *material, double rho, double u);
double tillRhoPU(TILLMATERIAL *material, double P, double u);
double tilldudrho(TILLMATERIAL *material, double rho, double u);
int tillSolveBC(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double rho1, double u1, double *prho2, double *pu2);

// Moved to tillinitlookup.h
void tillInitColdCurve(TILLMATERIAL *material);
void tillInitLookup(TILLMATERIAL *material);
TILL_LOOKUP_ENTRY *tillSolveIsentrope(TILLMATERIAL *material, double v);
/* Use bsstep.c from the Numerical Recipes */
TILL_LOOKUP_ENTRY *tillSolveIsentropeBS(TILLMATERIAL *material, double v);
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2);
int tillIsInTable(TILLMATERIAL *material,double rho,double u);
int tillIsBelowColdCurve(TILLMATERIAL *material,double rho,double u);

// Moved to tillsplint.h
/* Stuff for the cubic spline interpolator */
void tillInitSplines(TILLMATERIAL *material);
void tillInitSplineRho(TILLMATERIAL *material);
void tillInitSplinev(TILLMATERIAL *material);
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

// A general version of tillLookupU() that can be used as an interface for different EOS
double eosLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);

void tillBSderivs(TILLMATERIAL *material, float x, float y[], float dydx[]);




/* Defines for the Numerical Recipes routines */

/*
** Modified midpoint method.
**
** y[]:			dependent variable (vector if we solve more dim. problems)
** dydx[]:		its first derivatives at the starting value x
** nvar:		number of dependent variables y1,...,yn
** xs:			Starting point x
** htot:		total step to be made
** nstep:		number of sub steps
** yout[]:		vector containing the result
** derivs:		function to calculate the right hand side derivatives
*/
void mmid(float y[], float dydx[], int nvar, float xs, float htot, int nstep,
	float yout[], void (*derivs)(float, float[], float[]));

/*
** Polynomial extrapolation routine from the Numerical Recipes.
**
** iest:		number of the call in the sequence of calls
** xest:		input values for x
** yest[]:		input value for y1,..,yn
** yz[]:		extrapolated function values at x=0
** dy[]:		extrapolation errors
** nv:			number of functions
*/
void pzextr(int iest, float xest, float yest[], float yz[], float dy[],
		int nv);

/*
** Bulirsch-Stoer method to integrate ODEs.
**
** y[]:			dependent variable
** dydx[]:		its first derivative at the starting value x
** nv:			number of variables y1,...,yn
** xx:			
** htry:		step size (the algorithm can use a smaller value if needed)
** eps:			required accuracy
** yscal[]:		vector to scale the error
** hdid:		return actual step size
** hnext:		return estimated next step size
** derivs:		function to calculate the right hand side derivatives
*/
void bsstep(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(TILLMATERIAL*, float, float [], float []));

#endif

