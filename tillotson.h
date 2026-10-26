#include <stdio.h>
/*
 * The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_HINCLUDED
#define TILLOTSON_HINCLUDED

/*
 * Version.
 */
#define TILL_VERSION_TEXT    "3.4.4"
#define TILL_VERSION_MAJOR   3
#define TILL_VERSION_MINOR   4
#define TILL_VERSION_PATCH   4

#include <stdlib.h>
#include <assert.h>

/*
 * Error codes.
 */
#define TILL_SUCCESS 0
#define TILL_FAIL   -1

/*
 * Error handling convention (the same as in ANEOSmaterial).
 *
 * Initialization errors (unknown material, invalid table parameters) are fatal. Every evaluation
 * that cannot deliver a value returns exactly TILL_SENTINEL (-1e50) for doubles or TILL_FAIL for
 * status codes, in any unit system, so callers can test for (value < -1e40). Each such failure is
 * counted in material->nLookupFail and reported to stderr only if EOSLIB_VERBOSE is defined.
 */
#define TILL_SENTINEL (-1e50)

#ifdef EOSLIB_VERBOSE
#define TILL_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define TILL_WARN(...) ((void)0)
#endif

#define TILL_FAIL_RETURN(material, ret, ...) do { (material)->nLookupFail++; TILL_WARN(__VA_ARGS__); return (ret); } while (0)

/*
 * Relative tolerance with which the lookup table is treated as a closed interval at its edges,
 * so that a state on an edge survives unit conversions and log/exp round trips.
 */
#ifndef TILL_EDGE_TOL
#define TILL_EDGE_TOL 1e-10
#endif

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
 * Material constants (redefined on 10.10.2017). 0 is reserved for the ideal gas, which is
 * provided by igeos in EOSlib and not by this library.
 */
#define GRANITE     1
#define IRON        2
#define BASALT      3
#define ICE         4
#define WATER       5
#define DUNITE      6
#define FORSTERITE  7

#define TILL_N_MATERIAL_MAX 10

/*
 * Negative pressures are cut to zero in tillPressure, tillPressureSound and the derivatives (see
 * tillPressureCore in tillotson.c). The *NP functions return the raw EOS.
 */
/* The look up table is a 1D array with Lookup(i, j) = Lookup[i*nTableV + j], i in rho, j in s. */
#define TILL_INDEX(m, i, j) (((i)*(m)->nTableV) + (j))

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

/* Define eps so that v <= v_max-eps. */
#define V_EPS 1e-8

/*
 * The isentropes are labelled by v = u(rho0). The table is uniform in the internal variable s,
 *
 *     v(s) = vmin*(exp(s) - 1),   s(v) = log(1 + v/vmin),   s_j = j*dv,
 *
 * which is linear in v for v << vmin and logarithmic above, so that the first cell is at
 * v ~ vmin*dv and the solid/liquid regime (u < us2) is resolved. vmin is set in tillInitLookup to
 * TILL_VMIN_FRAC*vmax. All routines of tillsplint.c take and return s, not v; only
 * tillInitLookup/tillSolveIsentropeLogRho see the physical v.
 */
#ifndef TILL_VMIN_FRAC
#define TILL_VMIN_FRAC 1e-6
#endif

/*
 * One entry of the lookup table: the isentrope with label s_j = j*dv (v = u(rho0) = v(s_j)) at
 * rho_i = rhomin*exp(i*dlogrho); see tillLookupRho, tillLookupS, tillLookupV.
 */
typedef struct till_lookup_entry
{
    double u;
    double u1;      // du/dlogrho = P/rho along the isentrope
    /*
     * Second derivatives of u and u1 with respect to the table variable s (cubic splines along the
     * s axis), set in tillInitSplines.
     */
    double udv2;    // d2u/ds2
    double u1dv2;   // d2u1/ds2
} TILL_LOOKUP_ENTRY;

typedef struct tillMaterial
{
	int iMaterial;	/* What material is it? */
	int nTableRho;	/* Number of entries in the look up table in rho */
	int nTableV;	/* and in s */
	int n;			/* index of rho0 in the rho grid */
	double rhomin;	/* Min value for the lookup table */
	double rhomax;	/* Max value for the lookup table */
	double vmax;	/* Max value for the lookup table */
	double vmin;	/* Scale of the v(s) mapping, see TILL_VMIN_FRAC */
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

	/* the specific heat capacity (assumed constant) */
	double cv;

	/* Grid spacing of the lookup table */
	double dlogrho;
	double dv;		/* step in s (not in v!) */

	/* The lookup table u(rho, s) along the isentropes, NULL until tillInitLookup */
	TILL_LOOKUP_ENTRY *Lookup;

	/* Diagnostics: number of failed evaluations since init (see TILL_FAIL_RETURN). */
	long nLookupFail;
} TILLMATERIAL;

/*
 * tillotson.c
 */
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit);
void tillFinalizeMaterial(TILLMATERIAL *material);

void tilliMatString(TILLMATERIAL *material, char *MatName);
void tillErrorString(int iError, char *ErrorString);
void tillPrintMat(TILLMATERIAL *material, FILE *fp);

/* Pressure, sound speed and derivatives of the cut EOS (negative pressures set to zero). */
double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressure(TILLMATERIAL *material, double rho, double u);
double tillSoundSpeed(TILLMATERIAL *material, double rho, double u);
double tilldPdrho(TILLMATERIAL *material, double rho, double u);
double tilldPdu(TILLMATERIAL *material, double rho, double u);
/* The raw EOS without cutoff (used for inversions). */
double tillPressureSoundNP(TILLMATERIAL *material, double rho, double u, double *pcSound);
double tillPressureNP(TILLMATERIAL *material, double rho, double u);
double tilldudrho(TILLMATERIAL *material, double rho, double u);
double tilldudlogrho(TILLMATERIAL *material, double logrho, double u);

/* Temperature: u(rho, T) = u_cold(rho) + cv*T. These need the lookup table (cold curve). */
double tillTempRhoU(TILLMATERIAL *material, double rho, double u);
double tillURhoTemp(TILLMATERIAL *material, double rho, double T);
/* rho(P, T) by root finding; TILL_SENTINEL if P <= 0 or not bracketed. */
double tillRhoPTemp(TILLMATERIAL *material, double P, double T);

int tillPhase(TILLMATERIAL *material, double rho, double u);

/*
 * Brent root finder for f(material, x, ctx) = 0 on the bracket [a, b] with fa = f(a), fb = f(b)
 * of opposite sign. Converges to |x - root| <= tol_abs + tol_rel*|x|. Returns TILL_SUCCESS and
 * the root in *px, or TILL_FAIL (bracket invalid or max_iter exceeded) with *px untouched.
 */
int tillBrent(double (*f)(TILLMATERIAL *, double, const void *), TILLMATERIAL *material,
              const void *ctx, double a, double b, double fa, double fb,
              double tol_abs, double tol_rel, int max_iter, double *px);

/*
 * tillinitlookup.c
 */
void tillInitLookup(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmax);
TILL_LOOKUP_ENTRY *tillSolveIsentropeLogRho(TILLMATERIAL *material, double v);
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2);
int tillIsInTable(TILLMATERIAL *material,double rho,double u);
/* 1 if (rho, u) is below the cold curve, 0 if not, TILL_FAIL if rho is outside the table. */
int tillIsBelowColdCurve(TILLMATERIAL *material,double rho,double u);

/*
 * tillsplint.c
 */
void tillInitSplines(TILLMATERIAL *material);
void tillInitSplineU(TILLMATERIAL *material);
void tillInitSplineU1(TILLMATERIAL *material);
double tillSplineIntU(TILLMATERIAL *material, double s, int i);
double tillSplineIntU1(TILLMATERIAL *material, double s, int i);

/* u(rho, s) from the table; tillCubicIntD also returns du/ds. rho and s must be inside the table. */
double tillCubicInt(TILLMATERIAL *material, double rho, double s);
double tillCubicIntD(TILLMATERIAL *material, double rho, double s, double *pduds);
/* u(rho, s_j) on grid line j for rho in interval i at fraction e (see tillLocateRho). */
double tillCubicIntGrid(TILLMATERIAL *material, int i, double e, int j);
void tillLocateRho(TILLMATERIAL *material, double rho, int *pi, double *pe);

double tillFindUonIsentrope(TILLMATERIAL *material,double s,double rho);
/* Returns the table variable s of the isentrope through (rho, u), or TILL_SENTINEL. */
double tillFindEntropyCurve(TILLMATERIAL *material,double rho,double u,int iOrder);
/* Returns u2 on the isentrope through (rho1, u1) at rho2, or TILL_SENTINEL. */
double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);
/* u_cold(rho) = u(rho, s=0), or TILL_SENTINEL if rho is outside the table. */
double tillColdULookup(TILLMATERIAL *material,double rho);

int tillLookupIndexLogRho(TILLMATERIAL *material, double logrho);
int tillLookupIndexV(TILLMATERIAL *material, double s);

double tillLookupRho(TILLMATERIAL *material, int i);
double tillLookupLogRho(TILLMATERIAL *material, int i);
double tillLookupV(TILLMATERIAL *material, int j);
double tillLookupS(TILLMATERIAL *material, int j);
double tillVofS(TILLMATERIAL *material, double s);
double tillSofV(TILLMATERIAL *material, double v);

#endif