/*
 * The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_HINCLUDED
#define TILLOTSON_HINCLUDED

#include <stdio.h>   /* FILE in tillPrintMat */

/*
 * Version history
 *   3.4.x  GSL splines and NR Brent; isentrope table uniform in v; TILL_VERBOSE/TILL_OUTPUT_ALL_WARNINGS.
 *   4.0.0  Review release. Fixes: u1 spline was built from u (tillInitSplineU1); the downward
 *          RK4 stages were evaluated at logrho + h/2; the v grid was uniform in v, which put the
 *          whole condensed regime (v < u_s2) in its first cells; the cubic interpolation weight
 *          B used dv instead of j*dv; tilldTdrho/tilldTdu/tilldPdrho_s (wrong, removed); the
 *          derivatives ignored the negative-pressure cut; a spurious -eta term in dP_e/drho and
 *          the missing blend term in dP/du; tillPressureSound returned negative P where
 *          tillPressure returned 0.
 *          Changes: one analytic core tillPressureCore for P, c, dP/drho, dP/du (cut and uncut
 *          variants); the isentrope table is uniform in s = log(1 + v/vmin) (see below) with an
 *          explicit vmin and the standard table defined by TILL_TABLE_* (temperatures at rho0,
 *          densities in units of rho0; tillInitLookupDefault); natural cubic splines along s and
 *          cubic Hermite in logrho, no GSL, no NR code (own Brent); no mallocs in the lookups;
 *          the EOSlib failure contract (TILL_FAIL_VALUE, nLookupFail, TILL_WARN under
 *          EOSLIB_VERBOSE, no exit in any evaluation, tillTempRhoU falls back to T = 0);
 *          status codes TILL_SUCCESS 0 / TILL_FAIL -99 / TILL_LOOKUP_OUTSIDE_* -1..-4;
 *          3-argument TILL_INDEX; ideal gas removed (igeos); tillinitlookup.h, tillsplint.h and
 *          interpol/brent.c removed; test suite testtillotson.c.
 *          tillIsInTable returns TILL_FAIL for a non-finite u (NaN was reported as
 *          OUTSIDE_VMAX, which EOSIsUsable accepts for Tillotson); tillLookupU,
 *          tillFindEntropyCurve and tillIsBelowColdCurve reject a non-finite u as well.
 *          tillLookupU clamps a state below the cold curve to u_cold(rho2) (counted in
 *          nLookupFail) instead of returning TILL_FAIL_VALUE, like ANEOSisentropicU.
 *          No tolerance depends on the unit system: the energy tolerance of tillIsInTable and
 *          tillFindEntropyCurve is TILL_EDGE_TOL*(|u| + vmin) (was |u| + 1 in code units), the
 *          lower bracket of tillRhoPTemp is TILL_RHOPT_RHOMIN_FAC*rho0 (was 1e-10 code units), and
 *          the index of rho0 in the rho grid no longer depends on roundoff in the unit conversion
 *          (it was 799 instead of 800 in some unit systems, which changed the whole grid).
 *          Documented why the cold curve is zero below rho0 in both pressure modes
 *          (TILL_PRESS_MELOSH included), and that the total pressure is cut at zero in both.
 *          tillCalcU rejects a non-positive or non-finite rho1, rho2 and a non-finite u1 with
 *          TILL_FAIL_VALUE (a NaN or infinite rho was undefined behaviour in the step count) and
 *          clamps u >= 0 after every step like the table integration; tillLookupU rejects an
 *          invalid rho2 and counts its direct-integration fallback in nLookupFail.
 *          The table keeps u1 = P/rho at node 0 (rhomin) instead of copying u1 from node 1, so the
 *          first Hermite interval in log(rho) is consistent with the data.
 *          tillPressureCore evaluates w0 = u/(u0 eta^2) + 1 with max(u, 0), so a negative u no
 *          longer divides by zero at u = -u0 eta^2 (the pressure of a state with u < 0, below
 *          the cold curve, uses Gamma = a + b instead of a growing b/w0). dGasConst is k_B/m_H in
 *          cgs as well (was k_B). tillColdULookup fails for rho <= 0 (returned 0, the expanded
 *          cold curve), so u(rho, T) fails and T(rho, u) falls back to 0 there.
 *          Documentation: the unused V_EPS and the MIN/MAX macros are gone from this public
 *          header (MAX was used once, now fmax), <stdio.h> is included inside the include
 *          guard, the spline and index routines say s where they said v, and
 *          tillPrintMat prints nLookupFail.
 *          tillTempRhoU returns T = 0 (counted) below the cold curve and for a non-finite u
 *          instead of the negative (u - u_cold)/cv; within tillUTol of the cold curve T = 0,
 *          not counted. An inconsistent unit pair (exactly one of dKpcUnit, dMsolUnit <= 0, NaN,
 *          inf) is fatal with a message also under -DNDEBUG (it was only an assert).
 *          Hot path: tillInitLookup stores logrhomin, 1/dlogrho and 1/dv, so the index and
 *          spline routines need no log(rhomin) and no division per lookup; the Newton iteration of
 *          tillFindEntropyCurve reuses the location of rho (one log() per call, not one per
 *          iteration); tillPressureCore skips the derivatives when only P is wanted (P is
 *          bit-identical) and evaluates exp(-beta z^2) once.
 *          Isentrope integration: one adaptive RK4 integrator (step doubling, local tolerance
 *          TILL_ODE_RTOL = 1e-11 relative to |u| + vmin) for the table and for tillCalcU, instead
 *          of TILL_ISENTROPE_NSUB = 10 fixed substeps per table interval (table) and fixed steps of
 *          dlogrho/100 (tillCalcU, ~1e5 steps for a long expansion). The table nodes are ~1e4x
 *          more accurate (1e-8 instead of 1e-4 relative) and the table builds in 0.6 s instead
 *          of 2 s per material; TILL_ISENTROPE_NSUB is gone. The remaining lookup error is the
 *          interpolation across the kinks of P/rho (see the standard table below).
 *          tillInitColdCurve(Default) builds only the grid and the two edge isentropes of the
 *          table (cold curve and top isentrope, 77 kB and two integrations instead of 46 MB and
 *          ~0.6 s), which is all T(rho, u), u(rho, T), rho(P, T), the cold curve and tillIsInTable
 *          need; results are bit-identical to the full table. The isentropic lookups need the full
 *          table (tillInitLookup) and are fatal without it (tillRequireLookup). New members
 *          ColdCurve and TopCurve of TILLMATERIAL. EOSlib builds the edge curves in EOSinitMaterial
 *          and the table in EOSinitIsentropicLookup (it built the table in EOSinitMaterial: with one
 *          instance per thread 6 GB and over a minute of startup for 64 threads x 2 materials).
 *          tillURhoTemp returns TILL_FAIL_VALUE (counted) for T < 0 or a non-finite T, like u(rho, T)
 *          of the other libraries below T_min; it returned u_cold + cv T, a state below the cold
 *          curve, which P(rho, T) in EOSlib then evaluated. tillRhoPTemp fails for such a T at once,
 *          counted once (the bracket evaluations would have counted it three times).
 *          tillIsBelowColdCurve fails exactly where tillColdULookup does (rho > rhomax, rho <= 0,
 *          NaN); below rhomin, where the cold curve is zero, it returns u < 0. It failed for every
 *          rho below rhomin, while tillColdULookup and EOSisbelowColdCurve gave a valid answer there.
 */
#define TILL_VERSION_TEXT    "4.0.0"
#define TILL_VERSION_MAJOR   4
#define TILL_VERSION_MINOR   0
#define TILL_VERSION_PATCH   0

#include <stdlib.h>
#include <assert.h>
#include <math.h>

/*
 * Status codes (EOSlib convention: 0 is success, -99 a generic failure; the values -1..-4 are
 * reserved for the OUTSIDE_* codes, see enum till_error_lookup).
 */
#define TILL_SUCCESS 0
#define TILL_FAIL  -99

/*
 * Failure contract (the EOSlib convention, shared by all EOS libraries).
 *
 * A failed evaluation never exits. Every function returning a double returns exactly
 * TILL_FAIL_VALUE (-1e50), in any unit system, when it cannot deliver a value (rho outside the
 * lookup table for the cold curve and u(rho, T); T < 0 or not finite for u(rho, T) and rho(P, T),
 * below the lowest temperature T = 0 of this EOS; (rho, u) below the cold curve or above the table
 * for tillFindEntropyCurve; rho(P, T) not bracketed); status functions return TILL_FAIL. Two
 * functions fall back instead: tillTempRhoU never fails, outside the table in rho, below the cold
 * curve and for a non-finite u it returns T = 0, the cold curve, which is the "lowest tabulated
 * temperature" of this EOS (same contract as ANEOSTofRhoU), and it is never negative; tillLookupU clamps a state below the cold curve to u_cold(rho2), the end point of
 * the coldest isentrope (same convention as ANEOSisentropicU), and integrates the ODE directly
 * (tillCalcU) where the start or end point is outside the table. All of these count the fallback.
 * tillCalcU itself counts only its failures (invalid rho1, rho2 or u1, non-finite result).
 * The pressure, sound speed and their derivatives are analytic and never fail.
 *
 * Every failure is counted in material->nLookupFail and reported on stderr only when the library
 * is compiled with -DEOSLIB_VERBOSE (TILL_WARN). Test a returned value with tillFailed(x) (true
 * for the sentinel and for NaN); EOSlib exposes the same test as EOSFailed(). Initialization
 * errors (unknown material, invalid table parameters) are fatal.
 */
#define TILL_FAIL_VALUE (-1e50)
#define TILL_SENTINEL TILL_FAIL_VALUE   /* name used by 3.x */

static inline int tillFailed(double x) {
    return !(x > -1e40);
}

/*
 * Diagnostic messages of the lookups, only with -DEOSLIB_VERBOSE. The silent variant still
 * type-checks the arguments but compiles to nothing.
 */
#ifdef EOSLIB_VERBOSE
#define TILL_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define TILL_WARN(...) do { if (0) fprintf(stderr, __VA_ARGS__); } while (0)
#endif

#define TILL_FAIL_RETURN(material, ret, ...) do { (material)->nLookupFail++; TILL_WARN(__VA_ARGS__); return (ret); } while (0)

/*
 * Relative tolerance with which the lookup table is treated as a closed interval at its edges,
 * so that a state on an edge survives unit conversions and log/exp round trips.
 */
#ifndef TILL_EDGE_TOL
#define TILL_EDGE_TOL 1e-10
#endif

/*
 * Lower end of the density bracket of tillRhoPTemp in units of rho0 (the raw pressure diverges
 * for rho -> 0 at u = 0). In units of rho0, so the inversion covers the same physical range in
 * every unit system.
 */
#ifndef TILL_RHOPT_RHOMIN_FAC
#define TILL_RHOPT_RHOMIN_FAC 1e-10
#endif

/*
 * Return codes of tillIsInTable (EOSlib maps them one by one). Besides these, tillIsInTable
 * returns TILL_FAIL for a non-finite u (no state at all, mapped to EOS_FAIL by EOSlib).
 */
enum till_error_lookup
{
    TILL_LOOKUP_SUCCESS        =  0,
    TILL_LOOKUP_OUTSIDE_RHOMIN = -1,
    TILL_LOOKUP_OUTSIDE_RHOMAX = -2,
    TILL_LOOKUP_OUTSIDE_VMIN   = -3,
    TILL_LOOKUP_OUTSIDE_VMAX   = -4
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

/*
 * The isentropes are labelled by v = u(rho0). The table is uniform in the internal variable s,
 *
 *     v(s) = vmin*(exp(s) - 1),   s(v) = log(1 + v/vmin),   s_j = j*dv,
 *
 * which is linear in v for v << vmin and logarithmic above (constant relative resolution
 * dv in v/v), so that the solid/liquid regime (u < us2) is resolved as well as the vapour.
 * vmin is the scale of the mapping, not a lower bound: row 0 is the cold curve (v = 0). All
 * routines of tillsplint.c take and return s, not v; only tillInitLookup and
 * tillSolveIsentropeLogRho see the physical v.
 *
 * Standard table (tillInitLookupDefault, used by EOSlib): since u(rho0, T) = cv*T, the v range
 * is given as temperatures at rho0, vmin = cv*TILL_TABLE_TMIN and vmax = cv*TILL_TABLE_TMAX,
 * and the density range in units of rho0. The table is then the same physical table for every
 * material and unit system: 1200/ln(1e7) = 74 rows per decade in T and 1200/ln(1e6) = 87 per
 * decade in rho; 1200 x 1200 entries of 32 bytes = 46 MB per material. Override with -D.
 *
 * Accuracy: the nodes are integrated adaptively (TILL_ODE_RTOL in tillinitlookup.c) and agree with
 * a converged ODE solution to ~1e-8 relative. The lookup error (tillLookupU against the ODE) is set
 * by the interpolation, not by the integration: u(logrho) has a kink where P/rho has one inside an
 * interval, i.e. where the cold term is cut at P = 0 (expanded states in tension, e.g. just below
 * rho0 at low T) and, with TILL_PRESS_MELOSH, where P jumps at 0.8 rho0. There the cubic Hermite
 * interpolation in logrho is only first/second order: the maximum lookup error is ~1e-4 (forsterite,
 * granite 3e-5) in the default mode and ~1e-3 at the Melosh jump; the mean error is ~1e-7. Twice as
 * many rho rows (TILL_TABLE_NRHO 2400, 92 MB) only reduce the maximum by ~1.6x.
 */
#ifndef TILL_TABLE_NRHO
#define TILL_TABLE_NRHO 1200
#endif
#ifndef TILL_TABLE_NV
#define TILL_TABLE_NV 1200
#endif
#ifndef TILL_TABLE_RHOMIN_FAC
#define TILL_TABLE_RHOMIN_FAC 1e-4    /* rho_min = 1e-4 rho0 */
#endif
#ifndef TILL_TABLE_RHOMAX_FAC
#define TILL_TABLE_RHOMAX_FAC 100.0   /* rho_max = 100 rho0 */
#endif
#ifndef TILL_TABLE_TMIN
#define TILL_TABLE_TMIN 1.0           /* K: vmin = cv*T_min */
#endif
#ifndef TILL_TABLE_TMAX
#define TILL_TABLE_TMAX 1e7           /* K: vmax = cv*T_max, the top isentrope */
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
	double vmin;	/* Scale of the v(s) mapping (see above), in code units */
	/* Unit convertion factors */
	double dKpcUnit;
	double dMsolUnit;
	
	double dGasConst;      /* k_B/m_H in code units (specific energy per K), as in igeos */
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
	/* Derived from the grid by tillInitLookup, so that a lookup needs no log() or division for it */
	double logrhomin;	/* log(rhomin) */
	double invdlogrho;	/* 1/dlogrho */
	double invdv;		/* 1/dv */

	/* The lookup table u(rho, s) along the isentropes, NULL until tillInitLookup */
	TILL_LOOKUP_ENTRY *Lookup;
	/*
	 * The two edge isentropes of the table, nTableRho entries each: row j = 0 (the cold curve,
	 * v = 0) and row j = nTableV-1 (the top isentrope, v = vmax). They are all that T(rho, u),
	 * u(rho, T), rho(P, T), the cold curve and tillIsInTable need, so tillInitColdCurve builds only
	 * these (77 kB, two isentrope integrations, instead of the 46 MB table) and the full table is
	 * built on demand by tillInitLookup. Integrated by the same tillSolveIsentropeLogRho call as the
	 * table rows, so they are bit-identical to rows 0 and nTableV-1 of the table. NULL before.
	 */
	TILL_LOOKUP_ENTRY *ColdCurve;
	TILL_LOOKUP_ENTRY *TopCurve;

	/* Diagnostics: number of failed evaluations since init (see TILL_FAIL_RETURN). */
	long nLookupFail;
} TILLMATERIAL;

/*
 * Roundoff tolerance on an internal energy u (code units) at the table edges and in the isentrope
 * root finder: relative to |u|, with the floor vmin = cv*TILL_TABLE_TMIN (the energy scale of the
 * table, 1 K at rho0 for the standard table) instead of an absolute constant, so that the
 * tolerance is the same physical energy in every unit system. Needs the table (vmin > 0).
 */
static inline double tillUTol(const TILLMATERIAL *material, double u)
{
    return TILL_EDGE_TOL*(fabs(u) + material->vmin);
}

/*
 * tillotson.c
 */
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit);
void tillFinalizeMaterial(TILLMATERIAL *material);

void tilliMatString(TILLMATERIAL *material, char *MatName);
const char *tillLookupStatusString(int iError);
void tillErrorString(int iError, char *ErrorString);   /* copies tillLookupStatusString into the buffer */
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

/* Temperature: u(rho, T) = u_cold(rho) + cv*T. These need the cold curve (tillInitColdCurve or tillInitLookup). tillTempRhoU never fails and is >= 0; tillURhoTemp fails for T < 0 or a non-finite T. */
double tillTempRhoU(TILLMATERIAL *material, double rho, double u);
double tillURhoTemp(TILLMATERIAL *material, double rho, double T);
/* rho(P, T) by root finding; TILL_FAIL_VALUE if P <= 0, T < 0 or not finite, or not bracketed. */
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
/*
 * Build the isentrope table: nTableRho log-spaced densities in [rhomin, rhomax] (rhomax is moved
 * onto the grid so that rho0 is a grid point), nTableV isentropes with v = u(rho0) from 0 to
 * vmax, spaced in s = log(1 + v/vmin). All in code units. Replaces an existing table and builds
 * the edge curves (see tillInitColdCurve) as well. 46 MB and ~0.6 s for the standard table.
 */
void tillInitLookup(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmin, double vmax);
/* The standard table TILL_TABLE_* (see above) for this material. */
void tillInitLookupDefault(TILLMATERIAL *material);
/*
 * The grid of the table and its two edge isentropes only, the cold curve (row 0) and the top
 * isentrope (row nTableV-1), with the same arguments as tillInitLookup and the same grid: enough
 * for T(rho, u), u(rho, T), rho(P, T), tillColdULookup, tillIsInTable and tillIsBelowColdCurve,
 * bit-identical to the full table. The isentropic lookups (tillLookupU, tillFindEntropyCurve,
 * tillFindUonIsentrope, tillCubicInt*) need the full table and stop the program (fatal, with a
 * message) until tillInitLookup is called. Frees an existing full table.
 */
void tillInitColdCurve(TILLMATERIAL *material, int nTableRho, int nTableV, double rhomin, double rhomax, double vmin, double vmax);
/* The edge curves of the standard table TILL_TABLE_*. */
void tillInitColdCurveDefault(TILLMATERIAL *material);
/* Fatal (message, assert(0), exit(1)) if the full isentrope table of tillInitLookup is not built; chFunc names the caller. */
void tillRequireLookup(const TILLMATERIAL *material, const char *chFunc);
TILL_LOOKUP_ENTRY *tillSolveIsentropeLogRho(TILLMATERIAL *material, double v);
/* RK4 integration of du/dlogrho = P/rho from (rho1, u1) to rho2, u >= 0; TILL_FAIL_VALUE for invalid input. */
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2);
/* TILL_LOOKUP_SUCCESS, a TILL_LOOKUP_OUTSIDE_* code, or TILL_FAIL if u is not finite. */
int tillIsInTable(TILLMATERIAL *material,double rho,double u);
/* 1 if (rho, u) is below the cold curve, 0 if not, TILL_FAIL where tillColdULookup fails (rho > rhomax, rho <= 0, NaN) or u is not finite; below rhomin the cold curve is 0. */
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
/* Returns the table variable s of the isentrope through (rho, u), or TILL_FAIL_VALUE (also for a non-finite u). */
double tillFindEntropyCurve(TILLMATERIAL *material,double rho,double u,int iOrder);
/* Returns u2 on the isentrope through (rho1, u1) at rho2; u_cold(rho2) if (rho1, u1) is below the cold curve, TILL_FAIL_VALUE for a non-finite u1. */
double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);
/* u_cold(rho) = u(rho, s=0); 0 for 0 < rho < rho0 (also below rhomin), TILL_FAIL_VALUE for rho > rhomax, rho <= 0 or NaN. */
double tillColdULookup(TILLMATERIAL *material,double rho);

int tillLookupIndexLogRho(TILLMATERIAL *material, double logrho);
/* Index j with s_j <= s < s_j+1 (clamped to [0, nTableV-2]); the argument is s, not v (historical name). */
int tillLookupIndexV(TILLMATERIAL *material, double s);

double tillLookupRho(TILLMATERIAL *material, int i);
double tillLookupLogRho(TILLMATERIAL *material, int i);
double tillLookupV(TILLMATERIAL *material, int j);
double tillLookupS(TILLMATERIAL *material, int j);
double tillVofS(TILLMATERIAL *material, double s);
double tillSofV(TILLMATERIAL *material, double v);

#endif