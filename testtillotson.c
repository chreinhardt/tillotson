/*
 * testtillotson.c: invariant tests for the Tillotson EOS library.
 *
 * No reference values are used: every check is an identity the library must satisfy (consistency
 * between functions, round trips, monotonicity, derivatives against finite differences, the
 * lookup table against a direct integration of the same ODE). Each failed check prints one line;
 * the exit status is the number of failures (0 = pass). Known limitations of the EOS itself
 * (non-monotonic P(rho, T) of some materials) are reported and skipped, not counted.
 *
 *   make testtillotson && ./testtillotson [iMat ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tillotson.h"

static int nFail = 0, nCheck = 0, nSkip = 0;
static const char *MatName = "";

#define CHECK(cond, ...) do { nCheck++; if (!(cond)) { nFail++; \
    fprintf(stderr, "FAIL %s %s:%d: ", MatName, __func__, __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); } } while (0)
#define CHECK_REL(a, b, tol, ...) do { double a_ = (a), b_ = (b); nCheck++; \
    if (!(fabs(a_ - b_) <= (tol)*(fabs(a_) + fabs(b_)) )) { nFail++; \
    fprintf(stderr, "FAIL %s %s:%d: %.10g vs %.10g (rel %.2e > %.1e): ", MatName, __func__, __LINE__, a_, b_, fabs(a_-b_)/(fabs(a_)+fabs(b_)+1e-300), (double)(tol)); \
    fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); } } while (0)
#define SKIP(...) do { nSkip++; fprintf(stderr, "SKIP %s %s: ", MatName, __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); } while (0)

/* Reproducible uniform random numbers in [0, 1) (do not depend on the platform's rand()). */
static unsigned long long rngState = 88172645463325252ULL;
static double urand(void) { rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17; return (rngState >> 11)*(1.0/9007199254740992.0); }
static double lograndom(double a, double b) { return a*pow(b/a, urand()); }

/* The isentrope ODE du/dlogrho = P/rho integrated with RK4 at high resolution, both directions. */
static double refU(TILLMATERIAL *material, double rho1, double u1, double rho2, int nsteps)
{
    double x = log(rho1), h = (log(rho2) - x)/nsteps, u = u1;
    int n;
    for (n = 0; n < nsteps; n++) {
        double k1 = h*tilldudlogrho(material, x, u);
        double k2 = h*tilldudlogrho(material, x + 0.5*h, u + 0.5*k1);
        double k3 = h*tilldudlogrho(material, x + 0.5*h, u + 0.5*k2);
        double k4 = h*tilldudlogrho(material, x + h, u + k3);
        u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        x += h;
    }
    return u;
}

/* ---------------------------------------------------------------------------------------- */

static void testUnits(int iMat)
{
    /* The same physical state in cgs and in code units gives the same physical pressure. */
    const double dKpcUnit = 6.371e8/3.085678e21, dMsolUnit = 5.972e27/1.99e33;
    TILLMATERIAL *cgs = tillInitMaterial(iMat, 0.0, 0.0);
    TILLMATERIAL *code = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    int k;

    CHECK(cgs->dErgPerGmUnit == 1.0 && cgs->dGmPerCcUnit == 1.0, "cgs unit factors");
    CHECK(code->rho0 > 0.0 && code->A > 0.0 && code->us2 > code->us && code->cv > 0.0, "material parameters");
    CHECK_REL(cgs->rho0, code->rho0*code->dGmPerCcUnit, 1e-14, "rho0 conversion");

    for (k = 0; k < 200; k++) {
        double rho_cgs = cgs->rho0*lograndom(0.1, 10.0), u_cgs = lograndom(1e8, 1e13);
        double rho = rho_cgs/code->dGmPerCcUnit, u = u_cgs/code->dErgPerGmUnit;
        double c_cgs, c;
        double P_cgs = tillPressureSound(cgs, rho_cgs, u_cgs, &c_cgs);
        double P = tillPressureSound(code, rho, u, &c);
        CHECK_REL(P_cgs, P*code->dGmPerCcUnit*code->dErgPerGmUnit, 1e-12, "pressure rho= %g u= %g", rho_cgs, u_cgs);
        CHECK_REL(c_cgs, c*sqrt(code->dErgPerGmUnit), 1e-12, "sound speed rho= %g u= %g", rho_cgs, u_cgs);
    }
    tillFinalizeMaterial(cgs);
    tillFinalizeMaterial(code);
}

static void testPressure(TILLMATERIAL *material)
{
    const double cmin = sqrt(material->A/material->rho0);
    const double umax = 50.0*material->us2;
    int a, b;

    /* Reference state. */
    CHECK_REL(tillPressure(material, material->rho0, 0.0) + 1.0, 1.0, 1e-14, "P(rho0, 0) = 0");
    CHECK(tillPhase(material, material->rho0, 0.0) == TILL_PHASE_CONDENSED, "phase at rho0");

    for (a = 0; a < 50; a++) for (b = 0; b < 50; b++) {
        double rho = material->rho0*pow(10.0, -1.5 + 2.5*a/49.0);
        double u = 1e-4*material->us2*pow(umax/(1e-4*material->us2), b/49.0);
        double c, c2, P = tillPressureSound(material, rho, u, &c);
        double Pnp = tillPressureSoundNP(material, rho, u, &c2);
        double h = 1e-5, du = h*fmax(u, 1e-2*material->us2), fdr, fdu, dPdrho = tilldPdrho(material, rho, u), dPdu = tilldPdu(material, rho, u);
        int phase = tillPhase(material, rho, u);

        CHECK(phase >= 1 && phase <= 4, "phase rho= %g u= %g", rho, u);
        CHECK(P == tillPressure(material, rho, u), "tillPressure == tillPressureSound rho= %g u= %g", rho, u);
        CHECK(c == tillSoundSpeed(material, rho, u), "tillSoundSpeed == tillPressureSound rho= %g u= %g", rho, u);
        CHECK(P >= 0.0 && isfinite(P), "P >= 0 rho= %g u= %g P= %g", rho, u, P);
        CHECK(c >= cmin*(1.0 - 1e-14) && isfinite(c), "c >= sqrt(A/rho0) rho= %g u= %g c= %g", rho, u, c);
        CHECK(Pnp <= P*(1.0 + 1e-14) + 1e-300, "P_NP <= P rho= %g u= %g", rho, u);
        if (phase == TILL_PHASE_CONDENSED) CHECK_REL(Pnp, P, 1e-14, "condensed: P_NP == P rho= %g u= %g", rho, u);
        if (P > 0.0 && phase != TILL_PHASE_EXP_COLD) CHECK_REL(Pnp, P, 1e-14, "P > 0: P_NP == P rho= %g u= %g", rho, u);
        if (P == 0.0) CHECK(dPdrho == 0.0 && dPdu == 0.0, "cut: derivatives zero rho= %g u= %g", rho, u);

        /* Derivatives against central differences, away from the region boundaries and the cut. */
        if (fabs(u/material->us - 1.0) < 2.0*h || fabs(u/material->us2 - 1.0) < 2.0*h) continue;
        if (u - du <= 0.0 || fabs((u-du)/material->us - 1.0) < h || fabs((u+du)/material->us - 1.0) < h) continue;
        if (!(P > 0.0 && tillPressure(material, rho*(1-h), u) > 0.0 && tillPressure(material, rho, u - du) > 0.0)) continue;
        fdr = (tillPressure(material, rho*(1+h), u) - tillPressure(material, rho*(1-h), u))/(2*h*rho);
        fdu = (tillPressure(material, rho, u + du) - tillPressure(material, rho, u - du))/(2*du);
        CHECK(fabs(fdr - dPdrho) <= 1e-6*(fabs(P)/rho + material->A/material->rho0), "dP/drho rho= %g u= %g: fd %g analytic %g", rho, u, fdr, dPdrho);
        CHECK(fabs(fdu - dPdu) <= 1e-5*rho, "dP/du rho= %g u= %g: fd %g analytic %g", rho, u, fdu, dPdu);
        {
            double c2fd = fdr + P/(rho*rho)*fdu;
            if (c2fd > cmin*cmin) CHECK_REL(c*c, c2fd, 1e-6, "c^2 identity rho= %g u= %g", rho, u);
        }
    }

    /* Continuity in u across the region boundaries. */
    for (a = 0; a < 20; a++) {
        double rho = material->rho0*lograndom(0.1, 0.999), eps = 1e-9;
        double Pa = tillPressure(material, rho, material->us*(1-eps)), Pb = tillPressure(material, rho, material->us*(1+eps));
        double Pc = tillPressure(material, rho, material->us2*(1-eps)), Pd = tillPressure(material, rho, material->us2*(1+eps));
        CHECK(fabs(Pa - Pb) <= 1e-6*(Pa + Pb) + 1e-10*material->A, "P continuous at us, rho= %g: %g vs %g", rho, Pa, Pb);
        CHECK(fabs(Pc - Pd) <= 1e-6*(Pc + Pd) + 1e-10*material->A, "P continuous at us2, rho= %g: %g vs %g", rho, Pc, Pd);
    }
}

static void testTable(TILLMATERIAL *material)
{
    const int N = material->nTableRho, M = material->nTableV;
    int i, j, k;

    CHECK(material->Lookup != NULL && N > 2 && M > 2, "table allocated");
    CHECK(material->n >= 1 && material->n <= N-2, "rho0 inside the rho grid (n= %d)", material->n);
    CHECK_REL(tillLookupRho(material, material->n), material->rho0, 1e-12, "rho(n) == rho0");
    CHECK_REL(tillLookupRho(material, N-1), material->rhomax, 1e-12, "rho(N-1) == rhomax");
    CHECK_REL(tillLookupRho(material, 0), material->rhomin, 1e-12, "rho(0) == rhomin");
    CHECK_REL(tillLookupV(material, M-1), material->vmax, 1e-10, "v(M-1) == vmax");
    CHECK(tillLookupV(material, 0) == 0.0, "v(0) == 0");

    /* s <-> v mapping. */
    for (k = 0; k < 100; k++) {
        double v = lograndom(1e-9*material->vmax, material->vmax);
        CHECK_REL(tillVofS(material, tillSofV(material, v)), v, 1e-12, "v(s(v))");
    }

    /* Row n is the initial condition, u(rho0, s_j) = v_j; the cold curve is 0 below rho0. */
    for (j = 0; j < M; j++) {
        CHECK_REL(material->Lookup[TILL_INDEX(material, material->n, j)].u + 1.0, tillLookupV(material, j) + 1.0, 1e-12, "u(rho0, s_j) = v_j, j= %d", j);
    }
    for (i = 0; i < material->n; i++) CHECK(material->Lookup[TILL_INDEX(material, i, 0)].u == 0.0, "cold curve zero at i= %d", i);
    for (i = material->n; i < N-1; i++) {
        CHECK(material->Lookup[TILL_INDEX(material, i+1, 0)].u >= material->Lookup[TILL_INDEX(material, i, 0)].u, "cold curve increasing i= %d", i);
    }
    /* u increases with s along every row (isentropes do not cross). */
    for (i = 0; i < N; i += 7) for (j = 0; j < M-1; j++) {
        CHECK(material->Lookup[TILL_INDEX(material, i, j+1)].u > material->Lookup[TILL_INDEX(material, i, j)].u, "u increasing in s, i= %d j= %d", i, j);
    }
    /* u1 = P/rho on the table. */
    for (k = 0; k < 500; k++) {
        i = (int) (urand()*N); j = (int) (urand()*M);
        if (i == N) i--;
        if (j == M) j--;
        if (i == 0) continue;   /* u1(0) is copied from u1(1) by construction */
        {
            double rho = tillLookupRho(material, i), u = material->Lookup[TILL_INDEX(material, i, j)].u;
            /* P/rho is a difference of large terms at high compression (100 rho0): roundoff ~1e-11. */
            CHECK_REL(material->Lookup[TILL_INDEX(material, i, j)].u1 + 1e-300, tillPressure(material, rho, u)/rho + 1e-300, 1e-10, "u1 = P/rho at i= %d j= %d", i, j);
        }
    }

    /*
     * The spline second derivatives solve the natural-spline system of the field they belong to:
     * (y2[j-1] + 4 y2[j] + y2[j+1])/6 = (y[j+1] - 2 y[j] + y[j-1])/h^2, y2[0] = y2[M-1] = 0.
     */
    for (i = 0; i < N; i += 13) {
        const TILL_LOOKUP_ENTRY *row = &material->Lookup[TILL_INDEX(material, i, 0)];
        const double h2 = material->dv*material->dv;
        double scale_u = 0.0, scale_u1 = 0.0;
        CHECK(row[0].udv2 == 0.0 && row[M-1].udv2 == 0.0 && row[0].u1dv2 == 0.0 && row[M-1].u1dv2 == 0.0, "natural b.c. i= %d", i);
        for (j = 1; j < M-1; j++) {
            scale_u = fmax(scale_u, fabs(row[j+1].u - 2.0*row[j].u + row[j-1].u)/h2);
            scale_u1 = fmax(scale_u1, fabs(row[j+1].u1 - 2.0*row[j].u1 + row[j-1].u1)/h2);
        }
        for (j = 1; j < M-1; j++) {
            double lhs = (row[j-1].udv2 + 4.0*row[j].udv2 + row[j+1].udv2)/6.0, rhs = (row[j+1].u - 2.0*row[j].u + row[j-1].u)/h2;
            double lhs1 = (row[j-1].u1dv2 + 4.0*row[j].u1dv2 + row[j+1].u1dv2)/6.0, rhs1 = (row[j+1].u1 - 2.0*row[j].u1 + row[j-1].u1)/h2;
            CHECK(fabs(lhs - rhs) <= 1e-8*scale_u, "spline system of u, i= %d j= %d", i, j);
            CHECK(fabs(lhs1 - rhs1) <= 1e-8*scale_u1, "spline system of u1, i= %d j= %d", i, j);
        }
    }

    /* The interpolant reproduces the grid. */
    for (k = 0; k < 500; k++) {
        i = 1 + (int) (urand()*(N-2)); j = (int) (urand()*(M-1));
        {
            double rho = tillLookupRho(material, i), s = tillLookupS(material, j);
            double u = material->Lookup[TILL_INDEX(material, i, j)].u;
            CHECK_REL(tillCubicInt(material, rho, s) + 1e-300, u + 1e-300, 1e-9, "interpolant on grid i= %d j= %d", i, j);
            CHECK_REL(tillColdULookup(material, rho) + 1.0, material->Lookup[TILL_INDEX(material, i, 0)].u + 1.0, 1e-12, "cold curve on grid i= %d", i);
            CHECK(tillIsInTable(material, rho, u) == TILL_LOOKUP_SUCCESS, "grid point in table i= %d j= %d", i, j);
        }
    }

    /* Edges: closed with the tolerance, open beyond it. */
    {
        double u = tillColdULookup(material, material->rhomax);
        CHECK(tillIsInTable(material, material->rhomin, 0.0) == TILL_LOOKUP_SUCCESS, "rhomin inside (on the cold curve)");
        CHECK(tillIsInTable(material, material->rhomax, u) == TILL_LOOKUP_SUCCESS, "rhomax inside");
        CHECK(tillIsInTable(material, material->rhomin*(1.0 - 3.0*TILL_EDGE_TOL), 1.0) == TILL_LOOKUP_OUTSIDE_RHOMIN, "below rhomin");
        CHECK(tillIsInTable(material, material->rhomax*(1.0 + 3.0*TILL_EDGE_TOL), u) == TILL_LOOKUP_OUTSIDE_RHOMAX, "above rhomax");
        CHECK(tillIsInTable(material, material->rho0, 0.0) == TILL_LOOKUP_SUCCESS, "(rho0, 0) inside");
        CHECK(tillIsInTable(material, material->rho0, material->vmax) == TILL_LOOKUP_SUCCESS, "(rho0, vmax) inside");
        CHECK(tillIsInTable(material, material->rho0, material->vmax*(1.0 + 1e-6)) == TILL_LOOKUP_OUTSIDE_VMAX, "(rho0, >vmax) outside");
        CHECK(tillIsInTable(material, 2.0*material->rho0, 0.5*tillColdULookup(material, 2.0*material->rho0)) == TILL_LOOKUP_OUTSIDE_VMIN, "below cold curve");
        CHECK(tillIsInTable(material, NAN, 1.0) == TILL_LOOKUP_OUTSIDE_RHOMIN, "NaN rho");
        CHECK(tillIsInTable(material, material->rho0, NAN) == TILL_LOOKUP_OUTSIDE_VMAX, "NaN u");
    }
}

static void testIsentropes(TILLMATERIAL *material)
{
    double maxerr = 0.0;
    int k, n = 0;

    for (k = 0; k < 400; k++) {
        /* A state on the isentrope v at rho1, and its lookup. */
        double rho1 = material->rho0*lograndom(0.3, 3.0);
        double v = (k % 2) ? lograndom(1e-6*material->vmax, 0.9*material->vmax) : 0.9*material->vmax*urand();
        double u1 = refU(material, material->rho0, v, rho1, 2000);
        double s, rho2, ref, lk, back;

        if (tillIsInTable(material, rho1, u1) != TILL_LOOKUP_SUCCESS) continue;
        n++;
        /* The isentrope label is recovered. */
        s = tillFindEntropyCurve(material, rho1, u1, k);
        CHECK(!tillFailed(s), "entropy curve found rho= %g u= %g", rho1, u1);
        CHECK_REL(tillVofS(material, s) + 1e-3*material->vmax, v + 1e-3*material->vmax, 1e-5, "isentrope label rho1= %g v= %g", rho1, v);

        /* Lookup against the ODE, compression and expansion. */
        rho2 = (k % 4 < 2) ? 2.0*rho1 : 0.5*rho1;
        ref = refU(material, rho1, u1, rho2, 2000);
        lk = tillLookupU(material, rho1, u1, rho2, k);
        CHECK(!tillFailed(lk), "lookup rho1= %g u1= %g rho2= %g", rho1, u1, rho2);
        CHECK(fabs(lk - ref) <= 1e-3*fabs(ref) + 1e-6*v, "lookup vs ODE rho1= %g v= %g rho2= %g: %g vs %g", rho1, v, rho2, lk, ref);
        if (fabs(lk - ref)/(fabs(ref) + 1e-300) > maxerr) maxerr = fabs(lk - ref)/(fabs(ref) + 1e-300);

        /* Round trip rho1 -> rho2 -> rho1. */
        if (tillIsInTable(material, rho2, lk) == TILL_LOOKUP_SUCCESS) {
            back = tillLookupU(material, rho2, lk, rho1, k);
            CHECK(fabs(back - u1) <= 1e-6*fabs(u1) + 1e-9*v, "round trip rho1= %g v= %g: %g vs %g", rho1, v, back, u1);
        }

        /* Direct integration agrees with the ODE and with the lookup. */
        if (k % 8 == 0) {
            double d = tillCalcU(material, rho1, u1, rho2);
            CHECK(fabs(d - ref) <= 1e-6*fabs(ref) + 1e-9*v, "tillCalcU vs ODE rho1= %g v= %g rho2= %g", rho1, v, rho2);
        }
    }
    CHECK(n > 100, "enough in-table samples (%d)", n);
    fprintf(stdout, "  %s: lookup vs ODE max rel err %.2e (%d samples)\n", MatName, maxerr, n);

    /* u(rho, s) is increasing in s and the lookup is the identity for rho2 = rho1. */
    for (k = 0; k < 100; k++) {
        double rho = material->rho0*lograndom(0.3, 3.0), s1 = urand()*tillLookupS(material, material->nTableV-1), s2 = s1*(1.0 + 1e-3);
        if (s2 > tillLookupS(material, material->nTableV-1)) continue;
        CHECK(tillCubicInt(material, rho, s2) > tillCubicInt(material, rho, s1), "u increasing in s at rho= %g s= %g", rho, s1);
        {
            double u = tillCubicInt(material, rho, s1);
            if (tillIsInTable(material, rho, u) == TILL_LOOKUP_SUCCESS)
                /* The entropy solver converges to TILL_EDGE_TOL*(|u| + 1): relative 1e-8 plus absolute 1e-10. */
                CHECK_REL(tillLookupU(material, rho, u, rho, k) + 1e-2, u + 1e-2, 1e-8, "lookup identity rho= %g", rho);
        }
    }
}

static void testTemperature(TILLMATERIAL *material)
{
    int k, nskip = 0, n = 0;

    for (k = 0; k < 300; k++) {
        double rho = material->rho0*lograndom(0.2, 5.0), T = lograndom(10.0, 1e5);
        double u = tillURhoTemp(material, rho, T), P;
        CHECK_REL(tillTempRhoU(material, rho, u), T, 1e-10, "T(u(T)) rho= %g T= %g", rho, T);
        CHECK(fabs(tillTempRhoU(material, rho, tillColdULookup(material, rho))) <= 1e-12*T, "T on the cold curve is zero, rho= %g", rho);
        CHECK(u >= tillColdULookup(material, rho), "u(T) above the cold curve rho= %g T= %g", rho, T);

        /* rho(P, T) round trip where P(rho, T) is monotonic in rho. */
        P = tillPressure(material, rho, u);
        if (P > 0.0) {
            double r = tillRhoPTemp(material, P, T);
            if (tillFailed(r)) { nskip++; continue; }
            if (fabs(r - rho) > 1e-8*rho) {
                /* Several densities give this P at this T: a property of the EOS, not an error. */
                double ur = tillURhoTemp(material, r, T);
                if (fabs(tillPressureNP(material, r, ur) - P) <= 1e-8*P) { nskip++; continue; }
            }
            CHECK_REL(r, rho, 1e-8, "rho(P(rho,T),T) rho= %g T= %g", rho, T);
            n++;
        }
    }
    if (nskip > 0) SKIP("rho(P, T) not unique or not bracketed for %d of %d states (non-monotonic P(rho, T))", nskip, nskip + n);
}

static void testSentinels(TILLMATERIAL *material)
{
    double ucold = tillColdULookup(material, 2.0*material->rho0);
    long nfail0 = material->nLookupFail;

    CHECK(tillLookupU(material, 2.0*material->rho0, 0.5*ucold, material->rho0, 0) == TILL_FAIL_VALUE, "lookup below cold curve");
    CHECK(tillFindEntropyCurve(material, 2.0*material->rho0, 0.5*ucold, 0) == TILL_FAIL_VALUE, "entropy curve below cold curve");
    CHECK(tillFindEntropyCurve(material, material->rho0, 2.0*material->vmax, 0) == TILL_FAIL_VALUE, "entropy curve above table");
    CHECK(tillColdULookup(material, 2.0*material->rhomax) == TILL_FAIL_VALUE, "cold curve above rhomax");
    CHECK(tillColdULookup(material, 0.5*material->rhomin) == 0.0, "cold curve below rhomin is zero (rho < rho0)");
    CHECK(tillTempRhoU(material, 2.0*material->rhomax, 1.0) == 0.0, "T above rhomax falls back to T = 0 (cold curve)");
    CHECK(tillURhoTemp(material, 2.0*material->rhomax, 1.0) == TILL_FAIL_VALUE, "u(T) above rhomax");
    CHECK(tillRhoPTemp(material, -1.0, 300.0) == TILL_FAIL_VALUE, "rho(P<0, T)");
    CHECK(tillRhoPTemp(material, 0.0, 300.0) == TILL_FAIL_VALUE, "rho(P=0, T)");
    CHECK(tillIsBelowColdCurve(material, 2.0*material->rhomax, 1.0) == TILL_FAIL, "isbelow above rhomax");
    CHECK(tillIsBelowColdCurve(material, 2.0*material->rho0, 0.5*ucold) == 1, "isbelow true");
    CHECK(tillIsBelowColdCurve(material, 2.0*material->rho0, 2.0*ucold) == 0, "isbelow false");
    CHECK(material->nLookupFail == nfail0 + 9, "every failure is counted (%ld)", material->nLookupFail - nfail0);
    /* Outside the table in rho the lookup falls back to direct integration (and is finite, not a failure). */
    CHECK(isfinite(tillLookupU(material, 0.5*material->rhomin, 1.0, material->rhomin, 0)), "lookup below rhomin integrates");
    CHECK(material->nLookupFail == nfail0 + 9, "direct integration is not a failure");
}

static double fcube(TILLMATERIAL *m, double x, const void *ctx) { (void) m; (void) ctx; return x*x*x - 2.0; }

static void testBrent(void)
{
    double x;
    MatName = "brent";
    CHECK(tillBrent(fcube, NULL, NULL, 0.0, 2.0, -2.0, 6.0, 1e-14, 1e-14, 100, &x) == TILL_SUCCESS, "converges");
    CHECK_REL(x, cbrt(2.0), 1e-13, "root of x^3 - 2");
    CHECK(tillBrent(fcube, NULL, NULL, 2.0, 3.0, 6.0, 25.0, 1e-14, 1e-14, 100, &x) == TILL_FAIL, "rejects an unbracketed root");
    CHECK(tillBrent(fcube, NULL, NULL, 0.0, 2.0, -2.0, 6.0, 0.0, 0.0, 2, &x) == TILL_FAIL, "reports max_iter");
}

int main(int argc, char **argv)
{
    const double dKpcUnit = 6.371e8/3.085678e21, dMsolUnit = 5.972e27/1.99e33;
    const int allMat[] = {GRANITE, IRON, BASALT, ICE, WATER, DUNITE, FORSTERITE};
    int nMat = (int) (sizeof(allMat)/sizeof(allMat[0])), m;
    static char name[256];

    testBrent();

    for (m = 0; m < (argc > 1 ? argc-1 : nMat); m++) {
        int iMat = (argc > 1) ? atoi(argv[m+1]) : allMat[m];
        TILLMATERIAL *material = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
        int nFail0 = nFail;

        tilliMatString(material, name);
        MatName = name;
        fprintf(stdout, "Material %d (%s)\n", iMat, name);

        testUnits(iMat);
        testPressure(material);

        /* The standard table, as EOSlib initializes it. */
        tillInitLookupDefault(material);
        CHECK(material->nTableRho == TILL_TABLE_NRHO && material->nTableV == TILL_TABLE_NV, "standard table sizes");
        CHECK(material->vmin == material->cv*TILL_TABLE_TMIN && material->vmax == material->cv*TILL_TABLE_TMAX, "standard table v range");
        CHECK(material->rhomin == TILL_TABLE_RHOMIN_FAC*material->rho0, "standard table rhomin");
        testTable(material);
        testIsentropes(material);
        testTemperature(material);
        testSentinels(material);

        /* Re-initialization with other parameters works (no leak, no stale state). */
        tillInitLookup(material, 300, 200, 1e-3, 50.0, 1e-4, 100.0);
        CHECK(material->nTableRho == 300 && material->nTableV == 200, "re-init sizes");
        CHECK(material->rhomax >= 50.0*exp(-material->dlogrho) && material->rhomax <= 50.0*exp(material->dlogrho), "re-init rhomax= %g adjusted to the grid", material->rhomax);
        testTable(material);

        tillFinalizeMaterial(material);
        fprintf(stdout, "  %s: %d failures\n", name, nFail - nFail0);
    }

    fprintf(stdout, "%d checks, %d failures, %d skipped\n", nCheck, nFail, nSkip);
    return nFail;
}