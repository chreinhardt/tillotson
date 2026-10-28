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
 *
 * Build it with -DTILL_PRESS_MELOSH as well to test that pressure mode (same invariants; the
 * isentrope accuracy is lower there, see TILL_TEST_LABEL_TOL).
 */
/* fork() and waitpid() for the fatal paths (POSIX). */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
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
    /* dGasConst is k_B/m_H in every unit system (it was k_B in cgs). */
    CHECK_REL(cgs->dGasConst, 1.38e-16/1.67e-24, 1e-14, "dGasConst = k_B/m_H in cgs");
    CHECK_REL(code->dGasConst*code->dErgPerGmUnit, cgs->dGasConst, 1e-12, "dGasConst conversion");

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

/*
 * The same physical states in three unit systems must give the same physical answers: cgs, Earth
 * units (u ~ 1e-2..1e2 in code units) and a "heavy" system with GM/L ~ 1e18 erg/g, where code
 * energies are ~1e-8 and code densities ~1e9 g/cc. No tolerance may be an absolute constant in
 * code units. A reduced table with the physical parameters of the standard table keeps it fast.
 */
static void testUnitInvariance(int iMat)
{
    const double units[3][2] = {
        {0.0, 0.0},                                         /* cgs */
        {6.371e8/3.085678e21, 5.972e27/1.99e33},            /* Earth radius and mass */
        {1.3e8/3.085678e21, 1.0},                           /* one solar mass in 1300 km: GM/L ~ 1e18 erg/g */
    };
    TILLMATERIAL *m[3];
    int k, iu, nRhoPT = 0;
    /* Own random stream, so that the samples of the other tests do not depend on this one. */
    const unsigned long long rngSaved = rngState;
    rngState = 0x9E3779B97F4A7C15ULL;

    for (iu = 0; iu < 3; iu++) {
        m[iu] = tillInitMaterial(iMat, units[iu][0], units[iu][1]);
        tillInitLookup(m[iu], 300, 300, TILL_TABLE_RHOMIN_FAC*m[iu]->rho0, TILL_TABLE_RHOMAX_FAC*m[iu]->rho0,
                       m[iu]->cv*TILL_TABLE_TMIN, m[iu]->cv*TILL_TABLE_TMAX);
    }
    CHECK(m[2]->dErgPerGmUnit > 1e17, "heavy unit system: dErgPerGmUnit= %g", m[2]->dErgPerGmUnit);
    /* The same grid in every unit system: rho0 at the same index, the same spacing in log(rho) and s. */
    for (iu = 1; iu < 3; iu++) {
        CHECK(m[iu]->n == m[0]->n, "index of rho0 in unit system %d: %d vs %d", iu, m[iu]->n, m[0]->n);
        CHECK_REL(m[iu]->dlogrho, m[0]->dlogrho, 1e-12, "dlogrho in unit system %d", iu);
        CHECK_REL(m[iu]->rhomax*m[iu]->dGmPerCcUnit, m[0]->rhomax, 1e-12, "rhomax in unit system %d", iu);
        CHECK_REL(m[iu]->dv, m[0]->dv, 1e-12, "dv in unit system %d", iu);
    }

    for (k = 0; k < 300; k++) {
        /* A physical state (cgs) on the condensed or expanded side, and a target density. */
        double rho_cgs = m[0]->rho0*lograndom(0.3, 3.0), T = lograndom(10.0, 3e4);
        double u_cgs = tillURhoTemp(m[0], rho_cgs, T), rho2_cgs = rho_cgs*lograndom(0.5, 2.0);
        double s[3], u2[3];
        int st[3];

        for (iu = 0; iu < 3; iu++) {
            double rho = rho_cgs/m[iu]->dGmPerCcUnit, u = u_cgs/m[iu]->dErgPerGmUnit, rho2 = rho2_cgs/m[iu]->dGmPerCcUnit;
            st[iu] = tillIsInTable(m[iu], rho, u);
            s[iu] = (st[iu] == TILL_LOOKUP_SUCCESS) ? tillFindEntropyCurve(m[iu], rho, u, k) : 0.0;
            u2[iu] = tillLookupU(m[iu], rho, u, rho2, k)*m[iu]->dErgPerGmUnit;
        }
        for (iu = 1; iu < 3; iu++) {
            CHECK(st[iu] == st[0], "table status in unit system %d: %d vs %d (rho= %g u= %g)", iu, st[iu], st[0], rho_cgs, u_cgs);
            /* s is dimensionless; the root finder converges to TILL_EDGE_TOL in u. */
            CHECK(fabs(s[iu] - s[0]) <= 1e-8*(1.0 + s[0]), "isentrope label in unit system %d: %.12g vs %.12g (rho= %g T= %g)", iu, s[iu], s[0], rho_cgs, T);
            CHECK_REL(u2[iu] + m[0]->vmin, u2[0] + m[0]->vmin, 1e-8, "isentrope lookup in unit system %d (rho= %g T= %g)", iu, rho_cgs, T);
        }

        /* rho(P, T) on the same physical bracket, also for strongly expanded states (down to 1e-3 rho0). */
        {
            double r_cgs = m[0]->rho0*lograndom(1e-3, 3.0);
            double P_cgs = tillPressure(m[0], r_cgs, tillURhoTemp(m[0], r_cgs, T));
            double r[3];
            if (!(P_cgs > 0.0)) continue;
            for (iu = 0; iu < 3; iu++) {
                double rr = tillRhoPTemp(m[iu], P_cgs/(m[iu]->dGmPerCcUnit*m[iu]->dErgPerGmUnit), T);
                r[iu] = tillFailed(rr) ? rr : rr*m[iu]->dGmPerCcUnit;
            }
            for (iu = 1; iu < 3; iu++) {
                CHECK(tillFailed(r[iu]) == tillFailed(r[0]), "rho(P, T) fails in one unit system only (%d): rho= %g T= %g", iu, r_cgs, T);
                if (!tillFailed(r[iu]) && !tillFailed(r[0])) CHECK_REL(r[iu], r[0], 1e-8, "rho(P, T) in unit system %d (rho= %g T= %g)", iu, r_cgs, T);
            }
            if (!tillFailed(r[0])) nRhoPT++;
        }
    }
    CHECK(nRhoPT > 50, "enough rho(P, T) samples (%d)", nRhoPT);
    for (iu = 0; iu < 3; iu++) tillFinalizeMaterial(m[iu]);
    rngState = rngSaved;
}

/*
 * tillColdULookup returns 0 below rho0 without reading the table. That holds for any material in
 * both pressure modes (P >= 0 can only lower u when expanding from u(rho0) = 0, and u is clamped
 * at 0), also where the cold term A mu + B mu^2 is positive: test it with B = 3A.
 */
static void testColdCurveExpanded(int iMat)
{
    const double dKpcUnit = 6.371e8/3.085678e21, dMsolUnit = 5.972e27/1.99e33;
    TILLMATERIAL *m = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    int i, k;

    m->B = 3.0*m->A;   /* cold term positive at u = 0 for rho < 2/3 rho0 */
#ifndef TILL_PRESS_MELOSH
    CHECK(tillPressure(m, 0.5*m->rho0, 0.0) > 0.0, "B = 3A: positive cold pressure at 0.5 rho0");   /* cut below 0.8 rho0 with the flag */
#endif
    tillInitLookup(m, 300, 300, TILL_TABLE_RHOMIN_FAC*m->rho0, TILL_TABLE_RHOMAX_FAC*m->rho0,
                   m->cv*TILL_TABLE_TMIN, m->cv*TILL_TABLE_TMAX);
    for (i = 0; i <= m->n; i++) CHECK(m->Lookup[TILL_INDEX(m, i, 0)].u == 0.0, "B = 3A: cold curve zero at i= %d", i);
    for (k = 0; k < 100; k++) {
        double rho = m->rho0*lograndom(1e-4, 0.999);
        CHECK(tillColdULookup(m, rho) == 0.0, "B = 3A: u_cold(rho= %g) = 0", rho);
        CHECK(fabs(tillCubicInt(m, rho, 0.0)) <= 1e-14*m->vmin, "B = 3A: table row 0 at rho= %g is %g", rho, tillCubicInt(m, rho, 0.0));
    }
    /*
     * Direct integration clamps u >= 0 like the table: expanding from (rho0, 0) or from just above
     * the cold curve ends on the cold curve, u = 0. Without the clamp the positive cold pressure
     * drives u below zero (the unclamped reference shows it). With TILL_PRESS_MELOSH the pressure of
     * expanded cold states is zero below 0.8 rho0, so a state above the cold curve stays above it.
     */
#ifndef TILL_PRESS_MELOSH
    CHECK(refU(m, m->rho0, 0.0, 0.5*m->rho0, 2000) < 0.0, "B = 3A: the unclamped ODE goes below u = 0");
#endif
    for (k = 0; k < 20; k++) {
        double rho2 = 0.6*m->rho0*pow(1e-3/0.6, k/19.0);   /* below 2/3 rho0; deterministic, keeps the RNG stream of the other tests */
        CHECK(tillCalcU(m, m->rho0, 0.0, rho2) == 0.0, "B = 3A: tillCalcU from (rho0, 0) to rho2= %g is %g", rho2, tillCalcU(m, m->rho0, 0.0, rho2));
#ifndef TILL_PRESS_MELOSH
        CHECK(tillCalcU(m, m->rho0, 1e-3*m->vmin, rho2) == 0.0, "B = 3A: tillCalcU from (rho0, 1e-3 vmin) to rho2= %g", rho2);
#else
        /* The cold-region pressure is zero below 0.8 rho0 (Melosh cut), so u only decreases between rho0 and 0.8 rho0. */
        {
            double u2 = tillCalcU(m, m->rho0, 1e-3*m->vmin, rho2);
            CHECK(u2 >= 0.0 && u2 < 1e-3*m->vmin, "B = 3A: tillCalcU from (rho0, 1e-3 vmin) to rho2= %g is %g", rho2, u2);
        }
#endif
        CHECK(tillCalcU(m, m->rho0, 0.0, rho2) == tillLookupU(m, m->rho0, 0.0, rho2, k), "B = 3A: direct integration == lookup on the cold curve");
    }
    tillFinalizeMaterial(m);
}

/*
 * u < 0 is below the cold curve and not a Tillotson state, but a hydro code can produce it. w0 =
 * u/(u0 eta^2) + 1 is evaluated with max(u, 0), so Gamma = a + b for u <= 0 and the raw pressure is
 * P_NP = (a + b) u rho + A mu + B mu^2 for every rho (u < 0 < us is always the cold region). Before,
 * w0 reached 0 at u = -u0 eta^2 (division by zero) and became negative below (Gamma = a + b/w0 < a).
 */
static void testNegativeU(TILLMATERIAL *material)
{
    const double rho0 = material->rho0, A = material->A, B = material->B, ab = material->a + material->b;
    int a, b;

    for (a = 0; a < 40; a++) for (b = 0; b < 8; b++) {
        static const double f[] = {1e-12, 1e-3, 0.5, 1.0 - 1e-12, 1.0, 1.0 + 1e-12, 2.0, 10.0};
        double rho = rho0*pow(10.0, -2.0 + 2.5*a/39.0), eta = rho/rho0, mu = eta - 1.0;
        double u = -f[b]*material->u0*eta*eta;   /* f = 1: the former singularity w0 = 0 */
        double Gu = ab*u*rho, Pnp = tillPressureNP(material, rho, u), c, P = tillPressureSound(material, rho, u, &c);
        double scale = fabs(Gu) + A*fabs(mu) + B*mu*mu;
        CHECK(isfinite(Pnp) && fabs(Pnp - (Gu + A*mu + B*mu*mu)) <= 1e-12*scale, "u < 0: P_NP = (a+b) u rho + A mu + B mu^2 at rho= %g u= %g: %g", rho, u, Pnp);
        CHECK(isfinite(P) && P >= 0.0 && isfinite(c) && c > 0.0, "u < 0: P= %g c= %g finite at rho= %g u= %g", P, c, rho, u);
        CHECK(isfinite(tilldPdrho(material, rho, u)) && isfinite(tilldPdu(material, rho, u)), "u < 0: derivatives finite at rho= %g u= %g", rho, u);
#ifdef TILL_PRESS_MELOSH
        if (rho < 0.8*rho0) CHECK(P == 0.0, "u < 0, Melosh cut: P = 0 at rho= %g", rho); else
#endif
        CHECK(P == fmax(Pnp, 0.0), "u < 0: P == max(P_NP, 0) at rho= %g u= %g", rho, u);
    }

    /* Continuous at u = 0, and the derivatives are exact for u < 0 (compressed, P > 0, no cut). */
    for (a = 0; a < 20; a++) {
        double rho = rho0*(1.2 + 0.1*a), mu = rho/rho0 - 1.0;
        double u = -0.01*(A*mu + B*mu*mu)/(ab*rho), du = 1e-6*fabs(u), dr = 1e-6*rho;
        double P = tillPressure(material, rho, u);
        CHECK(P > 0.0, "compressed state with u < 0 has P > 0 (rho= %g)", rho);
        {
            /* Continuity: the change over -eps is first order, bounded by eps*dP/du(u = 0) (+ roundoff). */
            double eps = 1e-8*material->u0, P0 = tillPressure(material, rho, 0.0);
            CHECK(fabs(tillPressure(material, rho, -eps) - P0) <= 1.01*eps*tilldPdu(material, rho, 0.0) + 1e-14*P0, "P continuous at u = 0, rho= %g", rho);
        }
        CHECK_REL((tillPressure(material, rho, u + du) - tillPressure(material, rho, u - du))/(2.0*du), tilldPdu(material, rho, u), 1e-6, "dP/du at u < 0, rho= %g", rho);
        CHECK_REL((tillPressure(material, rho + dr, u) - tillPressure(material, rho - dr, u))/(2.0*dr), tilldPdrho(material, rho, u), 1e-6, "dP/drho at u < 0, rho= %g", rho);
    }
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
        /* The cut pressure is the raw one cut at zero; with TILL_PRESS_MELOSH expanded cold states below 0.8 rho0 have P = 0. */
#ifdef TILL_PRESS_MELOSH
        if (phase == TILL_PHASE_EXP_COLD && rho < 0.8*material->rho0) CHECK(P == 0.0, "Melosh cut: P = 0 rho= %g u= %g", rho, u);
        else
#endif
        CHECK(P == fmax(Pnp, 0.0), "P == max(P_NP, 0) rho= %g u= %g: %g vs %g", rho, u, P, Pnp);
        if (phase == TILL_PHASE_CONDENSED) CHECK_REL(Pnp, P, 1e-14, "condensed: P_NP == P rho= %g u= %g", rho, u);
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
        /* The blend weight changes by ~eps*us/(us2-us) across us, so the jump is of that order times P_e (<~ A). */
        CHECK(fabs(Pa - Pb) <= 1e-6*(Pa + Pb) + 1e-8*material->A, "P continuous at us, rho= %g: %g vs %g", rho, Pa, Pb);
        CHECK(fabs(Pc - Pd) <= 1e-6*(Pc + Pd) + 1e-8*material->A, "P continuous at us2, rho= %g: %g vs %g", rho, Pc, Pd);
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
        {
            double rho = tillLookupRho(material, i), u = material->Lookup[TILL_INDEX(material, i, j)].u;
            /* Where the downward integration clamps u at 0 the stored slope is 0 (flat curve). */
            double u1ref = (i < material->n && u == 0.0) ? 0.0 : tillPressure(material, rho, u)/rho;
            /* P/rho is a difference of large terms at high compression (100 rho0): roundoff ~1e-11. */
            CHECK_REL(material->Lookup[TILL_INDEX(material, i, j)].u1 + 1e-300, u1ref + 1e-300, 1e-10, "u1 = P/rho at i= %d j= %d", i, j);
        }
    }
    /* Node 0 (rhomin) has its own slope on every isentrope, not a copy of node 1. */
    for (j = 0; j < M; j++) {
        double rho = tillLookupRho(material, 0), u = material->Lookup[TILL_INDEX(material, 0, j)].u;
        double u1ref = (u == 0.0) ? 0.0 : tillPressure(material, rho, u)/rho;
        CHECK(isfinite(material->Lookup[TILL_INDEX(material, 0, j)].u1), "u1 finite at rhomin, j= %d", j);
        CHECK_REL(material->Lookup[TILL_INDEX(material, 0, j)].u1 + 1e-300, u1ref + 1e-300, 1e-10, "u1 = P/rho at rhomin, j= %d", j);
    }
    /*
     * The first interval [rhomin, rhomin e^dlogrho] interpolates the isentrope as well as any other:
     * against the ODE from node 1 to within 3e-6 (|u| + vmin). With u1(0) copied from node 1 the
     * error was up to ~1e-5 for the built-in materials (2e-7 to 7e-7 with the own slope).
     */
    if (N == TILL_TABLE_NRHO && M == TILL_TABLE_NV) {
        const double lr0 = log(tillLookupRho(material, 0));
        double maxerr = 0.0;
        for (j = 0; j < M; j += 3) for (k = 1; k < 10; k++) {
            double s = tillLookupS(material, j), u1n = material->Lookup[TILL_INDEX(material, 1, j)].u;
            double rho = exp(lr0 + 0.1*k*material->dlogrho);
            double ref = refU(material, tillLookupRho(material, 1), u1n, rho, 1000);
            double err = fabs(tillCubicInt(material, rho, s) - ref)/(fabs(ref) + material->vmin);
            if (err > maxerr) maxerr = err;
            CHECK(err <= 3e-6, "first interval rho= %g j= %d: %g vs %g", rho, j, tillCubicInt(material, rho, s), ref);
        }
        fprintf(stdout, "  %s: first rho interval max rel err %.2e\n", MatName, maxerr);
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
        /* A non-finite u is no state at all (TILL_FAIL), not "above the table": EOSIsUsable tolerates OUTSIDE_VMAX. */
        CHECK(tillIsInTable(material, material->rho0, NAN) == TILL_FAIL, "NaN u");
        CHECK(tillIsInTable(material, material->rho0, INFINITY) == TILL_FAIL, "+inf u");
        CHECK(tillIsInTable(material, material->rho0, -INFINITY) == TILL_FAIL, "-inf u");
    }
}

/*
 * With TILL_PRESS_MELOSH the pressure jumps at 0.8 rho0 for expanded cold and intermediate states
 * (the positive cold term is cut there), so the RK4 integration of the table and of the reference
 * is only first order across the jump: the isentrope label and the lookup are then accurate to
 * ~1e-4 and ~1e-3 instead of ~1e-5 and ~1e-4.
 */
#ifdef TILL_PRESS_MELOSH
#define TILL_TEST_LABEL_TOL  3e-4
#define TILL_TEST_LOOKUP_TOL 3e-3
#else
#define TILL_TEST_LABEL_TOL  1e-5
#define TILL_TEST_LOOKUP_TOL 1e-3
#endif

/*
 * The table nodes are integrated adaptively (TILL_ODE_RTOL). Against a fixed-step RK4 reference
 * from rho0 with 2e5 steps and the same clamp u >= 0 (the reference converges to ~1e-7 relative
 * across the kinks of P/rho; with 10 fixed substeps per interval the nodes were off by ~1e-4).
 * With TILL_PRESS_MELOSH P/rho jumps at 0.8 rho0, where a fixed-step reference is only first
 * order (it converges like 1/nsteps to the table value, 1e-6 at 2e5 steps for basalt), hence the
 * wider tolerance there.
 */
#ifdef TILL_PRESS_MELOSH
#define TILL_TEST_NODE_TOL 1e-5
#else
#define TILL_TEST_NODE_TOL 1e-6
#endif
static void testTableNodes(TILLMATERIAL *material)
{
    double maxerr = 0.0;
    int k;

    for (k = 0; k < 8; k++) {
        int i = (int) (urand()*material->nTableRho), j = (int) (urand()*material->nTableV), n;
        double v, rho, x, h, u, err;
        if (i >= material->nTableRho) i = material->nTableRho - 1;
        if (j >= material->nTableV) j = material->nTableV - 1;
        v = tillLookupV(material, j);
        rho = tillLookupRho(material, i);
        x = log(material->rho0);
        h = (log(rho) - x)/200000;
        u = v;
        for (n = 0; n < 200000; n++) {
            double k1 = h*tilldudlogrho(material, x, u);
            double k2 = h*tilldudlogrho(material, x + 0.5*h, u + 0.5*k1);
            double k3 = h*tilldudlogrho(material, x + 0.5*h, u + 0.5*k2);
            double k4 = h*tilldudlogrho(material, x + h, u + k3);
            u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
            if (u < 0.0) u = 0.0;
            x += h;
        }
        err = fabs(material->Lookup[TILL_INDEX(material, i, j)].u - u)/(fabs(u) + material->vmin);
        if (err > maxerr) maxerr = err;
        CHECK(err <= TILL_TEST_NODE_TOL, "table node i= %d (rho/rho0= %g) j= %d: %.12e vs %.12e", i, rho/material->rho0, j, material->Lookup[TILL_INDEX(material, i, j)].u, u);
    }
    fprintf(stdout, "  %s: table nodes vs ODE max rel err %.2e\n", MatName, maxerr);
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
        CHECK_REL(tillVofS(material, s) + 1e-3*material->vmax, v + 1e-3*material->vmax, TILL_TEST_LABEL_TOL, "isentrope label rho1= %g v= %g", rho1, v);

        /* Lookup against the ODE, compression and expansion. */
        rho2 = (k % 4 < 2) ? 2.0*rho1 : 0.5*rho1;
        ref = refU(material, rho1, u1, rho2, 2000);
        lk = tillLookupU(material, rho1, u1, rho2, k);
        CHECK(!tillFailed(lk), "lookup rho1= %g u1= %g rho2= %g", rho1, u1, rho2);
        CHECK(fabs(lk - ref) <= TILL_TEST_LOOKUP_TOL*fabs(ref) + 1e-6*v, "lookup vs ODE rho1= %g v= %g rho2= %g: %g vs %g", rho1, v, rho2, lk, ref);
        if (fabs(lk - ref)/(fabs(ref) + 1e-300) > maxerr) maxerr = fabs(lk - ref)/(fabs(ref) + 1e-300);

        /* Round trip rho1 -> rho2 -> rho1. */
        if (tillIsInTable(material, rho2, lk) == TILL_LOOKUP_SUCCESS) {
            back = tillLookupU(material, rho2, lk, rho1, k);
            CHECK(fabs(back - u1) <= 1e-6*fabs(u1) + 1e-9*v, "round trip rho1= %g v= %g: %g vs %g", rho1, v, back, u1);
        }

        /* Direct integration agrees with the ODE and with the lookup. */
        if (k % 8 == 0) {
            double d = tillCalcU(material, rho1, u1, rho2);
            CHECK(fabs(d - ref) <= 1e-6*fabs(ref) + 1e-9*v, "tillCalcU vs ODE rho1= %g v= %g rho2= %g: %.10e vs %.10e (u1= %g, us= %g)", rho1, v, rho2, d, ref, u1, material->us);
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
                /* The entropy solver converges to tillUTol = TILL_EDGE_TOL*(|u| + vmin): relative, with the floor vmin. */
                CHECK_REL(tillLookupU(material, rho, u, rho, k) + material->vmin, u + material->vmin, 1e-8, "lookup identity rho= %g", rho);
        }
    }
}

/*
 * tillCalcU is adaptive (step doubling, TILL_ODE_RTOL): long paths outside the table, where it is
 * the fallback of tillLookupU, against a fixed-step RK4 reference with 10^6 steps. The paths cross
 * the kinks of P/rho (pressure cut, us, us2) where the fixed-step reference itself converges only
 * like h^2, so the tolerance is 1e-6 (the reference is good to ~1e-8). Round trips there and back.
 */
static void testCalcUOutside(TILLMATERIAL *material)
{
    const double v[] = {0.0, 1e2*material->vmin, 1e-2*material->vmax, 0.5*material->vmax};
    int k;

    for (k = 0; k < 4; k++) {
        /* Expansion from rho0 far below the table, and compression from rhomax far above it. */
        double rho2 = 1e-3*material->rhomin;
        double u2 = tillCalcU(material, material->rho0, v[k], rho2);
        double ref = refU(material, material->rho0, v[k], rho2, 1000000);
        CHECK(!tillFailed(u2) && fabs(u2 - ref) <= 1e-6*fabs(ref) + 1e-6*material->vmin,
              "tillCalcU expansion to %g rho0, v= %g: %.10e vs %.10e", rho2/material->rho0, v[k], u2, ref);

        {
            double u1 = tillColdULookup(material, material->rhomax) + v[k];
            double rho3 = 10.0*material->rhomax;
            double u3 = tillCalcU(material, material->rhomax, u1, rho3);
            double ref3 = refU(material, material->rhomax, u1, rho3, 1000000);
            double back = tillCalcU(material, rho3, u3, material->rhomax);
            CHECK(!tillFailed(u3) && fabs(u3 - ref3) <= 1e-6*fabs(ref3), "tillCalcU compression to 10 rhomax, v= %g: %.10e vs %.10e", v[k], u3, ref3);
            CHECK(fabs(back - u1) <= 1e-7*fabs(u1) + 1e-6*material->vmin, "tillCalcU round trip rhomax -> 10 rhomax -> rhomax, v= %g: %.10e vs %.10e", v[k], back, u1);
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

    /* Below the cold curve the lookup follows the coldest isentrope: u2 = u_cold(rho2), counted. */
    CHECK(tillLookupU(material, 2.0*material->rho0, 0.5*ucold, 3.0*material->rho0, 0) == tillColdULookup(material, 3.0*material->rho0), "lookup below cold curve clamps to u_cold(rho2) (compression)");
    CHECK(tillLookupU(material, 2.0*material->rho0, 0.5*ucold, 0.5*material->rho0, 0) == 0.0, "lookup below cold curve clamps to u_cold(rho2) = 0 (expansion)");
    CHECK(tillLookupU(material, 2.0*material->rho0, 0.5*ucold, 0.5*material->rhomin, 0) == 0.0, "lookup below cold curve, rho2 below rhomin");
    {
        /* rho2 above the table: the cold curve continued by integration from rhomax. */
        double rho2 = 2.0*material->rhomax;
        double ref = tillCalcU(material, material->rhomax, tillColdULookup(material, material->rhomax), rho2);
        double lk = tillLookupU(material, 2.0*material->rho0, 0.5*ucold, rho2, 0);
        CHECK(isfinite(lk) && lk == ref && lk > tillColdULookup(material, material->rhomax), "lookup below cold curve, rho2 above rhomax: %g vs %g", lk, ref);
    }
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
    /*
     * Within roundoff below the cold curve the state is on it: not below, in the table, T = 0 without
     * a fallback. Clearly below (10 tolerances) it is below, outside and a counted fallback. The same
     * below rhomin, where the cold curve is 0. Up to 4.0.0 tillIsBelowColdCurve used u < u_cold.
     */
    {
        const long nf = material->nLookupFail;
        const double rhos[2] = { 2.0*material->rho0, 0.5*material->rhomin };
        for (int q = 0; q < 2; q++) {
            const double uc = tillColdULookup(material, rhos[q]);
            const double uOn = uc - 0.1*tillUTol(material, uc), uOff = uc - 10.0*tillUTol(material, uc);
            CHECK(tillIsBelowColdCurve(material, rhos[q], uOn) == 0, "isbelow within roundoff of the cold curve (rho= %g)", rhos[q]);
            CHECK(tillTempRhoU(material, rhos[q], uOn) == 0.0, "T within roundoff of the cold curve (rho= %g)", rhos[q]);
            CHECK(material->nLookupFail == nf, "state within roundoff of the cold curve counted (rho= %g)", rhos[q]);
            CHECK(tillIsBelowColdCurve(material, rhos[q], uOff) == 1, "isbelow 10 tolerances below the cold curve (rho= %g)", rhos[q]);
            if (q == 0) {
                /* tillIsInTable is defined in [rhomin, rhomax] only. */
                CHECK(tillIsInTable(material, rhos[q], uOn) == TILL_LOOKUP_SUCCESS, "in table within roundoff of the cold curve");
                CHECK(tillIsInTable(material, rhos[q], uOff) == TILL_LOOKUP_OUTSIDE_VMIN, "outside 10 tolerances below the cold curve");
            }
        }
        CHECK(material->nLookupFail == nf, "isbelow and isintable counted (%ld)", material->nLookupFail - nf);
    }
    /* A non-finite u is rejected by every function that takes one, never passed on as NaN. */
    CHECK(tillIsBelowColdCurve(material, 2.0*material->rho0, NAN) == TILL_FAIL, "isbelow NaN u");
    /*
     * Below rhomin the cold curve is zero (tillColdULookup), so the status is u < 0 there, as in
     * EOSisbelowColdCurve; it fails only where tillColdULookup fails (up to 4.0.0 it failed for every
     * rho below rhomin). The failures are counted once each; nothing else is counted.
     */
    {
        const long nf = material->nLookupFail;
        CHECK(tillIsBelowColdCurve(material, 0.5*material->rhomin, -1.0) == 1, "isbelow below rhomin, u < 0");
        CHECK(tillIsBelowColdCurve(material, 0.5*material->rhomin, 0.0) == 0, "isbelow below rhomin, u = 0");
        CHECK(tillIsBelowColdCurve(material, 0.5*material->rhomin, 1.0) == 0, "isbelow below rhomin, u > 0");
        /* Below rho0 the cold curve is 0: below it by more than the roundoff tolerance is below, -1e-30 is on it. */
        CHECK(tillIsBelowColdCurve(material, 0.5*material->rho0, -10.0*tillUTol(material, 0.0)) == 1, "isbelow below rho0, u < 0");
        CHECK(tillIsBelowColdCurve(material, 0.5*material->rho0, -1e-30) == 0, "isbelow below rho0, u = -1e-30 is on the cold curve");
        CHECK(material->nLookupFail == nf, "isbelow below rhomin counted %ld times", material->nLookupFail - nf);
        CHECK(tillIsBelowColdCurve(material, 0.0, 1.0) == TILL_FAIL, "isbelow at rho = 0");
        CHECK(tillIsBelowColdCurve(material, -1.0, 1.0) == TILL_FAIL, "isbelow at rho < 0");
        CHECK(tillIsBelowColdCurve(material, NAN, 1.0) == TILL_FAIL, "isbelow at a NaN rho");
        CHECK(tillIsBelowColdCurve(material, 2.0*material->rhomax, NAN) == TILL_FAIL, "isbelow above rhomax, NaN u");
        CHECK(material->nLookupFail == nf + 4, "isbelow failures counted %ld times (expected 4)", material->nLookupFail - nf);
        material->nLookupFail = nf;
    }
    CHECK(tillFindEntropyCurve(material, 2.0*material->rho0, NAN, 0) == TILL_FAIL_VALUE, "entropy curve NaN u");
    CHECK(tillLookupU(material, 2.0*material->rho0, NAN, material->rho0, 0) == TILL_FAIL_VALUE, "lookup NaN u");
    CHECK(tillLookupU(material, 2.0*material->rho0, INFINITY, material->rho0, 0) == TILL_FAIL_VALUE, "lookup +inf u");
    CHECK(material->nLookupFail == nfail0 + 16, "every failure and fallback is counted (%ld)", material->nLookupFail - nfail0);
    /* Outside the table in rho the lookup falls back to direct integration: finite, and counted as a fallback. */
    CHECK(isfinite(tillLookupU(material, 0.5*material->rhomin, 1.0, material->rhomin, 0)), "lookup below rhomin integrates");
    CHECK(material->nLookupFail == nfail0 + 17, "direct integration is counted as a fallback (%ld)", material->nLookupFail - nfail0);
    CHECK(isfinite(tillLookupU(material, material->rho0, 1.0, 2.0*material->rhomax, 0)), "lookup above rhomax integrates");
    CHECK(material->nLookupFail == nfail0 + 18, "direct integration is counted as a fallback (%ld)", material->nLookupFail - nfail0);
    /* An invalid end point is a failure in every path (also the clamp below the cold curve), counted once. */
    CHECK(tillLookupU(material, 2.0*material->rho0, 2.0*ucold, NAN, 0) == TILL_FAIL_VALUE, "lookup to NaN rho2");
    CHECK(tillLookupU(material, 2.0*material->rho0, 2.0*ucold, INFINITY, 0) == TILL_FAIL_VALUE, "lookup to +inf rho2");
    CHECK(tillLookupU(material, 2.0*material->rho0, 2.0*ucold, 0.0, 0) == TILL_FAIL_VALUE, "lookup to rho2 = 0");
    CHECK(tillLookupU(material, 2.0*material->rho0, 0.5*ucold, INFINITY, 0) == TILL_FAIL_VALUE, "lookup below cold curve to +inf rho2");
    CHECK(tillLookupU(material, 2.0*material->rho0, 0.5*ucold, -1.0, 0) == TILL_FAIL_VALUE, "lookup below cold curve to rho2 < 0");
    CHECK(tillLookupU(material, NAN, 1.0, material->rho0, 0) == TILL_FAIL_VALUE, "lookup from NaN rho1");
    CHECK(tillLookupU(material, INFINITY, 1.0, material->rho0, 0) == TILL_FAIL_VALUE, "lookup from +inf rho1");
    CHECK(material->nLookupFail == nfail0 + 25, "invalid rho counted once each (%ld)", material->nLookupFail - nfail0);
    /* tillCalcU: invalid input fails (no undefined behaviour in the step count), counted; a valid integration is not. */
    CHECK(tillCalcU(material, NAN, 1.0, material->rho0) == TILL_FAIL_VALUE, "tillCalcU NaN rho1");
    CHECK(tillCalcU(material, material->rho0, 1.0, NAN) == TILL_FAIL_VALUE, "tillCalcU NaN rho2");
    CHECK(tillCalcU(material, INFINITY, 1.0, material->rho0) == TILL_FAIL_VALUE, "tillCalcU +inf rho1");
    CHECK(tillCalcU(material, material->rho0, 1.0, INFINITY) == TILL_FAIL_VALUE, "tillCalcU +inf rho2");
    CHECK(tillCalcU(material, material->rho0, 1.0, 0.0) == TILL_FAIL_VALUE, "tillCalcU rho2 = 0");
    CHECK(tillCalcU(material, -1.0, 1.0, material->rho0) == TILL_FAIL_VALUE, "tillCalcU rho1 < 0");
    CHECK(tillCalcU(material, material->rho0, NAN, 2.0*material->rho0) == TILL_FAIL_VALUE, "tillCalcU NaN u1");
    CHECK(tillCalcU(material, material->rho0, INFINITY, 2.0*material->rho0) == TILL_FAIL_VALUE, "tillCalcU +inf u1");
    CHECK(tillCalcU(material, material->rho0, NAN, material->rho0) == TILL_FAIL_VALUE, "tillCalcU NaN u1, rho2 = rho1");
    CHECK(material->nLookupFail == nfail0 + 34, "tillCalcU failures counted (%ld)", material->nLookupFail - nfail0);
    CHECK(tillCalcU(material, material->rho0, 1.0, material->rho0) == 1.0, "tillCalcU identity");
    CHECK(isfinite(tillCalcU(material, material->rho0, 1.0, 2.0*material->rho0)), "tillCalcU valid");
    /* A state below the cold curve (u < 0 below rho0) is lifted to it by the clamp. */
    CHECK(tillCalcU(material, 0.5*material->rho0, -material->vmin, 0.25*material->rho0) == 0.0, "tillCalcU clamps u >= 0");
    CHECK(material->nLookupFail == nfail0 + 34, "a valid tillCalcU is not counted (%ld)", material->nLookupFail - nfail0);
    /* rho <= 0 is invalid input for the cold curve (it returned 0, the expanded cold curve). */
    CHECK(tillColdULookup(material, 0.0) == TILL_FAIL_VALUE, "cold curve at rho = 0");
    CHECK(tillColdULookup(material, -1.0) == TILL_FAIL_VALUE, "cold curve at rho < 0");
    CHECK(tillURhoTemp(material, -1.0, 300.0) == TILL_FAIL_VALUE, "u(rho < 0, T)");
    CHECK(tillTempRhoU(material, 0.0, 1.0) == 0.0, "T(rho = 0, u) falls back to T = 0");
    CHECK(tillColdULookup(material, 1e-3*material->rhomin) == 0.0, "cold curve below rhomin (rho > 0) is still zero");
    CHECK(material->nLookupFail == nfail0 + 38, "rho <= 0 counted (%ld)", material->nLookupFail - nfail0);
    /* T(rho, u) below the cold curve is T = 0 (the lowest temperature), counted, never negative. */
    CHECK(tillTempRhoU(material, 2.0*material->rho0, 0.5*ucold) == 0.0, "T below the cold curve (compressed) is 0");
    CHECK(tillTempRhoU(material, 0.5*material->rho0, -material->vmin) == 0.0, "T for u < 0 (expanded) is 0");
    CHECK(tillTempRhoU(material, 2.0*material->rho0, -1e3*ucold) == 0.0, "T far below the cold curve is 0");
    CHECK(tillTempRhoU(material, 2.0*material->rho0, NAN) == 0.0, "T for a NaN u is 0");
    CHECK(tillTempRhoU(material, 2.0*material->rho0, INFINITY) == 0.0, "T for an infinite u is 0");
    CHECK(material->nLookupFail == nfail0 + 43, "T below the cold curve counted (%ld)", material->nLookupFail - nfail0);
    /* On the cold curve and within roundoff below it: T = 0 exactly, not counted; just above it: T > 0. */
    CHECK(tillTempRhoU(material, 2.0*material->rho0, ucold) == 0.0, "T on the cold curve is 0");
    CHECK(tillTempRhoU(material, 2.0*material->rho0, ucold - 0.1*tillUTol(material, ucold)) == 0.0, "T within roundoff below the cold curve is 0");
    CHECK(tillTempRhoU(material, 0.5*material->rho0, 0.0) == 0.0, "T on the expanded cold curve is 0");
    CHECK(tillTempRhoU(material, 2.0*material->rho0, ucold + material->cv) > 0.0, "T just above the cold curve");
    CHECK(material->nLookupFail == nfail0 + 43, "T on the cold curve not counted (%ld)", material->nLookupFail - nfail0);
    /* u(rho, T) and rho(P, T) below the lowest temperature T = 0 (or for a non-finite T): the sentinel, counted once each. */
    {
        const double Ts[4] = { -1.0, -1e-300, NAN, INFINITY };
        const long nfT = material->nLookupFail;
        for (int q = 0; q < 4; q++) {
            long nf = material->nLookupFail;
            CHECK(tillURhoTemp(material, 2.0*material->rho0, Ts[q]) == TILL_FAIL_VALUE, "u(rho, T= %g)", Ts[q]);
            CHECK(material->nLookupFail == nf + 1, "u(rho, T= %g) counted %ld times", Ts[q], material->nLookupFail - nf);
            nf = material->nLookupFail;
            CHECK(tillRhoPTemp(material, material->A, Ts[q]) == TILL_FAIL_VALUE, "rho(P, T= %g)", Ts[q]);
            CHECK(material->nLookupFail == nf + 1, "rho(P, T= %g) counted %ld times", Ts[q], material->nLookupFail - nf);
        }
        /* T = 0 is the cold curve: valid, not counted. */
        const long nf = material->nLookupFail;
        CHECK(tillURhoTemp(material, 2.0*material->rho0, 0.0) == tillColdULookup(material, 2.0*material->rho0), "u(rho, T = 0) is the cold curve");
        CHECK(tillURhoTemp(material, 0.5*material->rho0, 0.0) == tillColdULookup(material, 0.5*material->rho0), "u(rho < rho0, T = 0) is the cold curve");
        CHECK(material->nLookupFail == nf, "u(rho, T = 0) counted");
        nfail0 += material->nLookupFail - nfT;   /* keep the running count of the checks below */
    }
    /* Just below the cold curve (roundoff in the hydro code): clamped, close to the in-table result. */
    {
        double rho1 = 2.0*material->rho0, rho2 = 1.5*material->rho0;
        double ulo = ucold*(1.0 - 1e-6);
        double lk = tillLookupU(material, rho1, ulo, rho2, 0);
        CHECK(lk == tillColdULookup(material, rho2), "lookup just below cold curve clamps");
        CHECK_REL(lk + 1.0, tillLookupU(material, rho1, ucold, rho2, 0) + 1.0, 1e-6, "clamp continuous with the cold isentrope");
    }
}

/*
 * Initialise a material with the unit pair (dKpcUnit, dMsolUnit) in a child process. Returns 0 if
 * the library accepted it, 1 if it stopped the program (assert(0), or exit(1) with -DNDEBUG), -1
 * for anything else.
 */
static int initUnitsRejected(int iMat, double dKpcUnit, double dMsolUnit)
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);   /* the expected diagnostic */
        TILLMATERIAL *m = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
        tillFinalizeMaterial(m);
        _exit(0);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    if ((WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) || (WIFEXITED(status) && WEXITSTATUS(status) == 1)) return 1;
    return -1;
}

/*
 * Call one of the lookups that need the full isentrope table in a child process on a material that
 * has only the cold curve. Returns 1 if it stopped the program (assert(0), or exit(1) with
 * -DNDEBUG), 0 if it returned, -1 otherwise.
 */
static int tableLookupFatal(TILLMATERIAL *m, int which)
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);   /* the expected diagnostic */
        const double rho = 1.5*m->rho0, u = tillURhoTemp(m, rho, 1000.0);
        switch (which) {
            case 0: (void) tillLookupU(m, rho, u, 2.0*rho, 0); break;
            case 1: (void) tillLookupU(m, 1e-9*m->rho0, u, 2.0*rho, 0); break;   /* a path that would integrate directly */
            case 2: (void) tillFindEntropyCurve(m, rho, u, 0); break;
            case 3: (void) tillFindUonIsentrope(m, 1.0, rho); break;
            case 4: (void) tillCubicInt(m, rho, 1.0); break;
            case 5: { double e; int i; tillLocateRho(m, rho, &i, &e); (void) tillCubicIntGrid(m, i, e, m->nTableV/2); } break;
        }
        _exit(0);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    if ((WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) || (WIFEXITED(status) && WEXITSTATUS(status) == 1)) return 1;
    return -1;
}

/*
 * tillInitColdCurve builds the grid and the two edge isentropes only (what EOSinitMaterial does);
 * the full table follows with tillInitLookup (EOSinitIsentropicLookup). The edge curves must be
 * rows 0 and nTableV-1 of the full table bit for bit, everything that needs only them (temperature,
 * u(rho, T), rho(P, T), cold curve, table status) must be bit-identical with and without the table,
 * and the isentropic lookups must be fatal without it. tillInitColdCurve after tillInitLookup frees
 * the table.
 */
static void testColdCurveOnly(int iMat)
{
    const double dKpcUnit = 6.371e8/3.085678e21, dMsolUnit = 5.972e27/1.99e33;
    TILLMATERIAL *c = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    TILLMATERIAL *t = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
    /* Own random states: the later tests keep the sequence they were written with. */
    const unsigned long long rngSaved = rngState;
    int k, nDiff = 0;

    tillInitColdCurveDefault(c);
    tillInitLookupDefault(t);
    CHECK(c->Lookup == NULL && c->ColdCurve != NULL && c->TopCurve != NULL, "tillInitColdCurve builds the edge curves only");
    CHECK(c->nTableRho == t->nTableRho && c->nTableV == t->nTableV && c->n == t->n && c->dlogrho == t->dlogrho && c->dv == t->dv &&
          c->rhomin == t->rhomin && c->rhomax == t->rhomax && c->vmin == t->vmin && c->vmax == t->vmax &&
          c->logrhomin == t->logrhomin && c->invdlogrho == t->invdlogrho && c->invdv == t->invdv, "same grid as the full table");
    for (k = 0; k < t->nTableRho; k++) {
        const TILL_LOOKUP_ENTRY *r0 = &t->Lookup[TILL_INDEX(t, k, 0)], *rn = &t->Lookup[TILL_INDEX(t, k, t->nTableV-1)];
        if (c->ColdCurve[k].u != r0->u || c->ColdCurve[k].u1 != r0->u1 || c->TopCurve[k].u != rn->u || c->TopCurve[k].u1 != rn->u1) nDiff++;
    }
    CHECK(nDiff == 0, "edge curves differ from rows 0 and nTableV-1 of the table at %d densities", nDiff);

    /* Same results and the same failure counts without and with the table. */
    nDiff = 0;
    for (k = 0; k < 20000; k++) {
        double rho = lograndom(1e-6*c->rho0, 300.0*c->rho0), T = lograndom(1.0, 1e8);
        double u = (k % 5 == 0) ? -lograndom(1e-6, 1.0)*c->cv*T : tillURhoTemp(t, rho, T);
        if (tillFailed(u)) u = c->cv*T;
        long nc = c->nLookupFail, nt = t->nLookupFail;
        double a1 = tillTempRhoU(c, rho, u), b1 = tillTempRhoU(t, rho, u);
        double a2 = tillURhoTemp(c, rho, T), b2 = tillURhoTemp(t, rho, T);
        double a3 = tillColdULookup(c, rho), b3 = tillColdULookup(t, rho);
        double a4 = tillRhoPTemp(c, 1e-3*c->A*(1 + k % 97), T), b4 = tillRhoPTemp(t, 1e-3*t->A*(1 + k % 97), T);
        int a5 = tillIsInTable(c, rho, u), b5 = tillIsInTable(t, rho, u);
        int a6 = tillIsBelowColdCurve(c, rho, u), b6 = tillIsBelowColdCurve(t, rho, u);
#define SAME(x, y) ((x) == (y) || (tillFailed(x) && tillFailed(y)))
        if (!SAME(a1, b1) || !SAME(a2, b2) || !SAME(a3, b3) || !SAME(a4, b4) || a5 != b5 || a6 != b6 ||
            c->nLookupFail - nc != t->nLookupFail - nt) nDiff++;
#undef SAME
    }
    CHECK(nDiff == 0, "%d of 20000 states differ between the cold curve only and the full table", nDiff);

    /* The isentropic lookups need the table: fatal without it, on every path. */
    for (k = 0; k <= 5; k++) CHECK(tableLookupFatal(c, k) == 1, "lookup %d without the table must be fatal", k);
    CHECK(tableLookupFatal(t, 0) == 0 && tableLookupFatal(t, 2) == 0, "with the table the lookups return");

    /* Building the table later on the cold-curve material gives the same table; tillInitColdCurve frees it again. */
    tillInitLookupDefault(c);
    nDiff = 0;
    for (k = 0; k < c->nTableRho*c->nTableV; k++)
        if (c->Lookup[k].u != t->Lookup[k].u || c->Lookup[k].u1 != t->Lookup[k].u1 || c->Lookup[k].udv2 != t->Lookup[k].udv2 || c->Lookup[k].u1dv2 != t->Lookup[k].u1dv2) nDiff++;
    CHECK(nDiff == 0, "table built after tillInitColdCurve differs at %d entries", nDiff);
    tillInitColdCurveDefault(c);
    CHECK(c->Lookup == NULL && c->ColdCurve != NULL, "tillInitColdCurve frees the table");

    tillFinalizeMaterial(c);
    tillFinalizeMaterial(t);
    rngState = rngSaved;
}

/*
 * One rule for the unit pair in EOSlib and every sublibrary: both units > 0 (and finite) are code
 * units, both <= 0 is cgs, anything else is fatal. (Before, exactly one unit <= 0 made igeos divide
 * by zero, Tillotson assert, and ANEOS, SCvH and REOS3 silently use cgs.)
 */
static void testUnitPair(int iMat)
{
    MatName = "units";
    static const struct { double dKpcUnit, dMsolUnit; int bRejected; const char *what; } cases[] = {
        {2.06701e-13, 4.80438e-08, 0, "both > 0: code units"},
        {0.0, 0.0, 0, "both 0: cgs"},
        {-1.0, -1.0, 0, "both < 0: cgs"},
        {0.0, -1.0, 0, "0 and < 0: cgs"},
        {2.06701e-13, 0.0, 1, "dMsolUnit = 0 only"},
        {0.0, 4.80438e-08, 1, "dKpcUnit = 0 only"},
        {-1.0, 4.80438e-08, 1, "dKpcUnit < 0 only"},
        {2.06701e-13, NAN, 1, "dMsolUnit NaN"},
        {INFINITY, 4.80438e-08, 1, "dKpcUnit inf"},
        {NAN, NAN, 1, "both NaN"},
    };
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        int r = initUnitsRejected(iMat, cases[k].dKpcUnit, cases[k].dMsolUnit);
        CHECK(r == cases[k].bRejected, "%s: %s (got %d)", cases[k].what, cases[k].bRejected ? "must be fatal" : "must be accepted", r);
    }

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
    testUnitPair(GRANITE);

    for (m = 0; m < (argc > 1 ? argc-1 : nMat); m++) {
        int iMat = (argc > 1) ? atoi(argv[m+1]) : allMat[m];
        TILLMATERIAL *material = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
        int nFail0 = nFail;

        tilliMatString(material, name);
        MatName = name;
        fprintf(stdout, "Material %d (%s)\n", iMat, name);

        testUnits(iMat);
        testUnitInvariance(iMat);
        testColdCurveExpanded(iMat);
        testPressure(material);
        testNegativeU(material);
        testColdCurveOnly(iMat);

        /* The standard table, as EOSlib initializes it. */
        tillInitLookupDefault(material);
        CHECK(material->nTableRho == TILL_TABLE_NRHO && material->nTableV == TILL_TABLE_NV, "standard table sizes");
        CHECK(material->vmin == material->cv*TILL_TABLE_TMIN && material->vmax == material->cv*TILL_TABLE_TMAX, "standard table v range");
        CHECK(material->rhomin == TILL_TABLE_RHOMIN_FAC*material->rho0, "standard table rhomin");
        testTable(material);
        testTableNodes(material);
        testIsentropes(material);
        testTemperature(material);
        testSentinels(material);
        testCalcUOutside(material);
        /* Grid constants derived by tillInitLookup for the lookups. */
        CHECK(material->logrhomin == log(material->rhomin), "logrhomin");
        CHECK(material->invdlogrho == 1.0/material->dlogrho && material->invdv == 1.0/material->dv, "invdlogrho, invdv");

        /* Re-initialization with other parameters works (no leak, no stale state). */
        tillInitLookup(material, 300, 200, 1e-3, 50.0, 1e-4, 100.0);
        CHECK(material->nTableRho == 300 && material->nTableV == 200, "re-init sizes");
        CHECK(material->rhomax >= 50.0*exp(-material->dlogrho) && material->rhomax <= 50.0*exp(material->dlogrho), "re-init rhomax= %g adjusted to the grid", material->rhomax);
        CHECK(material->logrhomin == log(material->rhomin) && material->invdlogrho == 1.0/material->dlogrho && material->invdv == 1.0/material->dv, "re-init grid constants");
        testTable(material);

        tillFinalizeMaterial(material);
        fprintf(stdout, "  %s: %d failures\n", name, nFail - nFail0);
    }

    fprintf(stdout, "%d checks, %d failures, %d skipped\n", nCheck, nFail, nSkip);
    return nFail;
}