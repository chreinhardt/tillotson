/*
 * Copyright (c) 2014-2018 Christian Reinhardt and Joachim Stadel.
 *
 * This file provides all the functions for the Tillotson EOS library.
 * The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 * and convenient to use equation of state that can describe matter over
 * a large range of pressures, densities and internal energies.
 *
 * tillotson.c:
 *
 * tillInitMaterial: set the material parameters and convert them to the given unit system.
 *
 * tillPressureSound: pressure and sound speed for a given density and internal energy.
 *
 * tillPressure, tillSoundSpeed: wrappers of tillPressureSound.
 *
 * tillRhoPTemp: invert P(rho, T) for rho.
 *
 * tillReadStrengthParameters, tillReadMeltCurve: optional material strength (see tillotson.h).
 *
 * tillFinalizeMaterial: free memory.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>   /* isspace, tolower in the strength readers */
#include <float.h>   /* DBL_MIN in tillPressureCore */
#include "tillotson.h"

/* Initialisation errors (defective strength files): fatal, also with -DNDEBUG. */
#define TILL_INIT_FATAL(...) do { fprintf(stderr, __VA_ARGS__); assert(0); exit(1); } while (0)

/*
 * Compile with -DTILL_PRESS_MELOSH to cut the cold term for rho/rho0 < 0.8 (Melosh 1989) instead
 * of where it is negative. The total pressure is cut at zero in both modes, so the flag does not
 * allow tension; it removes the (positive) cold term of expanded cold states below 0.8 rho0, which
 * makes P(rho) jump there.
 */

/*
 * Initialize a material from the Tillotson library
 *
 * We do:
 * Initialize variables
 * Convert quantities to code units
 * The memory for the look up table is allocated in tillInitLookup() (the cold curve and the top
 * isentrope alone in tillInitColdCurve())
 * Read the optional strength files <strengthBase>_strength and <strengthBase>_tmelt (none if
 * strengthBase is NULL)
 */
TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit)
{
    char base[256], name[64];
    TILLMATERIAL tmp;
    size_t k;

    /* The lower-case material name; an unknown material is rejected by tillInitMaterialFromFile. */
    tmp.iMaterial = iMaterial;
    tilliMatString(&tmp, name);
    for (k = 0; name[k] != '\0'; k++) name[k] = (char) tolower((unsigned char) name[k]);
    snprintf(base, sizeof(base), "%s%s.in", TILL_STRENGTH_FILE_PREFIX, name);

    return tillInitMaterialFromFile(iMaterial, base, dKpcUnit, dMsolUnit);
}

TILLMATERIAL *tillInitMaterialFromFile(int iMaterial, const char *strengthBase, double dKpcUnit, double dMsolUnit)
{
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    TILLMATERIAL *material;
	 
    material = (TILLMATERIAL *) calloc(1, sizeof(TILLMATERIAL));
    assert(material != NULL);

	material->iMaterial = iMaterial;

	/* This two parameters define the unit system we use. */
    material->dKpcUnit = dKpcUnit;
    material->dMsolUnit = dMsolUnit;

    material->rhomin = 0.0;
	material->rhomax = 0.0;
	material->vmax = 0.0;

    material->n = 0;
    material->dlogrho = 0.0;
    material->dv = 0.0;
    material->logrhomin = 0.0;
    material->invdlogrho = 0.0;
    material->invdv = 0.0;

	/*
	 * Unit system (the same rule in EOSlib and every sublibrary): both units > 0 and finite are code
	 * units, both <= 0 is cgs; anything else (exactly one of them <= 0, NaN, inf) is a caller error.
	 */
	if (!(((dKpcUnit > 0.0) && (dMsolUnit > 0.0) && isfinite(dKpcUnit) && isfinite(dMsolUnit)) || ((dKpcUnit <= 0.0) && (dMsolUnit <= 0.0)))) {
		fprintf(stderr, "tillInitMaterial: invalid unit system dKpcUnit= %g dMsolUnit= %g (give both > 0 for code units or both <= 0 for cgs).\n", dKpcUnit, dMsolUnit);
		assert(0);
		exit(1);
	}

	if (dKpcUnit <= 0.0)
	{
		/* Both <= 0 (checked above): units are not converted, so the code units are cgs. */
		material->dErgPerGmUnit = 1.0;
		material->dGmPerCcUnit = 1.0;
		material->dSecUnit = 1.0;
	} else {
		/* code energy per unit mass --> erg per g */
		material->dErgPerGmUnit = GCGS*material->dMsolUnit*MSOLG/(material->dKpcUnit*KPCCM);
		/* code density --> g per cc */
		material->dGmPerCcUnit = (material->dMsolUnit*MSOLG)/pow(material->dKpcUnit*KPCCM,3.0);
		/* code time --> seconds */
		material->dSecUnit = sqrt(1/(material->dGmPerCcUnit*GCGS));
	}

	/*
	 * k_B/m_H in code units (specific energy per K), as in igeos. It was k_B (not divided by m_H)
	 * in cgs; the code-unit value is unchanged (dKpcUnit*KPCCM/(GCGS*dMsolUnit*MSOLG) is
	 * 1/dErgPerGmUnit).
	 */
	material->dGasConst = (KBOLTZ/MHYDR)/material->dErgPerGmUnit;

	/*
	 * Set the Tillotson parameters for the material.
	 */
	switch(iMaterial)
	{
		case GRANITE:
			/*
			 * Material parameters from Benz et al. (1986).
			 */
			material->a = 0.5;
			material->b = 1.3;
			material->u0 = 1.6e11;		/* in ergs/g */
			material->rho0 = 2.7;		/* g/cc */
			material->A = 1.8e11;		/* ergs/cc */
			material->B = 1.8e11;		/* ergs/cc */
			material->us = 3.5e10;		/* ergs/g */
			material->us2 = 1.8e11;		/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.79e7;		/* ergs/g K (or 790 J/kg K) */ 
			break;
		case IRON:
			/*
			 * Material parameters from Benz et al. (1987).
			 */
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 9.5e10;		/* in ergs/g */
			material->rho0 = 7.86;		/* g/cc */
			material->A = 1.279e12;		/* ergs/cc */
			material->B = 1.05e12;		/* ergs/cc */
			material->us = 1.42e10;		/* ergs/g */
			material->us2 = 8.45e10;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.449e7;		/* ergs/g K */ 
			break;
		case BASALT:
			/*
			 ** Material parameters from Benz & Asphaug (1999).
			 */
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 4.87e12;		/* in ergs/g */
			material->rho0 = 2.7;		/* g/cc */
			material->A = 2.67e11;		/* ergs/cc */
			material->B = 2.67e11;		/* ergs/cc */
			material->us = 4.72e10;		/* ergs/g */
			material->us2 = 1.82e11;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.84e7;		/* ergs/g K */ 
			break;
		case ICE:
			/*
			 * Material parameters from Benz & Asphaug (1999).
			 */
			material->a = 0.3;
			material->b = 0.1;
			material->u0 = 1.0e11;		/* in ergs/g */
			material->rho0 = 0.917;		/* g/cc */
			material->A = 9.47e10;		/* ergs/cc */
			material->B = 9.47e10;		/* ergs/cc */
			material->us = 7.73e9;		/* ergs/g */
			material->us2 = 3.04e10;	/* ergs/g */
			/*
			 * Note that the expression for Pe provided in Benz et al. (1986) does not agree
			 * in alpha and beta with the one provided in Melosh (1989) or Benz & Asphaug (1999).
             * It seems that alpha and beta were switched so we use alpha=5 and beta=10.
			 */
			material->alpha = 5.0;
			material->beta = 10.0;
			// Have to look up in more details?
			material->cv = 1.0e7;		/* ergs/g K */ 
			break;
		case WATER:
			/*
			 * Material parameters from Woolfson (2007).
			 */
			material->a = 0.5;
			material->b = 0.9;
			material->u0 = 2.0e10;		/* in ergs/g */
			material->rho0 = 1.00;		/* g/cc */
			material->A = 2.00e11;		/* ergs/cc */
			material->B = 1.00e11;		/* ergs/cc */
			material->us = 4.00e9;		/* ergs/g */
			material->us2 = 2.04e10;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			// Have to look up in more details?
			material->cv = 4.1814e7;	/* ergs/g K */ 
			break;
		case DUNITE:
			/*
			 * Material parameters that imitate dunite in ANEOS (Thompson 1972).
			 */
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 4.87e12;		/* in ergs/g */
            // Key difference is the reference density
			material->rho0 = 3.32;		/* g/cc */
			material->A = 2.67e11;		/* ergs/cc */
			material->B = 2.67e11;		/* ergs/cc */
			material->us = 4.72e10;		/* ergs/g */
			material->us2 = 1.82e11;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.84e7;		/* ergs/g K */ 
			break;
		case FORSTERITE:
			/*
			 * Material parameters for Forsterite from S. Raducan.
			 */
			material->a = 0.1;
			material->b = 1.0;
			material->u0 = 4.87e12;		/* in ergs/g */
            // Key difference is the reference density
			material->rho0 = 3.3;		/* g/cc */
			material->A = 9.5e11;		/* ergs/cc */
			material->B = 9.5e11;		/* ergs/cc */
			material->us = 4.72e10;		/* ergs/g */
			material->us2 = 1.82e11;	/* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 1e7;		/* ergs/g K */ 
			break;
		default:
			fprintf(stderr, "tillInitMaterial: unknown material %i.\n", iMaterial);
			assert(0);
			exit(1);
	}

    /*
     * Convert energies and densities to code units!
     */
    material->u0 /= material->dErgPerGmUnit;
    material->us /= material->dErgPerGmUnit;
    material->us2 /= material->dErgPerGmUnit;
    material->rho0 /= material->dGmPerCcUnit;
    material->A /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
    material->B /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
 
    material->cv /= material->dErgPerGmUnit;

    /*
     * Optional material strength (see "Material strength" in tillotson.h): no strength file means
     * no yield stress (model -1). Models 1 and 2 soften the strength with T/T_melt, so a strength
     * file that asks for them without a melt curve is a misconfigured material and stops the run
     * here, as in ANEOSmaterial. Model 0 (no limiter) needs no melt curve.
     */
    material->yieldStrengthModel = -1;
    material->mud = -1.0;
    if (strengthBase != NULL) {
        char file[512];

        snprintf(file, sizeof(file), "%s_tmelt", strengthBase);
        (void) tillReadMeltCurve(material, file);

        snprintf(file, sizeof(file), "%s_strength", strengthBase);
        (void) tillReadStrengthParameters(material, file);

        if (material->yieldStrengthModel >= 1 && material->TMelt == NULL)
            TILL_INIT_FATAL("tillInitMaterial: %s_strength sets yieldStrengthModel = %d, which needs the melt curve "
                            "%s_tmelt, but that file was not found.\n", strengthBase, material->yieldStrengthModel, strengthBase);
    }

    return material;
}

/*
 * Free the memory.
 */
void tillFinalizeMaterial(TILLMATERIAL *material)
{

    if (material->Lookup != NULL) free(material->Lookup);
    if (material->ColdCurve != NULL) free(material->ColdCurve);
    if (material->TopCurve != NULL) free(material->TopCurve);
    if (material->rhoMelt != NULL) free(material->rhoMelt);
    if (material->TMelt != NULL) free(material->TMelt);

    free(material);
}

#define TILL_STRENGTH_LINE_MAX 256

/*
 * Read the strength parameters (text, the format of ANEOSReadStrengthParameters; see "Material
 * strength" in tillotson.h). Returns 1 if the file does not exist (model -1), 0 on success. A
 * malformed file or parameters the yield model cannot use are fatal, so that a bad file stops the
 * run at initialisation and not at the first strength evaluation:
 *   - every value finite (strtod accepts "nan" and "inf"), one value per line;
 *   - yieldStrengthModel 0, 1 or 2;
 *   - Y0 >= 0, YM >= 0, mui >= 0, xi >= 0, Gamma >= 0 (or the sentinel -1);
 *   - for models 1 and 2: YM > Y0 (otherwise the intact strength has a pole at
 *     P = (Y0 - YM)/mui), mui > 0, xi > 0 (xi = 0 means no strength at all) and a fixed Gamma > 0;
 *   - for model 2: mud >= 0.
 */
int tillReadStrengthParameters(TILLMATERIAL *material, const char *inputfile)
{
    FILE *file;
    char buffer[TILL_STRENGTH_LINE_MAX];
    char key[TILL_STRENGTH_LINE_MAX];
    char value[TILL_STRENGTH_LINE_MAX];
    int iLine = 0;

    if (material == NULL) TILL_INIT_FATAL("tillReadStrengthParameters: material not initialized.\n");

    if ((file = fopen(inputfile, "r")) == NULL) {
        material->yieldStrengthModel = -1;   /* no yield strength available */
        return 1;
    }

    material->yieldStrengthModel = -1;
    material->Y0 = -1.0;
    material->YM = -1.0;
    material->mui = -1.0;
    material->mud = -1.0;
    material->xi = -1.0;
    material->Gamma = -2.0;
    material->fixedGamma = 1;

    while (fgets(buffer, sizeof(buffer), file)) {
        int nUsed = 0, bLineComplete;
        char *hash, *end = NULL;

        iLine++;

        /*
         * Reject lines longer than the buffer. A line without '\n' is complete only at the end of
         * the file; feof is not yet set when the last line exactly fills the buffer, so the next
         * character decides.
         */
        bLineComplete = (strchr(buffer, '\n') != NULL) || feof(file);
        if (!bLineComplete) {
            int c = fgetc(file);
            if (c == EOF) bLineComplete = 1;
            else ungetc(c, file);
        }
        if (!bLineComplete) {
            fclose(file);
            TILL_INIT_FATAL("tillReadStrengthParameters: %s line %d exceeds %d characters.\n", inputfile, iLine, TILL_STRENGTH_LINE_MAX - 1);
        }

        hash = strchr(buffer, '#');
        if (hash) *hash = '\0';

        if (sscanf(buffer, " %255s", key) != 1) continue;   /* blank line */

        /* key = value; nothing but whitespace may follow the value. */
        if (sscanf(buffer, " %255[^= \t\r\n] = %255s%n", key, value, &nUsed) != 2) {
            fclose(file);
            TILL_INIT_FATAL("tillReadStrengthParameters: cannot parse %s line %d: '%s'\n", inputfile, iLine, buffer);
        }
        {
            const char *rest = buffer + nUsed;
            while (*rest != '\0' && isspace((unsigned char) *rest)) rest++;
            if (*rest != '\0') {
                fclose(file);
                TILL_INIT_FATAL("tillReadStrengthParameters: %s line %d: unexpected '%s' after the value of '%s'.\n",
                                inputfile, iLine, rest, key);
            }
        }

        if (strcmp(key, "yieldStrengthModel") == 0) {
            material->yieldStrengthModel = (int) strtol(value, &end, 10);
        } else {
            double v = strtod(value, &end);
            if      (strcmp(key, "Y0")  == 0) material->Y0  = v;
            else if (strcmp(key, "YM")  == 0) material->YM  = v;
            else if (strcmp(key, "mui") == 0) material->mui = v;
            else if (strcmp(key, "mud") == 0) material->mud = v;
            else if (strcmp(key, "xi")  == 0) material->xi  = v;
            else if (strcmp(key, "Gamma") == 0) {
                if (v == -1.0) {
                    /* Sentinel: shear modulus from the bulk modulus (tillGammaofRhoT). */
                    material->fixedGamma = 0;
                    material->Gamma = 0.0;
                } else {
                    material->fixedGamma = 1;
                    material->Gamma = v;
                }
            } else {
                fclose(file);
                TILL_INIT_FATAL("tillReadStrengthParameters: unknown parameter '%s' in %s line %d\n", key, inputfile, iLine);
            }
        }

        if (end == value || *end != '\0') {
            fclose(file);
            TILL_INIT_FATAL("tillReadStrengthParameters: invalid value '%s' for '%s' in %s line %d\n", value, key, inputfile, iLine);
        }
    }

    fclose(file);

    /* Missing parameters keep their negative initial values. NaN fails every comparison. */
    if ((material->yieldStrengthModel < 0) || (material->yieldStrengthModel > 2) ||
        !(material->Y0 >= 0.0) || !(material->YM >= 0.0) || !(material->mui >= 0.0) ||
        !(material->xi >= 0.0) || !(material->Gamma >= 0.0) ||
        !isfinite(material->Y0) || !isfinite(material->YM) || !isfinite(material->mui) ||
        !isfinite(material->xi) || !isfinite(material->Gamma) || !isfinite(material->mud))
        TILL_INIT_FATAL("tillReadStrengthParameters: missing or invalid parameter in %s "
                        "(yieldStrengthModel = %d (0, 1 or 2) Y0 = %g YM = %g mui = %g xi = %g Gamma = %g mud = %g)\n",
                        inputfile, material->yieldStrengthModel, material->Y0, material->YM,
                        material->mui, material->xi, material->Gamma, material->mud);

    if (material->yieldStrengthModel >= 1) {
        if (!(material->YM > material->Y0) || !(material->mui > 0.0) || !(material->xi > 0.0) ||
            (material->fixedGamma && !(material->Gamma > 0.0)))
            TILL_INIT_FATAL("tillReadStrengthParameters: yieldStrengthModel = %d in %s needs YM > Y0, "
                            "mui > 0, xi > 0 and a positive fixed Gamma (Y0 = %g YM = %g mui = %g xi = %g Gamma = %g)\n",
                            material->yieldStrengthModel, inputfile, material->Y0, material->YM,
                            material->mui, material->xi, material->Gamma);
        if (material->yieldStrengthModel == 2 && !(material->mud >= 0.0))
            TILL_INIT_FATAL("tillReadStrengthParameters: yieldStrengthModel = 2 in %s needs mud >= 0 "
                            "(damaged material), got mud = %g (missing?)\n", inputfile, material->mud);
    }

    return 0;
}

/*
 * Read the melt curve (binary: int n, double rho[n] in g/cc, double T_melt[n] in K; see "Material
 * strength" in tillotson.h). Returns 1 if the file does not exist, 0 on success. A short read,
 * trailing data, n < 2, a density that is not positive, finite and strictly increasing, or a melt
 * temperature that is negative or not finite are fatal. Replaces a curve read before.
 */
int tillReadMeltCurve(TILLMATERIAL *material, const char *inputfile)
{
    FILE *file;
    double *rho, *T;
    int n, i;

    if (material == NULL) TILL_INIT_FATAL("tillReadMeltCurve: material not initialized.\n");

    if ((file = fopen(inputfile, "rb")) == NULL) return 1;   /* optional file */

    if (fread(&n, sizeof(n), 1, file) != 1) {
        fclose(file);
        TILL_INIT_FATAL("tillReadMeltCurve: failed to read the number of nodes from %s.\n", inputfile);
    }
    if (n < 2 || n > 100000000) {
        fclose(file);
        TILL_INIT_FATAL("tillReadMeltCurve: invalid number of nodes n = %d in %s (need n >= 2).\n", n, inputfile);
    }

    rho = (double *) malloc((size_t) n*sizeof(double));
    T = (double *) malloc((size_t) n*sizeof(double));
    if (rho == NULL || T == NULL) TILL_INIT_FATAL("tillReadMeltCurve: could not allocate %d nodes for %s.\n", n, inputfile);

    if (fread(rho, sizeof(double), (size_t) n, file) != (size_t) n || fread(T, sizeof(double), (size_t) n, file) != (size_t) n) {
        fclose(file);
        TILL_INIT_FATAL("tillReadMeltCurve: %s is too short for n = %d nodes (int n, double rho[n], double T_melt[n]).\n", inputfile, n);
    }
    /* Nothing may follow: catches a file in another format (e.g. an ANEOS _tmelt file, which has no density axis). */
    if (fgetc(file) != EOF) {
        fclose(file);
        TILL_INIT_FATAL("tillReadMeltCurve: %s has data after the %d nodes (int n, double rho[n], double T_melt[n]).\n", inputfile, n);
    }
    fclose(file);

    for (i = 0; i < n; i++) {
        if (!isfinite(rho[i]) || !(rho[i] > 0.0))
            TILL_INIT_FATAL("tillReadMeltCurve: rho[%d] = %g g/cc in %s is not positive and finite.\n", i, rho[i], inputfile);
        if (i > 0 && !(rho[i] > rho[i-1]))
            TILL_INIT_FATAL("tillReadMeltCurve: rho in %s is not strictly increasing at %d (%.15e <= %.15e).\n", inputfile, i, rho[i], rho[i-1]);
        /* 0 where there is no condensed phase; negative or non-finite is a defect. */
        if (!isfinite(T[i]) || !(T[i] >= 0.0))
            TILL_INIT_FATAL("tillReadMeltCurve: T_melt[%d] = %g K in %s is negative or not finite.\n", i, T[i], inputfile);
    }

    if (material->rhoMelt != NULL) free(material->rhoMelt);
    if (material->TMelt != NULL) free(material->TMelt);
    material->nMelt = n;
    material->rhoMelt = rho;
    material->TMelt = T;

    return 0;
}

/*
 * Return the yield strength model and its parameters in code units (any pointer may be NULL).
 */
int tillYieldParameters(TILLMATERIAL *material, double *Y0, double *YM, double *mui, double *mud, double *xi)
{
    const double dPUnit = material->dErgPerGmUnit*material->dGmPerCcUnit;   /* code pressure -> erg/cm^3 */

    if (Y0) *Y0 = material->Y0/dPUnit;
    if (YM) *YM = material->YM/dPUnit;
    if (mui) *mui = material->mui;
    if (mud) *mud = material->mud;
    if (xi) *xi = material->xi;
    return material->yieldStrengthModel;
}

/*
 * The melt temperature T_melt(rho) in K (rho in code units), linear in rho between the nodes and
 * constant at the end values outside them (conservative for the strength: T_melt increases with
 * compression, so clamping above the curve underestimates it). Fatal without a melt curve (a
 * strength model was requested without its melt curve); TILL_FAIL_VALUE (counted) for a rho that
 * is not positive and finite.
 */
double tillTmeltofRho(TILLMATERIAL *material, double rho)
{
    const double rho_cgs = rho*material->dGmPerCcUnit;
    const int n = material->nMelt;
    int lo = 0, hi = n - 1;

    if (material->TMelt == NULL) {
        fprintf(stderr, "tillTmeltofRho: the melt curve of Tillotson material %i was not loaded.\n", material->iMaterial);
        assert(0);
        exit(1);
    }

    if (!(rho_cgs > 0.0) || !isfinite(rho_cgs))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillTmeltofRho: invalid rho= %g (iMat= %i).\n", rho, material->iMaterial);

    if (rho_cgs <= material->rhoMelt[0]) return material->TMelt[0];
    if (rho_cgs >= material->rhoMelt[n-1]) return material->TMelt[n-1];

    /* Bisection: rhoMelt[lo] <= rho_cgs < rhoMelt[hi]. */
    while (hi - lo > 1) {
        int mid = lo + (hi - lo)/2;
        if (material->rhoMelt[mid] <= rho_cgs) lo = mid; else hi = mid;
    }
    {
        const double x = (rho_cgs - material->rhoMelt[lo])/(material->rhoMelt[hi] - material->rhoMelt[lo]);
        return material->TMelt[lo] + x*(material->TMelt[hi] - material->TMelt[lo]);
    }
}

/*
 * Shear modulus in code units: the fixed Gamma of the strength file, or (Gamma = -1 there) the
 * bulk modulus K = rho c^2 times 3(1 - 2 nu)/(2(1 + nu)) with a guessed Poisson ratio nu = 0.25
 * ("normal" materials have nu in [0.2, 0.3]; the shear modulus changes by less than 25%), as in
 * ANEOSGammaofRhoT. c = c(rho, u(rho, T)) of the cut EOS, floored at sqrt(A/rho0). Returns
 * TILL_FAIL_VALUE (counted once, by u(rho, T) or the pressure) where u(rho, T) fails (T < 0 or not
 * finite, rho above the table or invalid).
 */
double tillGammaofRhoT(TILLMATERIAL *material, double rho, double T)
{
    const double nu = 0.25;
    double u, c;

    if (material->fixedGamma) return material->Gamma/(material->dErgPerGmUnit*material->dGmPerCcUnit);

    u = tillURhoTemp(material, rho, T);
    if (tillFailed(u)) return TILL_FAIL_VALUE;
    tillPressureSound(material, rho, u, &c);
    if (tillFailed(c)) return TILL_FAIL_VALUE;
    return rho*c*c*3.0*(1.0 - 2.0*nu)/(2.0*(1.0 + nu));
}

/*
 * This function stores the materials name in a string string.
 */
void tilliMatString(TILLMATERIAL *material, char *MatName)
{
    assert(MatName != NULL);

    switch(material->iMaterial)
    {
        case GRANITE:
            sprintf(MatName, "GRANITE");
            break;
        case IRON:
            sprintf(MatName, "IRON");
            break;
        case BASALT:
            sprintf(MatName, "BASALT");
            break;
        case ICE:
            sprintf(MatName, "ICE");
            break;
        case WATER:
            sprintf(MatName, "WATER");
            break;
        case DUNITE:
            sprintf(MatName, "DUNITE");
            break;
        case FORSTERITE:
            sprintf(MatName, "FORSTERITE");
            break;
        default:
            sprintf(MatName, "UNKNOWN");
    }
}

/*
 * Name of a tillIsInTable status code (static string).
 */
const char *tillLookupStatusString(int iError)
{
    switch(iError)
    {
        case TILL_LOOKUP_SUCCESS:        return "TILL_LOOKUP_SUCCESS";
        case TILL_LOOKUP_OUTSIDE_RHOMIN: return "TILL_LOOKUP_OUTSIDE_RHOMIN";
        case TILL_LOOKUP_OUTSIDE_RHOMAX: return "TILL_LOOKUP_OUTSIDE_RHOMAX";
        case TILL_LOOKUP_OUTSIDE_VMIN:   return "TILL_LOOKUP_OUTSIDE_VMIN";
        case TILL_LOOKUP_OUTSIDE_VMAX:   return "TILL_LOOKUP_OUTSIDE_VMAX";
        default:                         return "UNKNOWN ERROR";
    }
}

/*
 * Same, copied into the caller's buffer (at least 32 characters).
 */
void tillErrorString(int iError, char *ErrorString)
{
    assert(ErrorString != NULL);
    strcpy(ErrorString, tillLookupStatusString(iError));
}

/*
 * This function prints the material constants for a given material.
 */
void tillPrintMat(TILLMATERIAL *material, FILE *fp)
{
    char MatName[256];

    assert(material != NULL);
    tilliMatString(material, MatName);

    /*
     * Dimensional values in cgs with the unit, like ANEOSPrintMat, scvheosPrintMat and reos3PrintMat
     * (they were printed in code units without a label). The material is stored in code units.
     */
    const double dPUnit = material->dErgPerGmUnit*material->dGmPerCcUnit;   /* code pressure -> erg/cm^3 */
    fprintf(fp,"# Material: %i (%s)\n", material->iMaterial, MatName);
    fprintf(fp,"# dKpcUnit: %g dMsolUnit: %g\n", material->dKpcUnit, material->dMsolUnit);
    fprintf(fp,"# a:            %g\n", material->a);
    fprintf(fp,"# b:            %g\n", material->b);
    fprintf(fp,"# A (cgs units): %g erg/cm^3\n", material->A*dPUnit);
    fprintf(fp,"# B (cgs units): %g erg/cm^3\n", material->B*dPUnit);
    fprintf(fp,"# Reference density rho0 (cgs units): %g g/cc\n", material->rho0*material->dGmPerCcUnit);
    fprintf(fp,"# u0 (cgs units): %g erg/g\n", material->u0*material->dErgPerGmUnit);
    fprintf(fp,"# us (cgs units): %g erg/g\n", material->us*material->dErgPerGmUnit);
    fprintf(fp,"# us2 (cgs units): %g erg/g\n", material->us2*material->dErgPerGmUnit);
    fprintf(fp,"# alpha:        %g\n", material->alpha);
    fprintf(fp,"# beta:         %g\n", material->beta);
    fprintf(fp,"# cv (cgs units): %g erg/g/K\n", material->cv*material->dErgPerGmUnit);
    if (material->Lookup != NULL) {
        fprintf(fp,"# Lookup table (cgs units): %i x %i, rho in [%g, %g] g/cc, v in [0, %g] erg/g (vmin= %g erg/g)\n",
                material->nTableRho, material->nTableV, material->rhomin*material->dGmPerCcUnit, material->rhomax*material->dGmPerCcUnit,
                material->vmax*material->dErgPerGmUnit, material->vmin*material->dErgPerGmUnit);
    } else if (material->ColdCurve != NULL) {
        fprintf(fp,"# Lookup table (cgs units): cold curve and top isentrope only (%i x 2), rho in [%g, %g] g/cc, v in [0, %g] erg/g (vmin= %g erg/g); isentrope table not built\n",
                material->nTableRho, material->rhomin*material->dGmPerCcUnit, material->rhomax*material->dGmPerCcUnit,
                material->vmax*material->dErgPerGmUnit, material->vmin*material->dErgPerGmUnit);
    } else {
        fprintf(fp,"# Lookup table: not initialized\n");
    }
    /* Strength (as read, cgs). */
    if (material->yieldStrengthModel < 0) {
        fprintf(fp,"# Yield strength model: %d (no strength file)\n", material->yieldStrengthModel);
    } else {
        fprintf(fp,"# Yield strength model: %d, mui= %g mud= %g xi= %g\n", material->yieldStrengthModel,
                material->mui, material->mud, material->xi);
        fprintf(fp,"# Yield strength Y0, YM (cgs units): %g, %g erg/cm^3\n", material->Y0, material->YM);
        if (material->fixedGamma) fprintf(fp,"# Shear modulus Gamma (cgs units): %g erg/cm^3\n", material->Gamma);
        else fprintf(fp,"# Shear modulus: from the bulk modulus (nu = 0.25)\n");
    }
    if (material->TMelt != NULL) {
        fprintf(fp,"# Melt curve (cgs units): %d nodes, rho in [%g, %g] g/cc, T_melt in [%g, %g] K at the ends\n",
                material->nMelt, material->rhoMelt[0], material->rhoMelt[material->nMelt-1],
                material->TMelt[0], material->TMelt[material->nMelt-1]);
    } else {
        fprintf(fp,"# Melt curve: none\n");
    }
    fprintf(fp,"# Failed lookups and fallbacks (nLookupFail): %ld\n", material->nLookupFail);
}

/*
 * The Tillotson EOS: pressure, sound speed and the derivatives dP/drho (u = const) and dP/du
 * (rho = const), evaluated consistently in one place.
 *
 * Regions (Benz & Asphaug 1999):
 *   condensed (rho >= rho0) and expanded cold (rho < rho0, u < us):   P = Pc
 *   expanded hot (rho < rho0, u > us2):                                P = Pe
 *   intermediate (rho < rho0, us <= u <= us2):                         P = (1-y) Pc + y Pe
 *
 * Pressure cutoff (bClip != 0, the library default): the cold term Pc is set to zero where it is
 * negative (expanded cold states in tension), and the final pressure is cut at zero as well. The
 * derivatives are those of the cut pressure, i.e. zero where the cutoff is active, and the sound
 * speed c^2 = (dP/drho)_u + P/rho^2 (dP/du)_rho is floored at the bulk sound speed sqrt(A/rho0).
 * With TILL_PRESS_MELOSH the cold term is cut for rho/rho0 < 0.8 as suggested by Melosh (1989)
 * instead, also where it is positive (P jumps at 0.8 rho0); negative values above 0.8 rho0 are
 * still removed by the final cut. In both modes the cut pressure is P = max(P_NP, 0), except that
 * with TILL_PRESS_MELOSH it is 0 for expanded cold states below 0.8 rho0.
 *
 * bClip == 0 (tillPressureNP, tillPressureSoundNP) leaves the cold region raw, including negative
 * pressures, and applies no final cut; this is what the inversion tillRhoPTemp needs because the
 * raw P(rho) at fixed T is monotonic where the cut one is flat. The cut of Pc inside the
 * intermediate blend is kept in both modes (see below).
 *
 * Any of pcSound, pdPdrho and pdPdu may be NULL. If all three are NULL (tillPressure, the RK4
 * stages of the isentrope integration) the derivatives are not computed at all; P is the same
 * expression either way, so it is bit-identical with and without them.
 */
static double tillPressureCore(TILLMATERIAL *material, double rho, double u, int bClip,
                               double *pcSound, double *pdPdrho, double *pdPdu)
{
    const double rho0 = material->rho0, A = material->A, B = material->B;
    const double a = material->a, b = material->b, u0 = material->u0;
    const double alpha = material->alpha, beta = material->beta;
    double eta, mu, z, w0, ex;
    double Pc = 0.0, Pe = 0.0, dPcdrho = 0.0, dPedrho = 0.0, dPcdu = 0.0, dPedu = 0.0;
    double P, c2, dPdrho = 0.0, dPdu = 0.0, y;
    int bCold, bHot;
    const int bDeriv = (pcSound != NULL || pdPdrho != NULL || pdPdu != NULL);

    /*
     * An invalid density (rho <= 0, NaN, inf) has no pressure: eta = 0 makes z infinite and P NaN.
     * Like a density outside a table it fails, counted once, with every output set to the
     * sentinel. Up to 4.0.0 it returned NaN (rho = 0, NaN) or inf, uncounted. One comparison on the
     * hot path; the isentrope integration never gets here with an invalid rho (rho = exp(logrho)).
     */
    if (!(rho > 0.0) || !isfinite(rho)) {
        material->nLookupFail++;
        TILL_WARN("tillPressureCore: invalid density rho= %g (u= %g, iMat= %i).\n", rho, u, material->iMaterial);
        if (pcSound != NULL) *pcSound = TILL_FAIL_VALUE;
        if (pdPdrho != NULL) *pdPdrho = TILL_FAIL_VALUE;
        if (pdPdu != NULL) *pdPdu = TILL_FAIL_VALUE;
        return TILL_FAIL_VALUE;
    }

    eta = rho/rho0;
    mu = eta - 1.0;
    z = (1.0 - eta)/eta;
    /*
     * w0 = u/(u0 eta^2) + 1 enters as b/w0 and is only defined for u >= 0: for u < 0 it drops below
     * 1, the Grueneisen term b/w0 grows without bound and divides by zero at u = -u0 eta^2. A
     * negative u is below the cold curve and has no Tillotson meaning, so w0 is evaluated with
     * max(u, 0), i.e. w0 >= 1. P stays continuous at u = 0 (u still multiplies Gamma), and the
     * derivative formulas below remain exact because they reduce to the w0 = const case at w0 = 1.
     *
     * For u <= 0 w0 is set to 1 directly: below rho ~ 1e-162 rho0 (in Earth-like units; where
     * u0*eta^2 underflows to 0) max(u, 0)/(u0 eta^2) was 0/0 and P NaN. For u > 0 the expression
     * is the same as before; there w0 overflows to inf below ~1e-155 rho0, which gives b/w0 = 0
     * (the correct limit) and is guarded in the derivatives below. A NaN u gives w0 = 1 as before
     * (fmax(NaN, 0) = 0), and NaN in P through u itself.
     */
    w0 = (u > 0.0) ? u/(u0*eta*eta) + 1.0 : 1.0;

    bCold = (rho >= rho0 || u < material->us);  /* Pc only */
    bHot = (!bCold && u > material->us2);        /* Pe only */

    /* Cold term Pc: needed in the cold and the intermediate region. */
    if (!bHot) {
        double Gammac = a + b/w0;
        int bCut = 0;

        Pc = Gammac*u*rho + A*mu + B*mu*mu;

        /*
         * The cut of the cold term is applied in the intermediate region in both modes: there Pc
         * only serves as the lower end of the blend and a negative value has no meaning. Only the
         * pure cold region is left raw for the inversions.
         */
        if (bClip || !bCold) {
#ifdef TILL_PRESS_MELOSH
            bCut = (eta < 0.8);
#else
            /*
             * Pc == 0 below rho0 is cut too. Below eta ~ 1e-16, mu = eta - 1 rounds to -1 and for
             * A == B (granite, basalt, dunite, forsterite, ice) A mu + B mu^2 cancels to exactly 0,
             * absorbing Gamma u rho: the true Pc ~ eta (Gamma u rho0 - A) is negative, but Pc == 0
             * escaped the cut, and dP/drho was (A - 2B)/rho0 + a u < 0 for P = 0. Not at rho0, where
             * Pc = 0 at u = 0 is the reference state with the slope A/rho0. P is the same either way.
             */
            bCut = (Pc < 0.0) || (Pc == 0.0 && rho < rho0);
#endif
        }
        if (bCut) {
            Pc = 0.0;       /* dPcdrho = dPcdu = 0 */
        } else if (bDeriv) {
            /* w0 = inf (u > 0 at rho below ~1e-155 rho0): b u (3 w0 - 2)/w0^2 -> 0, written out as
             * 0*inf = NaN otherwise. Unchanged for every finite w0. */
            dPcdrho = a*u + (isinf(w0) ? 0.0 : b*u/(w0*w0)*(3.0*w0 - 2.0)) + (A + 2.0*B*mu)/rho0;
            dPcdu = a*rho + b*rho/(w0*w0);
        }
    }

    /* Expanded term Pe: needed in the hot and the intermediate region. */
    if (!bCold) {
        /* exp(-beta z^2) is computed once (it was written out three times, for Pe and both derivatives). */
        const double eb = exp(-beta*z*z);
        const double Gammae = a + b/w0*eb;

        ex = exp(-(alpha*z + beta*z*z));
        Pe = Gammae*u*rho + A*mu*ex;
        if (bDeriv) {
            /*
             * With z = 1/eta - 1: dz/deta = -1/eta^2, dw0/deta = -2(w0-1)/eta, so
             * d/deta [mu exp(-(alpha z + beta z^2))] = exp(...) (1 + mu (alpha + 2 beta z)/eta^2).
             * (The previous expression had a spurious "- eta" inside the bracket.)
             */
            /*
             * The exponentials underflow to 0 already at eta ~ 0.1 (z ~ 10), while the factors they
             * multiply grow like 1/eta^2 and 1/eta^3 and overflow below eta ~ 4e-103: 0*inf gave a
             * NaN dP/drho there, and c fell to its floor sqrt(A/rho0) because fmax drops the NaN.
             * The terms are 0 wherever their exponential is 0; unchanged otherwise.
             */
            dPedrho = (a + ((eb > 0.0) ? b/w0*eb*(2.0*beta*z/eta + 3.0 - 2.0/w0) : 0.0))*u
                    + ((ex > 0.0) ? A/rho0*ex*(1.0 + mu/(eta*eta)*(alpha + 2.0*beta*z)) : 0.0);
            dPedu = (a + b/(w0*w0)*eb)*rho;
        }
    }

    /* Combine. */
    if (bCold) {
        P = Pc; dPdrho = dPcdrho; dPdu = dPcdu;
    } else if (bHot) {
        P = Pe; dPdrho = dPedrho; dPdu = dPedu;
    } else {
        /* y depends on u, which contributes (Pe - Pc) dy/du to dP/du. */
        y = (u - material->us)/(material->us2 - material->us);
        P = (1.0 - y)*Pc + y*Pe;
        if (bDeriv) {
            dPdrho = (1.0 - y)*dPcdrho + y*dPedrho;
            dPdu = (1.0 - y)*dPcdu + y*dPedu + (Pe - Pc)/(material->us2 - material->us);
        }
    }

    /* Final cutoff of the total pressure (can act in the hot and intermediate region). */
    if (bClip && P < 0.0) {
        P = 0.0;
        dPdrho = 0.0;
        dPdu = 0.0;
    }

    if (!bDeriv) return P;

    if (pcSound != NULL) {
        /* c^2 = (dP/drho)_s = (dP/drho)_u + P/rho^2 (dP/du)_rho, from the derivatives above. */
        /*
         * rho^2 underflows below rho ~ 1e-154 (code units): then P/rho^2 is inf or NaN, so divide
         * twice. At denormal rho P/rho can still overflow while dP/du = a rho underflows to 0; the
         * term is 0 then (it was inf*0). Unchanged wherever rho^2 is a normal double.
         */
        const double rho2 = rho*rho;
        c2 = dPdrho + ((rho2 >= DBL_MIN) ? P/rho2*dPdu : ((P == 0.0 || dPdu == 0.0) ? 0.0 : (P/rho)*(dPdu/rho)));
        /*
         * Floor at the bulk sound speed of the uncompressed material (always with the cutoff). A NaN
         * (NaN u) is passed on like in P: fmax would replace it by the floor, a finite c for a NaN state.
         */
        if (c2 == c2) c2 = bClip ? fmax(c2, A/rho0) : fmax(c2, 0.0);
        *pcSound = sqrt(c2);
    }
    if (pdPdrho != NULL) *pdPdrho = dPdrho;
    if (pdPdu != NULL) *pdPdu = dPdu;

    return P;
}

/* Pressure and sound speed with the pressure cutoff. Set pcSound = NULL for the pressure only. */
double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
    return tillPressureCore(material, rho, u, 1, pcSound, NULL, NULL);
}

/* Pressure and sound speed of the raw EOS (no cutoff, no sound speed floor). */
double tillPressureSoundNP(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
    return tillPressureCore(material, rho, u, 0, pcSound, NULL, NULL);
}

double tillPressure(TILLMATERIAL *material, double rho, double u)
{
    return tillPressureCore(material, rho, u, 1, NULL, NULL, NULL);
}

double tillPressureNP(TILLMATERIAL *material, double rho, double u)
{
    return tillPressureCore(material, rho, u, 0, NULL, NULL, NULL);
}

double tillSoundSpeed(TILLMATERIAL *material, double rho, double u)
{
    double c;

    tillPressureCore(material, rho, u, 1, &c, NULL, NULL);
    return c;
}

/* dP/drho at u = const, of the cut pressure. */
double tilldPdrho(TILLMATERIAL *material, double rho, double u)
{
    double dPdrho;

    tillPressureCore(material, rho, u, 1, NULL, &dPdrho, NULL);
    return dPdrho;
}

/* dP/du at rho = const, of the cut pressure. */
double tilldPdu(TILLMATERIAL *material, double rho, double u)
{
    double dPdu;

    tillPressureCore(material, rho, u, 1, NULL, NULL, &dPdu);
    return dPdu;
}

/*
 * Calculate T(rho,u) for a material. As an approximation
 * we use u(rho,T) = uc(rho) + cv*T.
 *
 * Never fails and never negative (EOSlib contract for T(rho, u), same as ANEOSTofRhoU): where
 * u(rho, T) cannot be inverted it returns T = 0, the cold curve, which is the lowest temperature
 * of this EOS. That is the case outside the table in rho (counted by tillColdULookup) and for a
 * state below the cold curve, u < u_cold(rho), or a non-finite u (counted here). A state within
 * the roundoff tolerance tillUTol below the cold curve is on it: T = 0, not counted.
 */
double tillTempRhoU(TILLMATERIAL *material, double rho, double u)
{
    double ucold = tillColdULookup(material, rho);

    if (tillFailed(ucold)) return 0.0;

    if (!(u >= ucold - tillUTol(material, u)) || !isfinite(u))
        TILL_FAIL_RETURN(material, 0.0, "tillTempRhoU: (rho= %g, u= %g, iMat= %i) below the cold curve u_cold= %g (or u not finite), T = 0.\n",
                         rho, u, material->iMaterial, ucold);

    return fmax((u-ucold)/material->cv, 0.0);
}

/*
 * Calculate u(rho,T) for a material.
 *
 * Returns TILL_FAIL_VALUE (counted) for T below the lowest temperature of this EOS, T < 0, and
 * for a non-finite T, like u(rho, T) of the other libraries below T_min (igeosUofRhoT, ANEOS,
 * SCvH, REOS3); it returned u_cold + cv T, a state below the cold curve. Also TILL_FAIL_VALUE
 * if the cold curve fails at rho (counted there).
 */
double tillURhoTemp(TILLMATERIAL *material, double rho, double T)
{
    double ucold;

    if (!(T >= 0.0) || !isfinite(T))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillURhoTemp: T= %g is negative or not finite (rho= %g, iMat= %i).\n",
                         T, rho, material->iMaterial);

    ucold = tillColdULookup(material, rho);
    if (tillFailed(ucold)) return TILL_FAIL_VALUE;
    return(ucold + material->cv*T);
}

/*
 * Brent's method (inverse quadratic interpolation with bisection fallback), own implementation.
 */
int tillBrent(double (*f)(TILLMATERIAL *, double, const void *), TILLMATERIAL *material,
              const void *ctx, double a, double b, double fa, double fb,
              double tol_abs, double tol_rel, int max_iter, double *px)
{
    double c, fc, d, e, m, tol, p, q, r, s;
    int iter;

    if (!(fa*fb <= 0.0)) return TILL_FAIL;      /* not bracketed, or NaN */
    if (fa == 0.0) { *px = a; return TILL_SUCCESS; }
    if (fb == 0.0) { *px = b; return TILL_SUCCESS; }

    c = a; fc = fa;
    d = e = b - a;

    for (iter = 0; iter < max_iter; iter++) {
        /* Keep b the best estimate, c on the other side of the root. */
        if (fb*fc > 0.0) { c = a; fc = fa; d = e = b - a; }
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        tol = tol_abs + tol_rel*fabs(b);
        m = 0.5*(c - b);
        if (fabs(m) <= tol || fb == 0.0) { *px = b; return TILL_SUCCESS; }

        if (fabs(e) >= tol && fabs(fa) > fabs(fb)) {
            /* Try interpolation: secant if a == c, inverse quadratic otherwise. */
            s = fb/fa;
            if (a == c) {
                p = 2.0*m*s;
                q = 1.0 - s;
            } else {
                q = fa/fc;
                r = fb/fc;
                p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
                q = (q - 1.0)*(r - 1.0)*(s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;
            if (2.0*p < fmin(3.0*m*q - fabs(tol*q), fabs(e*q))) {
                e = d; d = p/q;         /* accept interpolation */
            } else {
                d = m; e = m;           /* bisect */
            }
        } else {
            d = m; e = m;
        }
        a = b; fa = fb;
        b += (fabs(d) > tol) ? d : (m > 0.0 ? tol : -tol);
        fb = f(material, b, ctx);
    }
    return TILL_FAIL;
}

/* f(rho) = P(rho, T) - P for the root finder in tillRhoPTemp (no pressure cutoff). */
struct tillRhoPTempCtx { double P; double T; };

static double tillRhoPTempFunc(TILLMATERIAL *material, double rho, const void *ctx)
{
    const struct tillRhoPTempCtx *p = (const struct tillRhoPTempCtx *) ctx;
    double u = tillURhoTemp(material, rho, p->T);
    if (tillFailed(u)) return NAN;     /* propagates as a bracket failure */
    return tillPressureSoundNP(material, rho, u, NULL) - p->P;
}

/*
 * Calculate rho(P,T) for a material by root finding on the uncut pressure (tillPressureSoundNP),
 * which is monotonic in rho at fixed T where (dP/drho)_T > 0.
 *
 * Returns TILL_FAIL_VALUE if P <= 0, T < 0 or not finite, or the root is not bracketed by
 * [rho_min, 0.999*rhomax], with rho_min = TILL_RHOPT_RHOMIN_FAC*rho0 (the same physical bracket in
 * every unit system). Counted once per call.
 */
double tillRhoPTemp(TILLMATERIAL *material, double P, double T)
{
    struct tillRhoPTempCtx ctx;
    const double rho_min = TILL_RHOPT_RHOMIN_FAC*material->rho0;   /* P diverges for rho = 0 and u = 0 */
    const double rho_max = 0.999*material->rhomax;
    double fa, fb, rho;

    /* For P <= 0 the inversion is not defined (pressure cutoff). */
    if (P <= 0.0)
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillRhoPTemp: Called for negative pressure (P= %g).\n", P);

    /* T < 0 (or not finite): u(rho, T) fails at every density. Checked here, so that it is counted
     * once and not by each of the bracket evaluations and the bracket failure. */
    if (!(T >= 0.0) || !isfinite(T))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillRhoPTemp: T= %g is negative or not finite (P= %g).\n", T, P);

    ctx.P = P;
    ctx.T = T;

    fa = tillRhoPTempFunc(material, rho_min, &ctx);
    fb = tillRhoPTempFunc(material, rho_max, &ctx);

    if (!(fa < 0.0))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillRhoPTemp: P <= P(rho_min, T) (P= %g, rho_min= %g, P_min= %g, T= %g).\n",
                         P, rho_min, fa + P, T);
    if (!(fb > 0.0))
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillRhoPTemp: P >= P(rho_max, T) (P= %g, rho_max= %g, P_max= %g, T= %g).\n",
                         P, rho_max, fb + P, T);

    if (tillBrent(tillRhoPTempFunc, material, &ctx, rho_min, rho_max, fa, fb, 0.0, 1e-10, 1000, &rho) != TILL_SUCCESS)
        TILL_FAIL_RETURN(material, TILL_FAIL_VALUE, "tillRhoPTemp: root finder did not converge (P= %g, T= %g).\n", P, T);

    return rho;
}

/*
 * Calculate the derivative of u with respect to rho.
 *
 * du/drho = P/rho^2
 *
 */
double tilldudrho(TILLMATERIAL *material, double rho, double u)
{
    const double P = tillPressure(material,rho,u);
    /*
     * The sentinel of an invalid rho (counted in tillPressureCore) is passed on, never divided by rho^2.
     * Only the sentinel: a NaN P (from a NaN u) stays NaN like P, c and dP/drho (tillFailed is true for
     * NaN as well, which made du/drho the sentinel there, uncounted).
     */
    if (P == TILL_FAIL_VALUE) return TILL_FAIL_VALUE;
    return (P/(rho*rho));
}

/*
 * Calculate the logarithmic derivative of u with respect to rho.
 *
 * du/dlogrho = P/rho
 *
 */
double tilldudlogrho(TILLMATERIAL *material, double logrho, double u)
{
    const double rho = exp(logrho);
    const double P = tillPressure(material,rho,u);
    if (P == TILL_FAIL_VALUE) return TILL_FAIL_VALUE;   /* the sentinel only, as in tilldudrho; rho = exp(logrho) is 0 or inf only for |logrho| > ~700 */
    return (P/(rho));
}

/*
 * Determine the region of the EOS.
 */
int tillPhase(TILLMATERIAL *material, double rho, double u)
{
    /*
     *  Evaluate which part of the equation of state we need.
     */
    if (rho >= material->rho0) {
        /*
         *  Condensed states (rho > rho0).
         */
        return TILL_PHASE_CONDENSED;
    } else if (u < material->us) {
        /*
         * Expanded cold states (rho < rho0 and u < us).
         */
        return TILL_PHASE_EXP_COLD;
    } else if (u > material->us2) {
        /*
         * Expanded hot states (rho < rho0 and u > us2).
         */
        return TILL_PHASE_EXP_HOT;
    } else {
        /*
         *  intermediate states (rho < rho0 and us < u < us2)
         */
        return TILL_PHASE_EXP_INT;
    }

    return -1;
}