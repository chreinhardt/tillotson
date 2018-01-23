/*
 * Copyright (c) 2018 Christian Reinhardt.
 *
 * This file provides a general interface that allows the use of different
 * equations of state (EOS) in hydro-codes.
 *
 * Currently supported are:
 *
 * Ideal gas EOS
 * Tillotson EOS (Tillotson 1962)
 * Van der Waals EOS
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "eoslib.h"

/*
 * Basic functions:
 *
 * eosInitMat: initialise a material of a given EOS from its material ID.
 *
 * eosFinalizeMat: finalize a material, free memory if needed.
 *
 * eosPressureSound: calculate the pressure and soundspeed for a given EOS/material.
 *
 * eosPressure: calculate the pressure only (uses eosPressureSound).
 *
 * eosSoundSpeed: calculate the sound speed for a given rho and u (uses eosPressureSound).
 *
 * eosFinalize: free memory.
 */

void *eosInitMat(int iMat, double dKpcUnit, double dMsolUnit, void *param)
{
	/*
	 * Initialize an EOS/material.
	 *
	 * Parameters:
     *
     * iMat:        material ID
     * params:      a general array containing different parameters that depend
     *              on the EOS.
	 */

    if (iMat == IDEALGAS)
    {
        // Ideal gas
        return((igMat *) igInit());
    } else if (EOS_TILL_BEGIN <= iMat <= EOS_TILL_END) {
        // Tillotson EOS
        return((TILLMATERIAL*) (tillInitMaterial(iMat, dKpcUnit, dMsolUnit, param->nTableRho, param->nTableV, param->rhomax, param->vmax, param->iExpV)));
    } else {
        fprintf(stderr, "iMat= %i undefined.\n", iMat);
        assert(0);
    }
}

void eosFinalizeMat(void *eosMat)
{
	/*
     * Free the memory.
     */

	free(eosMat);
}

double eosPressureSound(void *eosMat, double rho, double u, double *pcSound)
{
    /*
     * Calculate the pressure and sound speed for a given EOS and material.
     *
     * Input:
     * void *eosMat
     *
     * Output:
     * Pressure
     */
}

double eosPressure(void *eosMat, double rho, double u)
{
    /*
     * Calculate the pressure for a given EOS and material.
     */
	return (eosPressureSound(eosMat, rho, u, NULL));
}

double eosdPdrho(void *eosMat, double rho, double u)
{
    /*
     * Calculate dP/drho (at constant u).
     */
}

double eosdPdu(void *eosMat, double rho, double u)
{
    /*
     * Calculate dP/du (at constant rho).
     */
}

double eosTempRhoU(void *eosMat, double rho, double u)
{
    /*
     * Calculate T(rho, u).
     */
}

double eosRhoPU(void *eosMat, double P, double u)
{
    /*
     * Calculate rho(P, u).
     */
}

