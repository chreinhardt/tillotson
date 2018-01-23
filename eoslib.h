/*
 * The header file for eoslib.c.
 */
#ifndef EOSLIB_HINCLUDED
#define EOSLIB_HINCLUDED

/*
 * Define the range of material ID for each EOS.
 */
#define IDEALGAS    0
#define TILL_ID_MIN 1
#define TILL_ID_MAX 100

void *eosInitMat(int iMat, void *param);
void eosFinalizeMat(void *eosMat);

double eosPressureSound(void *eosMat, double rho, double u, double *pcSound);
double eosPressure(void *eosMat, double rho, double u);

double eosdPdrho(void *eosMat, double rho, double u);
double eosdPdu(void *eosMat, double rho, double u);

double eosTempRhoU(void *eosMat, double rho, double u);
double eosRhoPU(void *eosMat, double P, double u);
#endif

