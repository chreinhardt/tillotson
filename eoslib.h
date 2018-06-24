/*
 * The header file for eoslib.c.
 */
#ifndef EOSLIB_HINCLUDED
#define EOSLIB_HINCLUDED

/*
 * Define the range of material ID for each EOS.
 */
static int IDEALGAS      = 0
static int EOS_TILL_MIN  = 1
static int EOS_TILL_MAX  = 30
static int EOS_ANEOS_MIN = 31
static int EOS_ANEOS_MAX = 50
static int VDW_EOS       = 51

static int EOS_N_MAT_MAX = 51
if 0
#define IDEALGAS      0
#define EOS_TILL_MIN  1
#define EOS_TILL_MAX  20
#define EOS_ANEOS_MIN 21
#define EOS_ANEOS_MAX 50
#endif

void *eosInitMat(int iMat, void *param);
void eosFinalizeMat(void *eosMat);

double eosPressureSound(void *eosMat, double rho, double u, double *pcSound);
double eosPressure(void *eosMat, double rho, double u);

double eosdPdrho(void *eosMat, double rho, double u);
double eosdPdu(void *eosMat, double rho, double u);

double eosTempRhoU(void *eosMat, double rho, double u);
double eosRhoPU(void *eosMat, double P, double u);
#endif

