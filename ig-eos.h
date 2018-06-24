/*
 * The header file for the ideal gas EOS.
 */
#ifndef IG_EOS_HINCLUDED
#define IG_EOS_HINCLUDED

/*
 * Material id
 */
#define IDEALGAS     0
#define IG_EOS_N_MAT 1

typedef struct igMat
{
	int iMaterial;	/* What material is it? */
	double rhomax;	/* Max value (if b > 0) */
	
    /* Unit convertion factors */
	double dKpcUnit;
	double dMsolUnit;
	
	double dGasConst;
	double dErgPerGmUnit;
	double dGmPerCcUnit;
	double dSecUnit;
    
    /*
     * The constants needed for the ideal gas EOS.
     */
    double dConstGamma;
    double dMeanMolMass;        /* in units of m_hydrogen */
    double b;                   /* particle volume */
    double rho0;

	 /* the specific heat capacity (assumed constant) */
	double cv;
} IGMAT;

// Initialize a material
IGMAT *igInitMaterial(double dConstGamma, double dMeanMolMass, double b, double dKpcUnit, double dMsolUnit);
void igFinalizeMaterial(TILLMATERIAL *material);

double igPressureSound(IGMAT *material, double rho, double u, double *pcSound);
double igPressure(IGMAT *material, double rho, double u);
double igdPdrho(IGMAT *material, double rho, double u);
double igdPdu(IGMAT *material, double rho, double u);
double igTempRhoU(IGMAT *material, double rho, double u);
double igRhoPU(IGMAT *material, double P, double u);
double igURhoP(IGMAT *material, double rho, double P);
double igPhi(IGMAT *material, double rho, double u);
double igGamma(IGMAT *material, double rho, double u);
double igLookupU(IGMAT *material, double rho1, double u1, double rho2, int iOrder);
#endif

