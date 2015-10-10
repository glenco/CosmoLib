#include <iostream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <cstddef>
#ifdef ENABLE_GSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#endif

#ifndef pi
#define pi  3.141593
#endif

#ifndef Grav
#define Grav 4.7788e-20  // Newton's over c^2 in Mpc / Msun
#endif

#ifndef lightspeed
#define lightspeed 2.99792458e5
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::cout << "ERROR: file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#define CRITD 2.49783e18    /* critical density / Ho^2 solar/Mpc */
#define CRITD2 2.7752543e11 /* critical density / h^2 M_sun/Mpc^3 */

#ifndef cosmo_declare

enum CosmoParamSet {WMAP5yr,Millennium,Planck1yr};

/** \ingroup cosmolib
 *
 * \brief The cosmology and all the functions required to calculated quantities based on the cosmology.
 *
 * This class is used to store the cosmological parameters, calculate cosmological distances, calculate
 * the power spectrum of density fluctuations and the mass function of halos.
 *
 * As set now, there are no baryon acoustic oscillations in the power spectrum, but this can be changed
 * at the expense of not including neutrinos.
 */
class COSMOLOGY{
public:

	COSMOLOGY(CosmoParamSet cosmo_p = WMAP5yr);
	COSMOLOGY(double omegam,double omegal,double h, double w);
	~COSMOLOGY();

    void SetConcordenceCosmology(CosmoParamSet cosmo_p = WMAP5yr);
    void PrintCosmology(short physical = 0) const;

    // Lengths
    double rcurve() const;
    //double emptyDist(double zo,double z);
    double coorDist(double zo,double z) const;
    double coorDist(double z) const {return coorDist(0,z);}
    double radDist(double zo,double z) const;
    double radDist(double z) const {return radDist(0,z);}
    double angDist(double zo,double z) const;
    double angDist(double z) const {return angDist(0,z);}
    double lumDist(double zo,double z) const;
    double lumDist(double z) const {return lumDist(0,z);}
  
    double DRradius(double zo,double z,double pfrac);
    double DRradius2(double zo,double z);

    double invCoorDist(double d) const;
    double invRadDist(double d) const;

    double scalefactor(double rad) const;
    double Omegam(double z) const;
    double rho_crit(double z) const;
    //double drdz_empty(double x);
    double drdz(double x) const;
    double adrdz(double x) const;
    double drdz_dark(double x) const;
    double adrdz_dark(double x) const;
    double DeltaVir(double z,int caseunit=0) const;
    double Deltao(double m) const;
    double time(double z) const;
    double nonlinMass(double z) const;
  
    // Stuff having to do with the power spectrum
    double power_normalize(double sigma8);
    double power_linear(double k,double z);
    double Dgrowth(double z) const;
    double powerCDMz(double k,double z);
    double psdfdm(double z,double m,int caseunit=0);
    double halo_bias (double m, double z, int t=0);
    double stdfdm(double z,double m,int caseunit=0);
    double powerlawdfdm(double z,double m,double alpha,int caseunit=0);
    double haloNumberDensity(double m,double z,double a, int t,double alpha = 0.0);
	
    double haloNumberDensityOnSky (double m,double z1,double z2,int t,double alpha = 0.0);
    double haloNumberInBufferedCone (double mass ,double z1,double z2,double fov,double buffer ,int type ,double alpha=0.0);
    double haloMassInBufferedCone (double mass ,double z1,double z2,double fov,double buffer ,int type ,double alpha=0.0);
	
    double TopHatVariance(double m) const;
    double TopHatVarianceR(double R,double z);
    double TopHatVarianceM(double M,double z);
  //double gradius(double R,double rd);

    double getTimefromDeltaC(double dc);
    double getZfromDeltaC(double dc);
    
    /// accesser functions:

    /// Hubble paremters in units of 100 km/s/Mpc, renormalizes P(k) to keep sig8 fixed
    void sethubble(double ht){ h = ht; TFmdm_set_cosm(); power_normalize(sig8);}
    double gethubble() const {return h;}
    /// Hubble parameter in 1/Mpc units
    double getHubble() const {return 100*h/lightspeed;}
  
    /// Primordial spectral index, renormalizes P(k) to keep sig8 fixed
    void setindex(double nn){ n = nn;  power_normalize(sig8);}
    double getindex() const { return n;}
    
    /// Omega matter, renormalizes P(k) to keep sig8 fixed
    void setOmega_matter(double Omega_matter,bool FLAT = false){Omo = Omega_matter; if(FLAT) Oml = 1-Omo ; TFmdm_set_cosm(); power_normalize(sig8);}
    double getOmega_matter() const {return Omo;}
    
    /// Omega lambda, renormalizes P(k) to keep sig8 fixed
    void setOmega_lambda(double Omega_lambda,bool FLAT = false){Oml = Omega_lambda;  if(FLAT) Oml = 1-Omo ; TFmdm_set_cosm(); power_normalize(sig8);}
    double getOmega_lambda() const {return Oml;}

    /// Omega baryon, renormalizes P(k) to keep sig8 fixed
    void setOmega_baryon(double Omega_baryon){Omb = Omega_baryon; TFmdm_set_cosm(); power_normalize(sig8);}
    double getOmega_baryon() const {return Omb;}
    
    /// Omega neutrino, renormalizes P(k) to keep sig8 fixed
    void setOmega_neutrino(double Omega_neutrino){Omnu = Omega_neutrino; TFmdm_set_cosm(); power_normalize(sig8);}
    double getOmega_neutrino() const {return Omnu;}

    /// Number of neutrino species, renormalizes P(k) to keep sig8 fixed
    void setNneutrino(double Nneutrino){Nnu = Nneutrino; TFmdm_set_cosm(); power_normalize(sig8);}
    double getNneutrino() const {return Nnu;}
	
   /// Dark energy equation of state parameter p/rho = w + w_1 (1+z)
    void setW(double ww){w = ww;}
    double getW() const {return w;}
    void setW1(double ww){w1 = ww;}
    double getW1() const {return w1;}
    
    /// Running of primordial spectral index, P(k)_primordial \propto pow(k/h,n+dndlnk*log(k)), renormalizes P(k) to keep sig8 fixed
    void setdndlnk(double w){dndlnk = w; power_normalize(sig8);}
    double getdndlnk() const {return dndlnk;}

    /// Alternative to w for dark energy/ alt. grav. structure evolution
    void setgamma(double gamm){gamma = gamm;}
    double getgamma() const {return gamma;}
    // If physical = 1 all Omega are Omega*h^2 otherwise they have the usual definitions.
    // 2 gamma parameterization is used for dark energy
    void setDEtype(short tt){darkenergy = tt;}
    short getDEtype() const {return darkenergy;}
  
    void setSigma8(double my_sig8){power_normalize(my_sig8);}
    double getSigma8() const {return sig8;}


    void dzdangDist(double D,double z[],double dzdD[]);

    double totalMassDensityinHalos(int t,double alpha,double m_min,double z,double z1,double z2);

    /// set interpolation range
    void setInterpolation(double z_interp);
    /// set interpolation range and number of points
    void setInterpolation(double z_interp, std::size_t n_interp);

  /// The lensing critical density in Msun / Mpc^2
   double SigmaCrit(double zlens,double zsource) const;
    
protected:

	/// Hubble paremters in units of 100 km/s/Mpc
  double h;
    /// Primordial spectral index
  double n;
    /// Omega matter, dark matter + baryons
  double Omo;
    /// Omega lambda
  double Oml;
    /// Omega baryon
  double Omb;
    /// Omega neutrino
  double Omnu;
   /// Number of neutrino species
  double Nnu;
   /// Dark energy equation of state parameter p/rho = w + w_1 (1+z)
  double w;
  /// Dark energy equation of state parameter p/rho = w + w_1 (1+z)
  double w1;
   /// Running of primordial spectral index, P(k)_primordial \propto pow(k/h,n+dndlnk*log(k))
  double dndlnk;

  /// Alternative to w for dark energy/ alt. grav. structure evolution
  double gamma;
  // If physical = 1 all Omega are Omega*h^2 otherwise they have the usual definitions.
  //short physical;
  // 2 gamma parameterization is used for dark energy
  short darkenergy;

  /* table for growth parameter */
  //std::auto_ptr<double> *aa;
  //std::auto_ptr<double> *growth;
  //int Ntable;
  
  double A;
  double sig8;  /* do not access these normalization outside */
  
  double Rtophat;
  double ztmp;
  
  double powerEH(double k,double z);
  double powerEHv2(double k);
  double powerloc(double k,double z);
  double npow(double k);
  double normL(double lgk);
  double gradius(double R,double rd) const;

  int nbin;  // number of points when binning to interpolate
  std:: vector<double> vDeltaCz,vlz,vt;
  double DpropDz(double z);
  double dsigdM(double m);
  double timeEarly(double a) const;
  double dNdz(double z);
  
  typedef double (COSMOLOGY::*pt2MemFunc)(double) const;
  typedef double (COSMOLOGY::*pt2MemFuncNonConst)(double);
  
  double nintegrateDcos(pt2MemFunc func,double a,double b,double tols) const;
  double trapzdDcoslocal(pt2MemFunc func, double a, double b, int n, double *s2) const;
  double nintegrateDcos(pt2MemFuncNonConst func,double a,double b,double tols);
  double trapzdDcoslocal(pt2MemFuncNonConst func, double a, double b, int n, double *s2);
  double dfridrDcos(pt2MemFunc func, double x, double h, double *err);
  double f4(double u) const;

  int ni;
  std::vector<float> xf, wf;
  
  // The parameters used in Eisenstein & Hu power spectrum
    /* in powerEH.c */
  short TFmdm_set_cosm_change_z(double redshift);
  short TFmdm_set_cosm();
  double TFmdm_onek_mpc(double kk);
  double TFmdm_onek_hmpc(double kk);
  
  double f_baryon;	// Baryon fraction 
  double f_bnu;		// Baryon + Massive Neutrino fraction
  double f_cb;		// Baryon + CDM fraction
  double f_cdm;		// CDM fraction
  double f_hdm;		// Massive Neutrino fraction
  
  double alpha_gamma;	 // sqrt(alpha_nu)
  double alpha_nu;	     // The small-scale suppression
  double beta_c;		 // The correction to the log in the small-scale
  double num_degen_hdm;	 // Number of degenerate massive neutrino species
  double growth_k0;	     // D_1(z) -- the growth function as k->0
  double growth_to_z0;	 // D_1(z)/D_1(0) -- the growth relative to z=0
  double k_equality;	 // The comoving wave number of the horizon at equality
  double obhh;		     // Omega_baryon * hubble^2 
  double omega_curv;	 // = 1 - omega_matter - omega_lambda
  double omega_lambda_z; // Omega_lambda at the given redshift
  double omega_matter_z; // Omega_matter at the given redshift
  double omhh;		     // Omega_matter * hubble^2
  double onhh;		     // Omega_hdm * hubble^2
  double p_c;		     // The correction to the exponent before drag epoch
  double p_cb;		     // The correction to the exponent after drag epoch 
  double sound_horizon_fit;  // The sound horizon at the drag epoch 
  double theta_cmb;	     // The temperature of the CMB, in units of 2.7 K
  double y_drag;		 // Ratio of z_equality to z_drag 
  double z_drag;		 // Redshift of the drag epoch
  double z_equality;	 // Redshift of matter-radiation equality
    
    // in powerEHv2.c
  void TFset_parameters(double omega0hh, double f_baryon, double Tcmb);
  double TFfit_onek(double k, double *tf_baryon, double *tf_cdm);

	double R_drag;		// Photon-baryon ratio at drag epoch
	double R_equality;	// Photon-baryon ratio at equality epoch
	double sound_horizon;	// Sound horizon at drag epoch, in Mpc
	double k_silk;		// Silk damping scale, in Mpc^-1
	double alpha_c;	    // CDM suppression
	double alpha_b;	    // Baryon suppression
	double beta_b;		// Baryon envelope shift
	double beta_node;	// Sound horizon shift
	double k_peak;		// Fit to wavenumber of first peak, in Mpc^-1

  // temporary variables for doing interations
  int tmp_type;
  double tmp_alpha;
  double tmp_mass;
  double tmp_a;

  // interpolation of functions
  double z_interp;
  std::size_t n_interp;
  void calc_interp(double z_interp, std::size_t n_interp);
  double interp(const std::vector<double>& table, double z) const;
  double invert(const std::vector<double>& table, double f_z) const;
  std::vector<double> redshift_interp;
  std::vector<double> coorDist_interp;
  std::vector<double> radDist_interp;
};

typedef COSMOLOGY *CosmoHndl;
/**
 *  \brief Class for calculating properties of NFW halo profile.
 *   
 *   This class does not take into affect the cosmological correlation between concentration and mass.
 *   For this see the HALOCalculator class.
 */
class NFW_Utility {
public:
	NFW_Utility(){return;}
	~NFW_Utility(){};

	// methods for NFW profile
	double NFW_V200(double M200,double R200);
	double NFW_Vmax(double cons,double M200,double R200);
	double NFW_Vr(double x,double cons,double M200,double R200);
	double NFW_deltac(double cons);
	double NFW_Concentration(double Vmax,double M200,double R200);
	double NFW_rho(double cons,double x);

	void match_nfw(float Vmax,float R_half,float mass,float *cons,float *Rsize);
	float Rsize(float cons,float Vmax,float mass);
	float g_func(float x);

private:
	float Vmax,R_half,mass;       /// Mass (solar masses)

    typedef float (NFW_Utility::*MemFunc)(float);
    float zbrentD(MemFunc func, float a,float b,float tols);
    float nfwfunc(float cons);
    float funcforconcentration(float cons);
};


/// wrapper functions for gsl integration / differentiation
double drdz_wrapper(double x, void *params);
double drdz_dark_wrapper(double x, void *params);
double Deltao_wrapper(double m, void *params);

#define cosmo_declare
#endif

/*** in cosmo.c ***/
int cosmo_compare(CosmoHndl cos1,CosmoHndl cos2);
void cosmo_copy(CosmoHndl cos1,CosmoHndl cos2);
void ders(double z,double Da[],double dDdz[]);
void dir(double r,double a[],double dadr[]);
double arctanh(double x);
double fmini(double a,double b);
double fmaxi(double a,double b);
double **dmatrixcos(long nrl, long nrh, long ncl, long nch);
void free_dmatrixcos(double **m, long nrl, long nrh, long ncl, long nch);

/* in nfw.c */

/* in powerEHv2.c */
void TFset_parameters(double omega0hh, double f_baryon, double Tcmb);
double TFfit_onek(double k, double *tf_baryon, double *tf_cdm);
