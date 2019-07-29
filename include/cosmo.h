#include <iostream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <cstddef>
#ifdef ENABLE_GSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <cmath>

#endif

#ifndef PI
#define PI  3.141593
#endif

#ifndef Grav
#define Grav 4.7788e-20  // Newton's over c^2 in Mpc / Msun
#endif

#ifndef lightspeed
#define lightspeed 2.99792458e5  // km / s
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::cout << "ERROR: file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#define CRITD 2.49783e18    /* critical density / Ho^2 solar/Mpc */
#define CRITD2 2.7752543e11 /* critical density / h^2 M_sun/Mpc^3 */

#ifndef cosmo_declare

enum CosmoParamSet {WMAP5yr,Millennium,Planck1yr,Planck15,Planck18,BigMultiDark,none};

std::ostream &operator<<(std::ostream &os, CosmoParamSet p);

/**
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
  /// if justdistances== true
  COSMOLOGY(double omegam,double omegal,double h, double w = -1,double wa = 0
            ,bool justdistances=false /// if true the internals needed calculate structure formation will not be calculated making constructuion faster
  );
  COSMOLOGY(const COSMOLOGY &cosmo);
  ~COSMOLOGY();
  
  COSMOLOGY& operator=(const COSMOLOGY &cosmo);
  
  // returns the parameter set if any
  CosmoParamSet ParamSet(){return cosmo_set;}
  
  void PrintCosmology(short physical = 0) const;
  std::string print(short physical= 0) const;
  
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
  
  /// derivative of coordinate distance with respect to Omo in flat cosmology.  Useful for Fisher matrix calculations
  double d_coorDist_dOmo(double zo,double z) const;
  /// derivative of coordinate distance with respect to w.  Useful for Fisher matrix calculations
  double d_coorDist_dw(double zo,double z) const;
  /// derivative of coordinate distance with respect to w1.  Useful for Fisher matrix calculations
  double d_coorDist_dw1(double zo,double z) const;
  
  double DRradius(double zo,double z,double pfrac);
  double DRradius2(double zo,double z);
  
  double invCoorDist(double d) const;
  double invRadDist(double d) const;
  double invComovingDist(double d) const;
  
  double scalefactor(double rad) const;
  double Omegam(double z) const;
  double rho_crit_comoving(double z) const;
  //double drdz_empty(double x);
  double drdz(double x) const;
  //double adrdz(double x) const;
  double drdz_dark(double x) const;
  inline double dark_factor(double x) const;
  //double adrdz_dark(double x) const;
  
  /// the critical over density
  double delta_c() const;
  
  double DeltaVir(double z,int caseunit=0) const;
  double Deltao(double m) const;
  double time(double z) const;
  double nonlinMass(double z) const;

  // Stuff having to do with the power spectrum
  double power_normalize(double sigma8);
  /**
   * \brief Linear power spectrum P(k,z)/a^2
   */
  double power_linear(double k,double z);
  double Dgrowth(double z) const;
  
  /** 
   * \brief  powerCDM.c calculates the nonlinear P(k,z)/a(r)^2
   *
   * The method of Peacock & Dodds 1996 is used to convert the linear
   * power spectrum into the nonlinear one.
   * This could be updated to a more recent nonlinear power spectrum
   */
  double powerCDMz(double k,double z);
  double psdfdm(double z,double m,int caseunit=0);
  double halo_bias (double m, double z, int t=0);
  double stdfdm(double z,double m,int caseunit=0);
  double powerlawdfdm(double z,double m,double alpha,int caseunit=0);
  double haloNumberDensity(double m,double z,double a, int t,double alpha = 0.0);
  
  /** \brief Dark matter correlation function 
   *
   *  This integrates powerCDMz() to get the correlation function as a function of comoving radius.
   *  Care should be taken that the right range of k is integrated over if the radius is very small or large.
   */
  double CorrelationFunction(double radius,double redshift
                      ,double k_max = 100,double k_min = 1.0e-3);

  struct CorrFunctorType{
    CorrFunctorType(COSMOLOGY *cosmo,double radius,double redshift)
    : cosmology(cosmo),z(redshift),r(radius)
    {
      norm = 0.5/PI/PI/(1+z)/(1+z);
    };
    
    COSMOLOGY *cosmology;
    double z;
    double r;
    
    double norm;
    
    double operator () (double k) {
      
      double rk = r*k;
      double jo = sin(rk)/rk;
      if(rk < 1.0e-3) jo = 1.0;
      
      return norm*jo*k*k*cosmology->powerCDMz(k,z);
    }
  };

  
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
  void sethubble(double ht){ h = ht; TFmdm_set_cosm();
    power_normalize(sig8); 	calc_interp_dist(); cosmo_set = none;}
  double gethubble() const {return h;}
  /// Hubble parameter in 1/Mpc units
  double getHubble() const {return 100*h/lightspeed;}
  
  /// Primordial spectral index, renormalizes P(k) to keep sig8 fixed
  void setindex(double nn){ n = nn;  power_normalize(sig8); cosmo_set = none;}
  double getindex() const { return n;}
  
  /// Omega matter, renormalizes P(k) to keep sig8 fixed
  void setOmega_matter(double Omega_matter,bool FLAT = false){Omo = Omega_matter; if(FLAT) Oml = 1-Omo ; TFmdm_set_cosm(); power_normalize(sig8); calc_interp_dist(); cosmo_set = none;}
  double getOmega_matter() const {return Omo;}
  
  /// Omega lambda, renormalizes P(k) to keep sig8 fixed
  void setOmega_lambda(double Omega_lambda,bool FLAT = false){Oml = Omega_lambda;  if(FLAT) Oml = 1-Omo ; TFmdm_set_cosm(); power_normalize(sig8); calc_interp_dist(); cosmo_set = none;}
  double getOmega_lambda() const {return Oml;}
  
  /// Omega baryon, renormalizes P(k) to keep sig8 fixed
  void setOmega_baryon(double Omega_baryon){Omb = Omega_baryon; TFmdm_set_cosm(); power_normalize(sig8); calc_interp_dist(); cosmo_set = none;
  }
  double getOmega_baryon() const {return Omb;}
  
  /// Omega neutrino, renormalizes P(k) to keep sig8 fixed
  void setOmega_neutrino(double Omega_neutrino){Omnu = Omega_neutrino; TFmdm_set_cosm(); power_normalize(sig8); 	calc_interp_dist(); cosmo_set = none;
  }
  double getOmega_neutrino() const {return Omnu;}
  
  /// Number of neutrino species, renormalizes P(k) to keep sig8 fixed
  void setNneutrino(double Nneutrino){Nnu = Nneutrino; TFmdm_set_cosm(); power_normalize(sig8); cosmo_set = none;}
  double getNneutrino() const {return Nnu;}
  
  /// Dark energy equation of state parameter p/rho = w + w_1 (1+z)
  void setW(double w){ww = w; 	calc_interp_dist(); cosmo_set = none;}
  double getW() const {return ww;}
  void setW1(double w){ww1 = w; 	calc_interp_dist(); cosmo_set = none;}
  double getW1() const {return ww1;}
  
  /// Running of primordial spectral index, P(k)_primordial \propto pow(k/h,n+dndlnk*log(k)), renormalizes P(k) to keep sig8 fixed
  void setdndlnk(double w){dndlnk = w; power_normalize(sig8); cosmo_set = none;}
  double getdndlnk() const {return dndlnk;}
  
  /// Alternative to w for dark energy/ alt. grav. structure evolution
  void setgamma(double gamm){gamma = gamm; cosmo_set = none;}
  double getgamma() const {return gamma;}
  // If physical = 1 all Omega are Omega*h^2 otherwise they have the usual definitions.
  // 2 gamma parameterization is used for dark energy
  void setDEtype(short tt){darkenergy = tt; cosmo_set = none;}
  short getDEtype() const {return darkenergy;}
  
  void setSigma8(double my_sig8){power_normalize(my_sig8); cosmo_set = none;}
  double getSigma8() const {return sig8;}
  
  
  void dzdangDist(double D,double z[],double dzdD[]);
  
  double totalMassDensityinHalos(int t,double alpha,double m_min,double z,double z1,double z2);
  
  /// set interpolation range and number of points
  void setInterpolation(double z_interp, std::size_t n_interp);
  
  /// The lensing critical density in Msun / Mpc^2
  double SigmaCrit(double zlens,double zsource) const;
  
protected:
  void SetConcordenceCosmology(CosmoParamSet cosmo_p);
  
  CosmoParamSet cosmo_set;
  
  // structure rappers to make integration thread safe
  struct drdz_struct{
    drdz_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.drdz(x);}
    COSMOLOGY const &cos;
  };
  struct drdzdark_struct{
    drdzdark_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.drdz_dark(x);}
    COSMOLOGY const &cos;
  };
  struct adrdz_struct{
    adrdz_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.drdz(x)/x;}
    COSMOLOGY const &cos;
  };
  struct adrdzdark_struct{
    adrdzdark_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.drdz_dark(x)/x;}
    COSMOLOGY const &cos;
  };

  struct normL_struct{
    normL_struct(COSMOLOGY &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.normL(x);}
    COSMOLOGY &cos;
  };

  double ddrdzdOmo(double x) const;
  double ddrdzdw(double x) const;
  double ddrdzdw1(double x) const;
  
  struct ddrdzdOmo_struct{
    ddrdzdOmo_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.ddrdzdOmo(x);}
    COSMOLOGY const &cos;
  };
  struct ddrdzdw_struct{
    ddrdzdw_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.ddrdzdw(x);}
    COSMOLOGY const &cos;
  };
  struct ddrdzdw1_struct{
    ddrdzdw1_struct(COSMOLOGY const &cosmo):cos(cosmo){};
    double operator()(double x) const {return cos.ddrdzdw1(x);}
    COSMOLOGY const &cos;
  };

  bool init_structure_functions;
  
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
  /// Dark energy equation of state parameter p/rho = ww + ww_1 (1+z)
  double ww;
  /// Dark energy equation of state parameter p/rho = ww + ww_1 (1+z)
  double ww1;
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
  void polintD(double xa[], double ya[], int n, double x, double *y, double *dy) const;
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
  
  void setinternals();
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
  void calc_interp_dist();
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
  double NFW_M(double x,double cons,double M200);
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
