
#include <iostream>
using namespace std;

#ifndef pi
#define pi  3.141593
#endif

#ifndef Grav
#define Grav 4.7788e-20
#endif

#ifndef cosmo_declare


/** \ingroup cosmolib
 *
 * \brief The cosmology and all the functions required to calculated quantities based on the cosmology.
 */
class COSMOLOGY{
public:

    void SetConcordenceCosmology();
    void PrintCosmology(short physical = 0);
    double rcurve();
    double DRradius(double zo,double z,double pfrac);
    double DRradius2(double zo,double z);
    double Dgrowth(double z);
    double coorDist(double zo,double z);
    double angDist(double zo,double z);
    double lumDist(double zo,double z);
    double power_normalize(double sigma8);

    double power_linear(double k,double z);
    double powerCDMz(double k,double z);
    //double powerCDM(double k,double rt);
    double psdfdm(double z,double m);
    double stdfdm(double z,double m);
    double De(double rad);

    typedef double (COSMOLOGY::*pt2MemFunc)(double);

    /// accesser functions

	/// Hubble paremters in units of 100 km/s/Mpc
	void sethubble(double ht){ h = ht; TFmdm_set_cosm();}
	double gethubble(){return h;}

    /// Primordial spectral index
	void setindex(double nn){ n = nn;}
	double getindex(){ return n;}

    /// Omega matter
	void setOmega_matter(double Omega_matter){Omo = Omega_matter; TFmdm_set_cosm();}
	double getOmega_matter(){return Omo;}

    /// Omega lambda
	void setOmega_lambda(double Omega_lambda){Oml = Omega_lambda; TFmdm_set_cosm();}
	double getOmega_lambda(){return Oml;}

    /// Omega baryon
	void setOmega_baryon(double Omega_baryon){Omb = Omega_baryon; TFmdm_set_cosm();}
	double getOmega_baryon(){return Omb;}

    /// Omega neutrino
	void setOmega_neutrino(double Omega_neutrino){Omnu = Omega_neutrino; TFmdm_set_cosm();}
	double getOmega_neutrino(){return Omnu;}

	/// Number of neutrino species
	void setNneutrino(double Nneutrino){Nnu = Nneutrino; TFmdm_set_cosm();}
	double getNneutrino(){return Nnu;}

   /// Dark energy equation of state parameter p/rho = w + w_1 (1+z)
	void setW(double ww){w = ww;}
	double getW(){return w;}
	void setW1(double ww){w1 = ww;}
	double getW1(){return w1;}

	/// Running of primordial spectral index, P(k)_primordial \propto pow(k/h,n+dndlnk*log(k))
	void setdndlnk(double w){dndlnk = w;}
	double getdndlnk(){return dndlnk;}

	/// Alternative to w for dark energy/ alt. grav. structure evolution
	void setgamma(double gamm){gamma = gamm;}
	double getgamma(){return gamma;}
	// If physical = 1 all Omega are Omega*h^2 otherwise they have the usual definitions.
	// 2 gamma parameterization is used for dark energy
	void setDEtype(short tt){darkenergy = tt;}
	short getDEtype(){return darkenergy;}
    double getSigma8(){return sig8;}

    COSMOLOGY();
    ~COSMOLOGY();

private:

	/// Hubble paremters in units of 100 km/s/Mpc
  double h;
    /// Primordial spectral index
  double n;
    /// Omega matter
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
  double w1;
  //double Gamma;
   /// Running of primordial spectral index, P(k)_primordial \propto pow(k/h,n+dndlnk*log(k))
  double dndlnk;
  /// Alternative to w for dark energy/ alt. grav. structure evolution
  double gamma;
  // If physical = 1 all Omega are Omega*h^2 otherwise they have the usual definitions.
  //short physical;
  // 2 gamma parameterization is used for dark energy
  short darkenergy;

/* table for growth parameter */
    double *a;
    double *growth;
    int Ntable;

    double A;
    double sig8;  /* do not access these normalization outside */

    double powerEH(double k,double z);
    double powerEHv2(double k);
    double powerloc(double k,double z);
    double npow(double k);
    double normL(double lgk);
    double radiusm(double x);
    double radiusm_dark(double x);
    double gradius(double R,double rd);

    double nintegrateDcos(pt2MemFunc func,double a,double b,double tols);
    double trapzdDcoslocal(pt2MemFunc func, double a, double b, int n);

    /* in powerEH.c */
    short TFmdm_set_cosm_change_z(double redshift);
    short TFmdm_set_cosm();
    double TFmdm_onek_mpc(double kk);
    double TFmdm_onek_hmpc(double kk);

    // The parameters used in Eisenstein & Hu power spectrum

    	double alpha_gamma;	/* sqrt(alpha_nu) */
    	double alpha_nu;	/* The small-scale suppression */
    	double beta_c;		/* The correction to the log in the small-scale */
    	double num_degen_hdm;	/* Number of degenerate massive neutrino species */
    	double f_baryon;	/* Baryon fraction */
    	double f_bnu;		/* Baryon + Massive Neutrino fraction */
    	double f_cb;		/* Baryon + CDM fraction */
    	double f_cdm;		/* CDM fraction */
    	double f_hdm;		/* Massive Neutrino fraction */
    	double growth_k0;	/* D_1(z) -- the growth function as k->0 */
    	double growth_to_z0;	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
    	double k_equality;	/* The comoving wave number of the horizon at equality*/
    	double obhh;		/* Omega_baryon * hubble^2 */
    	double omega_curv;	/* = 1 - omega_matter - omega_lambda */
    	double omega_lambda_z; /* Omega_lambda at the given redshift */
    	double omega_matter_z;	/* Omega_matter at the given redshift */
    	double omhh;		/* Omega_matter * hubble^2 */
    	double onhh;		/* Omega_hdm * hubble^2 */
    	double p_c;		    /* The correction to the exponent before drag epoch */
    	double p_cb;		/* The correction to the exponent after drag epoch */
    	double sound_horizon_fit;  /* The sound horizon at the drag epoch */
    	double theta_cmb;	/* The temperature of the CMB, in units of 2.7 K */
    	double y_drag;		/* Ratio of z_equality to z_drag */
    	double z_drag;		/* Redshift of the drag epoch */
    	double z_equality;	/* Redshift of matter-radiation equality */

};

typedef COSMOLOGY *CosmoHndl;


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
double dsigdM(double m);
double f4(double u);
double Deltao(double m);

/* in powerEHv2.c */
void TFset_parameters(double omega0hh, double f_baryon, double Tcmb);
double TFfit_onek(double k, double *tf_baryon, double *tf_cdm);
