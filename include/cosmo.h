#include <iostream>
#include <vector>

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
  double h;
  double n;
  double Omo;
  double Oml;
  double Omb;
  double Omnu;
  double Nnu;
  double w;
  double w1;
  double Gamma;
  double dndlnk;
  double gamma;   
  short physical; /* if physical =1 all Omega are Omega*h^2 */
  short darkenergy; /* if 2 gamma parameterization is used for dark energy */

/* table for growth parameter */
    double *a;
    double *growth;
    int Ntable;

    void SetConcordenceCosmology();
    void PrintCosmology();
    double rcurve();
    double Omegam(double z);
    double DRradius(double zo,double z,double pfrac);
    double DRradius2(double zo,double z);
    double Dgrowth(double z);
    double psdfdm(double sig8,double z,double m,int caseunit=0);
    double stdfdm(double sig8,double z,double m,int caseunit=0);
    double numberDensity(double sig8,double m,double z,double a, int t);
    double number (double sig8,double m,double z1,double z2,int t);
    double coorDist(double zo,double z);
    double angDist(double zo,double z);
    double lumDist(double zo,double z);
    double gradius(double R,double rd);
    double power_normalize(double sigma8);
    double radiusm(double x);
    double radiusm_dark(double x);
    double Variance(double m);
    double DeltaV(double z,int caseunit=0);
    double powerEH(double k,double z);
    double powerEHv2(double k);
    double npow(double k);
    double powerloc(double k,double z);
    double power_linear(double k,double z);
    double powerCDMz(double k,double z);
    double powerCDM(double k,double rt);
    double De(double rad);
    double normL(double lgk);
    double getZfromDeltaC(double dc);
    double eH(double a);
    double time(double z);
    double getTimefromDeltaC(double z);
    typedef double (COSMOLOGY::*pt2MemFunc)(double);

    double nintegrateDcos(pt2MemFunc func,double a,double b,double tols);
    double trapzdDcoslocal(pt2MemFunc func, double a, double b, int n);

    // accesser functions
    double getSigma8(){return sig8;}

    COSMOLOGY();
    ~COSMOLOGY();

private:
    int nbin;  // number of points when binning to interpolate
    double A;
    double sig8;  // do not access these normalization outside
    std:: vector<double> vDeltaCz,vlz,vt;
    double DpropDz(double z);
    double dsigdM(double m);
    double timeEarly(double a);
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

/* in powerEH.c */
int TFmdm_set_cosm(double omega_matter, double omega_baryon, double omega_hdm,
	double degen_hdm, double omega_lambda, double hubble, double redshift);
double TFmdm_onek_mpc(double kk);
double TFmdm_onek_hmpc(double kk);

/* in powerEHv2.c */
void TFset_parameters(double omega0hh, double f_baryon, double Tcmb);
double TFfit_onek(double k, double *tf_baryon, double *tf_cdm);
