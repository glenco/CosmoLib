
#ifndef pi
#define pi  3.141593
#endif

#ifndef Grav
#define Grav 4.7788e-20
#endif

#ifndef cosmo_declare
typedef struct cosmology{
  double h;
  double n;
  double A;
  double sig8;  /* do not access these normalization outside */
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
} COSMOLOGY;

typedef COSMOLOGY *CosmoHndl;

#define cosmo_declare
#endif

/*** in cosmo.c ***/
void SetConcordenceCosmology(const CosmoHndl cosmo);
void PrintCosmology(struct cosmology cos1);
int cosmo_compare(struct cosmology *cos1,struct cosmology *cos2);
void cosmo_copy(struct cosmology *cos1,struct cosmology *cos2);
double rcurve(const CosmoHndl cosmo);
double radiusm(double x);
double DRradius(double zo,double z,double pfrac,const CosmoHndl cosmo);
double DRradius2(double zo,double z,const CosmoHndl cosmo);
void ders(double z,double Da[],double dDdz[]);
double arctanh(double x);
double fmini(double a,double b);
double fmaxi(double a,double b);
double Dgrowth(double z,const CosmoHndl cosmo);
double psdfdm(double sig8,double z,double m,const CosmoHndl cosmo);
double stdfdm(double sig8,double z,double m,const CosmoHndl cosmo);
double dsigdM(double m,const CosmoHndl cosmo);
double Deltao(double m,const CosmoHndl cosmo);
double f4(double u);
double Ex(double x,CosmoHndl cosmo);
double nintegrateDcos(double (*func)(double),double a,double b,double tols);
double trapzdDcoslocal(double (*func)(double), double a, double b, int n);
double coorDist(double zo,double z,const CosmoHndl cosmo);
double angDist(double zo,double z,const CosmoHndl cosmo);
double lumDist(double zo,double z,const CosmoHndl cosmo);
double gradius(double R,double rd,const CosmoHndl cosmo);
/* in powerCDM.c */
double powerCDMz(double k,double z,const CosmoHndl cosmo);
double powerCDM(double k,double rt,const CosmoHndl cosmo);
double De(double rad,const CosmoHndl cosmo);
void dir(double r,double a[],double dadr[]);
double npow(double k,const CosmoHndl cosmo);
double powerloc(double k,double z,const CosmoHndl cosmo);
double power_normalize(double sigma8,CosmoHndl cosmo);
double normL(double lgk);
double powerEH(double k,double z,const CosmoHndl cosmo);
double powerEHv2(double k,const CosmoHndl cosmo);
/* in powerCDMdark.c */
double powerCDMdark(double k,double z,const CosmoHndl cosmo);
double radiusm_dark(double x);
double Ex_dark(double x,CosmoHndl cosmo);
double dEx_darkda(double x,const CosmoHndl cosmo);
double growth_dark1(double a,const CosmoHndl cosmo);
double growth_dark2(double a,const CosmoHndl cosmo);
double dg2dla(double lna);
double growth_dark(double a,const CosmoHndl cosmo);
void dgrowth_dark(double lna,double y[],double dydlna[]);
double growth_table(double a,double *a_table,double *g_table,int Ntable);
void make_growth_table(double z,double *a_table,double *g_table,int Ntable,const CosmoHndl cosmo);

/* in powerEH.c */
int TFmdm_set_cosm(double omega_matter, double omega_baryon, double omega_hdm,
	double degen_hdm, double omega_lambda, double hubble, double redshift);
double TFmdm_onek_mpc(double kk);
double TFmdm_onek_hmpc(double kk);

/* in powerEHv2.c */
void TFset_parameters(double omega0hh, double f_baryon, double Tcmb);
double TFfit_onek(double k, double *tf_baryon, double *tf_cdm);
