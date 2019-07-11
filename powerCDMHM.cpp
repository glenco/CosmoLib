#include <powerCDMHM.h>
#include <nrD.h>
#include <nr.h>
#include <nrutil.h>
#include "utilities.h"

/**
  This  code reconstructs  the  non-linear dark  matter power  spectrum
  using the halo model. It considers 1 and 2 Halo terms and scatter in
  the c-m relation. 
*/

#ifdef ENABLE_GSL

const double tiny = 1.e-4;
const double CRITDD = 2.7752543e+11;

// integration parameters global
int nn  = 128;
const size_t limit = 128; 
const int key = 1;
gsl_integration_workspace * work = gsl_integration_workspace_alloc (nn);    
double err;
int nn2 = 16;
bool Init = true;

COSMOLOGY *co;   // cosmological model
double k;        // wave number 1/[Mpc/h]
double l;        // angular wave number 1/[arcmin]
double z;        // redshift
double rhob;     // back ground matter density
double minmass;  // minum halo mass in the integration of the mass function
const double maxmass = 1.e+25;  ;
double sigmalnC; // log-normal scatter in the c-m
int cmRelation;  // model for the c-m relation see halo.cpp
double slopeCM;  // slope of the c-m

/**
 * \brief Normalized fourier transform of the NFW profile
 */
double ukNFW(double m, double c){
  HALOCalculator ha(co,m,z);
  double Rvir = ha.getRvir()*(1+z)*pow(co->gethubble(),2./3.);
  double rs=Rvir/c;
  double rhos=m/( 4.*M_PI*pow( rs,3 )*(log( 1.+c )-c/( 1.+c )));
  // compute analytic Fourier transform
  double krs=k*rs;
  double ckrs=c*krs;
  double c1krs=ckrs+krs;
  double Si=gsl_sf_Si( krs );
  double Ci=gsl_sf_Ci( krs );
  double Sic=gsl_sf_Si( c1krs );
  double Cic=gsl_sf_Ci( c1krs );
  double ssin=sin( krs );
  double csin=sin( ckrs );
  double ccos=cos( krs );
  double graffa=ssin*(Sic-Si)-(csin/c1krs)+ccos*(Cic-Ci);
  return  4.*M_PI*rhos*pow(rs,3)*graffa/m;
}

// log nornal distribution normalization function
/*
double lognormal (double lx, void *ip){
  return 1/sqrt(2*M_PI*gsl_pow_2(sigmalnC))*exp(-gsl_pow_2(lx)/2/gsl_pow_2(sigmalnC));
}
*/

struct Glob{
  double dp1,dp2;
  int dp3;
};

/// log normal distribution times f function
double ftimeslognormal (double lx, void *ip){
  struct Glob* g=static_cast<struct Glob*> (ip);
  // dp1 - mass
  // dp2 - mean concentration
  // dp3 - exponent of uk in the integral

  double c = exp(lx)*g->dp2;
  double uk = ukNFW(g->dp1,c);
  return pow(uk,g->dp3)/sqrt(2*M_PI*gsl_pow_2(sigmalnC))*exp(-gsl_pow_2(lx)/2/gsl_pow_2(sigmalnC));
}

/// halo mass function normalization function
double intNORMmassfunction (double m, void *ip){
  m = pow(10.,m);
  double fnu = co->stdfdm(z,m,1)/pow(co->gethubble(),2);
  return log(10.)*m*m/rhob*fnu;  
}

/// halo mass function times bias normalization function
double intNORMmassfunctionbias (double m, void *ip){
  m = pow(10.,m);
  double bias = co->halo_bias (m,z,1);
  double fnu = co->stdfdm(z,m,1)/pow(co->gethubble(),2);
  return log(10.)*m*m/rhob*fnu*bias;
}

/// integral of the 1-Halo Component function
double int1H (double m, void *ip){
  m = pow(10.,m);
  HALOCalculator ha(co,m,z);  
  double c = ha.getConcentration(cmRelation,slopeCM);
  double uk = ukNFW(m,c);
  struct Glob g; 
  // dp1 - mass
  // dp2 - mean concentration
  // dp3 - exponent of uk in the integral

  g.dp1 = m;
  g.dp2 = c;
  g.dp3 = 2;

  gsl_function intflognorm;
  intflognorm.function = &ftimeslognormal;
  intflognorm.params = &g;
  double res;
  gsl_integration_workspace * workc = gsl_integration_workspace_alloc (nn);    
  gsl_integration_qag (&intflognorm,-2,2,tiny,tiny,limit,key,workc,&res,&err);    
  gsl_integration_workspace_free(workc);  

  double fnu = co->stdfdm(z,m,1)/pow(co->gethubble(),2);
  // return log(10.)*uk*uk*m/rhob*m*m/rhob*fnu;  // < --- no   scatter in c
  return log(10.)*m*m/rhob*fnu*res*m/rhob;     // < --- with scatter in c
}

/// integral of the 2-Halo Component function
double int2H (double m, void *ip){
  m = pow(10.,m);
  HALOCalculator ha(co,m,z);  
  double c = ha.getConcentration(cmRelation,slopeCM);
  double bias = co->halo_bias (m,z,1);
  double uk = ukNFW(m,c);

  struct Glob g; 
  // dp1 - mass
  // dp2 - mean concentration
  // dp3 - exponent of uk in the integral

  g.dp1 = m;
  g.dp2 = c;
  g.dp3 = 1;
  
  gsl_function intflognorm;
  intflognorm.function = &ftimeslognormal;
  intflognorm.params = &g;
  double res;
  gsl_integration_workspace * workc = gsl_integration_workspace_alloc (nn);    
  gsl_integration_qag (&intflognorm,-2,2,tiny,tiny,limit,key,workc,&res,&err);    
  gsl_integration_workspace_free(workc);  

  double fnu = co->stdfdm(z,m,1)/pow(co->gethubble(),2);
  // return log(10.)*m*m/rhob*fnu*bias*uk;     // no     scatter in c
  return log(10.)*m*m/rhob*fnu*bias*res;       // with   scatter in c
}

double POWERCDMHM::weight (double z1, double z2){
  double w1 = co->coorDist(0.,z1);
  double w2 =  co->coorDist(0.,z2);
  // we assume flat universe
  return (w2-w1)/w2;
}

double POWERCDMHM::weight (double z0){
  double a0=1.0/(1.0+z0);
  int i=Utilities::locate<double> (wgf.ai, a0);
  i=std::min(std::max(i,0),nn-2);
  return (wgf.wi[i+1]-wgf.wi[i])/(wgf.ai[i+1]-wgf.ai[i])*
      (a0-wgf.ai[i])+wgf.wi[i];
}

/** 
 * \brief Initialize cosmological models and integration to compute the non linear matter power spectrum using the halo model
 */
POWERCDMHM::POWERCDMHM(COSMOLOGY *_co      /// pointer to cosmology
		       ,double _z          /// redshift
		       ,double _minmass    /// minimum mass in the integral
		       ,double _sigmalnC   /// log-normal scatter in the c-m relation as fixed halo mass
		       ,int _cmRelation     /// c-m relation model to use, see halo.cpp
		       ,double _slopeCM    /// for power law cm relation it sets the slope
		       ){
  //co(_co), z(_z), minmass(_minmass), sigmalnC(_sigmalnC){
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
  co = _co;
  co->power_normalize(co->getSigma8());
  z=_z;
  minmass = _minmass;
  sigmalnC = _sigmalnC;
  cmRelation=_cmRelation;
  slopeCM=_slopeCM;
  rhob = CRITDD*co->Omegam(0.);
  // normalization of the 2Halo Term - mass function times bias
  gsl_function intNORMmfbias;
  intNORMmfbias.function = &intNORMmassfunctionbias;
  gsl_integration_qag (&intNORMmfbias,log10(minmass),log10(maxmass),
		       tiny,tiny,limit,key,work,&Pk20,&err);    
  Pk20 = gsl_pow_2(Pk20);
  intPk1.function = &int1H;
  intPk2.function = &int2H;

  // allocate step and weight for gauleg integration
  xf=new float[nn2];
  wf=new float[nn2];
}

/** 
 * \brief Return the non linear matter power spectrum calculated using the halo model
 */
double POWERCDMHM::nonlinpowerCDMHM(double _k){
  k=_k;
  Pklin = co->power_linear(k*co->gethubble(),0.)/(1+z)/(1+z)*pow(co->gethubble(),3);
  gsl_integration_qag (&intPk1,log10(minmass),log10(maxmass),tiny,tiny,limit,key,work,&Pk1,&err);
  gsl_integration_qag (&intPk2,log10(minmass),log10(maxmass),tiny,tiny,limit,key,work,&Pk2,&err);
  Pk2 = Pk2*Pk2*Pklin/Pk20;
  return Pk1+Pk2;
}

/** 
 * \brief Return the 1Halo term of the non linear matter power spectrum calculated using the halo model
 */
double POWERCDMHM::nonlinpowerCDMHM1Halo(double _k){
  k=_k;
  gsl_integration_qag (&intPk1,log10(minmass),log10(maxmass),tiny,tiny,limit,key,work,&Pk1,&err);
  return Pk1;
}
/** 
 * \brief Return the 2Halo term of the non linear matter power spectrum calculated using the halo model
 */
double POWERCDMHM::nonlinpowerCDMHM2Halo(double _k){
  k=_k;
  Pklin = co->power_linear(k*co->gethubble(),0.)/(1+z)/(1+z)*pow(co->gethubble(),3);
  gsl_integration_qag (&intPk2,log10(minmass),log10(maxmass),tiny,tiny,limit,key,work,&Pk2,&err);
  Pk2 = Pk2*Pk2*Pklin/Pk20;
  return Pk2;
}

void POWERCDMHM:: Initweight(double zs){
  Utilities::fill_linear (wgf.ai,nn, 0.1, 1.0-0.01);
  wgf.wi.resize (nn);
  for (int i=0;i<nn;i++){
    double zd=1.0/wgf.ai[i]-1.0;
    wgf.wi[i]=weight(zd,zs);
  }
  Init = false;
}

double POWERCDMHM:: nonlinKAPPApowerCDMHM(double l,double zs){
  if(Init) Initweight(zs);
  double as=(zs>0.0)?1.0/(1.0+zs):tiny;
  gauleg(as,1.,xf-1,wf-1,nn2);
  double p=0.0;
  for (int i=0;i<nn2;i++){
    double at=xf[i];
    double zt=1.0/at-1.0;
    double z0=0.0;
    double wt = co->coorDist(z0,zt);
    double dj=at*at/co->drdz_dark(1/at);
    double wwf=(zs>0.0)?weight(zt,zs)/at:weight(zt)/at;
    double kk=l/wt;
    double pk=nonlinpowerCDMHM(kk/co->gethubble());
    p+=wwf*wwf*pk*wf[i]/dj*co->gethubble();
  }
  p*=9.0*co->Omegam(0.)*co->Omegam(0.)/4.0/2.7e10;
  return p;
} 

double POWERCDMHM:: nonlinKAPPApowerCDMHM1Halo(double l, double zs){
  if(Init) Initweight(zs);
  double as=(zs>0.0)?1.0/(1.0+zs):tiny;
  gauleg(as,1.,xf-1,wf-1,nn2);
  double p=0.0;
  for (int i=0;i<nn2;i++){
    double at=xf[i];
    double zt=1.0/at-1.0;
    double z0=0.0;
    double wt = co->coorDist(z0,zt);
    double dj=at*at/co->drdz_dark(1/at);
    double wwf=(zs>0.0)?weight(zt,zs)/at:weight(zt)/at;
    double kk=l/wt;
    double pk=nonlinpowerCDMHM1Halo(kk/co->gethubble());
    p+=wwf*wwf*pk*wf[i]/dj*co->gethubble();
  }
  p*=9.0*co->Omegam(0.)*co->Omegam(0.)/4.0/2.7e10;
  return p;
}

double POWERCDMHM:: nonlinKAPPApowerCDMHM2Halo(double l, double zs){
  if(Init) Initweight(zs);
  double as=(zs>0.0)?1.0/(1.0+zs):tiny;
  gauleg(as,1.,xf-1,wf-1,nn2);
  double p=0.0;
  for (int i=0;i<nn2;i++){
    double at=xf[i];
    double zt=1.0/at-1.0;
    double z0=0.0;
    double wt = co->coorDist(z0,zt);
    double dj=at*at/co->drdz_dark(1/at);
    double wwf=(zs>0.0)?weight(zt,zs)/at:weight(zt)/at;
    double kk=l/wt;
    double pk=nonlinpowerCDMHM2Halo(kk/co->gethubble());
    p+=wwf*wwf*pk*wf[i]/dj*co->gethubble();
  }
  p*=9.0*co->Omegam(0.)*co->Omegam(0.)/4.0/2.7e10;
  return p;
}

double POWERCDMHM:: linKAPPApowerCDMHM(double l, double zs){
  if(Init) Initweight(zs);
  double as=(zs>0.0)?1.0/(1.0+zs):tiny;
  gauleg(as,1.,xf-1,wf-1,nn2);
  double p=0.0;
  for (int i=0;i<nn2;i++){
    double at=xf[i];
    double zt=1.0/at-1.0;
    double z0=0.0;
    double wt = co->coorDist(z0,zt);
    double dj=at*at/co->drdz_dark(1/at);
    double wwf=(zs>0.0)?weight(zt,zs)/at:weight(zt)/at;
    double kk=l/wt;
    double pk=co->power_linear(kk,0.)*at*at*pow(co->gethubble(),3);
    p+=wwf*wwf*pk*wf[i]/dj;
  }
  p*=9.0*co->Omegam(0.)*co->Omegam(0.)/4.0/2.7e10;
  return p;
}

double POWERCDMHM::nonlinfitKAPPApowerCDMHM(double l, double zs){
  if(Init) Initweight(zs);
  double as=(zs>0.0)?1.0/(1.0+zs):tiny;
  gauleg(as,1.,xf-1,wf-1,nn2);
  double p=0.0;
  for (int i=0;i<nn2;i++){
    double at=xf[i];
    double zt=1.0/at-1.0;
    double z0=0.0;
    double wt = co->coorDist(z0,zt);
    double dj=at*at/co->drdz_dark(1/at);
    double wwf=(zs>0.0)?weight(zt,zs)/at:weight(zt)/at;
    double kk=l/wt;
    double pk=co->powerCDMz(kk,zt)*at*at*pow(co->gethubble(),3);
    p+=wwf*wwf*pk*wf[i]/dj;
  }
  p*=9.0*co->Omegam(0.)*co->Omegam(0.)/4.0/2.7e10;
  return p;
}


POWERCDMHM:: ~POWERCDMHM(){
  delete[] xf;
  delete[] wf;
};

#endif
