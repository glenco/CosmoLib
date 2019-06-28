/** ******************************************

  cosmo.cpp calculates some useful comological parameters

*******************************************
*******************************************/
#include <math.h>
#include <nrD.h>
#include <nr.h>
#include <nrutil.h>
#include <utilities.h>
#include <cosmo.h>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <assert.h>
//#include <gsl/gsl_integration_glfixed_table_alloc.h>

#ifdef ENABLE_GSL
  #include <gsl/gsl_errno.h>
  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <gsl/gsl_math.h>
  #include <gsl/gsl_spline.h>
  #include <gsl/gsl_rng.h>
#endif

#include "utilities.h"

#define JMAX 34
#define JMAXP (JMAX+1)
#define K 6
#define NRANSI
#define Hubble_length 3.0e3 /* Mpc/h */

int kmax,kount;
double *xp,**yp,dxsav;
static double alph_static;  /* DR-distance parameter */
static double Omo_static, Oml_static;

std::ostream &operator<<(std::ostream &os, CosmoParamSet p){
  switch (p) {
    case WMAP5yr :
      os << "WMAP5yr";
      break;
    case Millennium :
      os << "Millennium";
      break;
    case Planck1yr :
      os << "Planck1yr";
      break;
    case Planck :
      os << "Planck";
      break;
    case  BigMultiDark :
      os << "BigMultiDark";
      break;
    default:
      os << "????";
      break;
  }
  
  return os;
}

using namespace std;

COSMOLOGY::COSMOLOGY(double omegam,double omegal,double hubble, double w, double wa,bool justdistances) :
		h(hubble), Omo(omegam), Oml(omegal), ww(w) , ww1(wa){
	n=1.0;
	Omnu=0;
	Nnu=3.0;
	dndlnk=0.0;
	gamma=0.55;
  sig8 = 0.812;
      
  Omb = 0.02225/h/h;
      
	darkenergy=1;

	/* if 2 gamma parameterization is used for dark energy */
	/* if 1 w,w_1 parameterization is used for dark energy */

  if(!justdistances){
      setinternals();
  }else{
    // interpolate functions
    z_interp = 10.0;
    n_interp = 1024;
    calc_interp_dist();
  }
}

COSMOLOGY::COSMOLOGY(const COSMOLOGY &cosmo){
  
  h = cosmo.h;
  Omo = cosmo.Omo;
  Oml = cosmo.Oml;
  ww = cosmo.ww;
  n = cosmo.n;
  Omnu = cosmo.Omnu;
  Nnu = cosmo.Nnu;
  dndlnk = cosmo.dndlnk;
  gamma = cosmo.gamma;
  ww1 = cosmo.ww1;
  sig8 = cosmo.sig8;
  
  Omb = 0.02225/h/h;

  darkenergy = cosmo.darkenergy;

  setinternals();
}

COSMOLOGY::COSMOLOGY(CosmoParamSet cosmo_p){
  
	SetConcordenceCosmology(cosmo_p);
  
  setinternals();
  
  cosmo_set = cosmo_p;
}

COSMOLOGY::~COSMOLOGY(){
}

COSMOLOGY& COSMOLOGY::operator=(const COSMOLOGY &cosmo){
  
  if(this == &cosmo) return *this;

  h = cosmo.h;
  Omo = cosmo.Omo;
  Oml = cosmo.Oml;
  ww = cosmo.ww;
  n = cosmo.n;
  Omnu = cosmo.Omnu;
  Nnu = cosmo.Nnu;
  dndlnk = cosmo.dndlnk;
  gamma = cosmo.gamma;
  ww1 = cosmo.ww1;
  
  Omb = 0.02225/h/h;
  
  sig8 = cosmo.sig8;
  
  darkenergy = cosmo.darkenergy;
  
  setinternals();
  
  return *this;
}

void COSMOLOGY::setinternals(){
  
  init_structure_functions = true;

  TFmdm_set_cosm();
  power_normalize(sig8);
  // allocate step and weight for gauleg integration
  ni = 64;
  xf.resize(ni);
  wf.resize(ni);
  gauleg(0.,1.,&xf[0]-1,&wf[0]-1,ni);
  // construct table of log(1+z), time, and \delta_c for interpolation
  Utilities::fill_linear(vlz,ni,0.,1.7);
  double dc;
  for(int i=0;i<ni;i++){
    double z = -1. + pow(10.,vlz[i]);
    double Omz=Omegam(z);
    if(Omo<1 && Oml==0) dc=1.68647*pow(Omz,0.0185);
    if(Omo+Oml==1) dc=1.68647*pow(Omz,0.0055);
    vDeltaCz.push_back(dc/Dgrowth(z));
    vt.push_back(time(z));
  }
  
  // interpolate functions
  z_interp = 10.0;
  n_interp = 1024;
  calc_interp_dist();
  
}
/** 
 * \brief Sets cosmology to WMAP 2009 model.  This is done automatically in the constructor.
 */
void COSMOLOGY::SetConcordenceCosmology(CosmoParamSet cosmo_p){


  cosmo_set = cosmo_p;
	if(cosmo_p == WMAP5yr){
		// set cosmological parameters to standard WMAP 5r values
		// Komatsu et al. 2009, ApJ 180, 330
		// does not set power spectrum normalization
		// which needs to be done separately

		Omo=0.1358;
		Omb=0.02267;
		h=0.705;

		Omo/=h*h;
		Omb/=h*h;
		Oml=1.0-Omo;

		ww=-1.0;
		ww1=0.0;
		n=1.0;
		Omnu=0;
		Nnu=3.0;
		dndlnk=0.0;
		gamma=0.55;
        sig8=0.812;
		darkenergy=1;

  }else if(cosmo_p == Planck1yr){
    
    Oml = 0.6817;
    Omo = 1-Oml;
    h = 0.6704;
    
    Omb = 0.022032/h/h;
    
    ww=-1.0;
    ww1=0.0;
    n=1.0;
    Omnu=0;
    Nnu=3.0;
    dndlnk=0.0;
    gamma=0.55;
    sig8 = 0.829;
    
    darkenergy=1;
    
    /* if 2 gamma parameterization is used for dark energy */
    /* if 1 w,w_1 parameterization is used for dark energy */
    
  }else if(cosmo_p == Planck){
    
    // Final Planck cosmology Ade et al. 2015
    
    Omo = 0.308;
    Oml = 1-Omo;
    h = 0.678;
    
    Omb = 0.02225/h/h;
    
    ww=-1.0;
    ww1=0.0;
    n=0.968;
    Omnu=0;
    Nnu=3.0;
    dndlnk=0.0;
    gamma=0.55;
    sig8 = 0.8347;
    
    darkenergy=1;
    
	}else if(cosmo_p == Millennium){

		// The cosmology used in the Millennium simulations

		Omo=0.25;
		Omb=0.02267;
		h=0.73;

		Omb/=h*h;
		Oml=1.0-Omo;

		ww=-1.0;
		ww1=0.0;
		n=1.0;
		Omnu=0.0;
		Nnu=3.0;
		dndlnk=0.0;
		gamma=0.55;
    sig8 = 0.9;

		darkenergy=1;

  }else if (cosmo_p == BigMultiDark){
    
    Oml = 0.692885;
    Omo = 1-Oml;
    h = 0.677700;
    
    Omb = 0.022032/h/h;
    
    ww=-1.0;
    ww1=0.0;
    n=1.0;
    Omnu=0;
    Nnu=3.0;
    dndlnk=0.0;
    gamma=0.55;
    sig8 = 0.829;
    
    darkenergy=1;

  }
}

/** 
 * \brief Print cosmological parameters 
 */

void COSMOLOGY::PrintCosmology(short physical) const {

  cout << "parameters set :" << cosmo_set << "\n";
	cout << "h: " << h << "\n";
	cout << "n: " << n << "\n";
	cout << "dndlnk: " << dndlnk << "\n";
	cout << "A: " << A << "\n";
	cout << "sig8: " << sig8 << "\n";

	if(physical==0){
		cout << "Omo: " << Omo << "\n";
		cout << "Oml: " << Oml << "\n";
		cout << "Omb: " << Omb << "\n";
		cout << "Omnu: " << Omnu << "\n";
		cout << "Nnu: " << Nnu << "\n";
  }
  if(physical==1){
		cout << "Omo hh: " << Omo*h*h << "\n";
		cout << "Oml hh: " << Oml*h*h << "\n";
		cout << "Omb hh: " << Omb*h*h << "\n";
		cout << "Omnu hh: " << Omnu*h*h << "\n";
		cout << "Nnu: " << Nnu << "\n";
  }
  if(darkenergy==2) cout << "darkenery=" << darkenergy << " gamma=" << gamma << "\n";
  else cout << "darkenery: "<< darkenergy << " w: " << ww << " w1: " << ww1 << "\n";
}

/** 
 * \brief Print cosmological parameters
 */

std::string COSMOLOGY::print(short physical) const {
  
  std::string out;
  out = "h: " + std::to_string(h) + " ";
  out += " n: " + std::to_string(n) + " ";
  out += " dndlnk: " + std::to_string(dndlnk) + " ";
  out += " A: " + std::to_string(A) + " ";
  out += " sig8: " + std::to_string(sig8) + " ";
  
  if(physical==0){
    out += " Omo: " + std::to_string(Omo) + " ";
    out += " Oml: " + std::to_string(Oml) + " ";
    out += " Omb: " + std::to_string(Omb) + " ";
    out += " Omnu: " + std::to_string(Omnu) + " ";
    out += " Nnu: " + std::to_string(Nnu) + " ";
  }
  if(physical==1){
    out += " Omo hh: " + std::to_string(Omo*h*h) + " ";
    out += " Oml hh: " + std::to_string(Oml*h*h) + " ";
    out += " Omb hh: " + std::to_string(Omb*h*h) + " ";
    out += " Omnu hh: " + std::to_string(Omnu*h*h) + " ";
    out += " Nnu: " + std::to_string(Nnu) + " ";
  }
  if(darkenergy==2) out += " darkenery=" + std::to_string(darkenergy) + " gamma="
    + std::to_string(gamma) + " ";
  else out += " darkenery: "+ std::to_string(darkenergy) + " w: " + std::to_string(ww)
    + " w1: " + std::to_string(ww1) + " ";
  
  return out;
}

/** 
 * \brief see if cosmologies are identical
 */

int cosmo_compare(COSMOLOGY *cos1, COSMOLOGY *cos2){

  return 1-(cos1->gethubble() == cos2->gethubble())*(cos1->getindex() == cos2->getindex())
		  *(cos1->getSigma8() == cos2->getSigma8())
		  *(cos1->getOmega_matter() == cos2->getOmega_matter())*(cos1->getOmega_lambda() == cos2->getOmega_lambda())
		  *(cos1->getOmega_baryon() == cos2->getOmega_baryon())*(cos1->getOmega_neutrino() == cos2->getOmega_neutrino())
		  *(cos1->getW() == cos2->getW())*(cos1->getW1() == cos2->getW1())
		  *(cos1->getNneutrino() == cos2->getNneutrino());
}

/** 
 * \brief Copy a cosmology
 */

void cosmo_copy(CosmoHndl cos1, CosmoHndl cos2){

	//cos1->physical=cos2->physical;
	cos1->setOmega_matter(cos2->getOmega_matter());
	cos1->setOmega_lambda(cos2->getOmega_lambda());
	cos1->setOmega_baryon(cos2->getOmega_baryon());
	cos1->setOmega_neutrino(cos2->getOmega_neutrino());
	cos1->setNneutrino(cos2->getNneutrino());
	cos1->sethubble(cos2->gethubble());
	cos1->setW(cos2->getW());
	cos1->setW1(cos2->getW1());
	cos1->setindex(cos2->getindex());
	cos1->setdndlnk(cos2->getdndlnk());
	cos1->setgamma(cos2->getgamma());
	cos1->setDEtype(cos2->getDEtype());
	cos1->power_normalize(cos2->getSigma8());
}

/** 
 * \brief Curvature radius in Mpc
 */

double COSMOLOGY::rcurve() const{
  if(Omo+Oml != 1.0){
    return Hubble_length/(h*sqrt(fabs(1-Omo-Oml)));  /** curviture scale **/
  }
  return 0;
}


double COSMOLOGY::DpropDz(double z){
	double a=1.0/(1.0+z);
	return a*drdz(1/a);
}


/** 
 * \brief Matter density parameter at redshift z
 */
double COSMOLOGY::Omegam(double z) const{
	return Omo*pow(1+z,3)/( Omo*pow(1+z,3)+(1-Omo-Oml)*pow(1+z,2)+Oml );
}

/** 
 * \brief Comoving Dyer-Roeder angular size distance for lambda=0
 */
double COSMOLOGY::DRradius(double zo,double z,double pfrac){
  double zl,Omr,Omz,Da[2],Ho;
  int nok,nbad;

  Omo_static=Omo;
  Oml_static=Oml;

  Ho=h/Hubble_length;
  alph_static=1.5*(1-pfrac);

    /*    if(Oml_static !=0.0){
        printf("ERROR in DRradius: Oml_static != 0\n");
        return 0.0;
        }*/
    zl=1+z;

    if(zo==0 && (Oml==0 && alph_static==0)){
      if(Omo==1.0){
	return 2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);

      }else if(Oml==0){
	Omr=sqrt(1-Omo);
	Omz=sqrt(1+Omo*(zl-1));
	return zl*( (2-5*Omo) + Omz*( Omo*(3*zl+2)-2 )/(zl*zl) + 1.5*Omo*Omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-Omo,2));/**/
      }
    }else{
      if((z-zo)<0.001) return (z-zo)/(Ho*(1+zo)*(1+zo)*sqrt(1+Omo*zo+Oml*(pow(1+zo,-2)-1)) );
      Da[0]=0.0;
      Da[1]=-(1+zo)/(zl*zl*sqrt(1+Omo*z+Oml*(pow(zl,-2)-1)));

      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
    return 0;
}

/** 
 * \brief Comoving Dyer-Roeder angular size distance for lambda=0 and pfrac = 1 (all matter in particles)
 */
double COSMOLOGY::DRradius2(double zo,double z){
  double zl,Omr,Omz,Da[2],ds,dl,Ho;
  int nok,nbad;

  Omo_static=Omo;
  Oml_static=Oml;

    Ho=h/Hubble_length;
    /*    if(Oml !=0.0){
        printf("ERROR in DRradius: Oml != 0\n");
        return 0.0;
        }*/
    alph_static=0;
    zl=1+z;

    if(Omo==1.0){
      ds=2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
      if(zo==0.0){ return ds;
      }else{
	zl=zo+1;
	dl=2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
 	return   ds-dl*(1+z)/(1+zo);
 	/*return   ds*(1+zo)-dl*(1+z);*/
      }
      /*        return 2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho); */

    }else if(Oml==0){
      Omr=sqrt(1-Omo);
      Omz=sqrt(1+Omo*(zl-1));
      ds= zl*( (2-5*Omo) + Omz*( Omo*(3*zl+2)-2 )/(zl*zl) + 1.5*Omo*Omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-Omo,2));/**/

      if(zo==0.0){ return ds;
      }else{
	zl=zo+1;

	dl=zl*( (2-5*Omo) + Omz*( Omo*(3*zl+2)-2 )/(zl*zl) + 1.5*Omo*Omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-Omo,2));/**/

 	return   ds-dl*(1+z)/(1+zo);
	/* 	return   ds*(1+zo)-dl*(1+z); */
	/*return  (1+zo)*( ds-dl );*/
      }
	/*        return zl*( (2-5*Omo) + Omz*( Omo*(3*zl+2)-2 )/(zl*zl) + 1.5*Omo*Omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-Omo,2));*/

    }else{
      if((z-zo)<0.001) return (z-zo)/(Ho*(1+zo)*(1+zo)*sqrt(1+Omo*zo+Oml*(pow(1+zo,-2)-1)) );
      Da[0]=0.0;
      Da[1]=-1./(zl*zl*sqrt(1+Omo*z+Oml*(pow(zl,-2)-1)));

      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
}

void ders(double z,double Da[],double dDdz[]){
  dDdz[1]=Da[2];
  dDdz[2]=-3*Da[2]/(1+z) - (0.5*(Omo_static-2*Oml_static*pow(1+z,-3))*Da[2]+alph_static*Omo_static*Da[1]/(1+z) )/(1+z*Omo_static+(pow(1+z,-2)-1)*Oml_static);
}

/** 
 * \brief Linear growth factor normalized to 1 at z=0
 */

double COSMOLOGY::Dgrowth(double z) const{
  double a,Omot,Omlt,g,go;
  assert(init_structure_functions);


  a=1./(1+z);
  if(Omo==1 && Oml==0){
	  return a;
  }else{
	  go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) );
	  Omot=Omo/(Omo+Oml*a*a*a-a*(Omo+Oml-1));
	  Omlt=a*a*a*Oml*Omot/Omo;
	  g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
	  return a*g/go;
  }
}
/** 
 * \brief comoving critical density in M_sun/Mpc^3
 */
double COSMOLOGY::rho_crit_comoving(double z) const {
  return CRITD2*h*h*( Omo*pow(1+z,3)+Oml-(Omo+Oml-1)*pow(1+z,2) )/pow(1+z,3);
}

/** 
 * \brief The derivative of the comoving radial distance with respect to redshift in units of Ho^-1, x=1+z
 */

double drdz_wrapper(double x, void *params){
  return static_cast<CosmoHndl>(params)->COSMOLOGY::drdz(x);
}

double COSMOLOGY::drdz(double x) const{
 double temp;

  temp=Omo*x*x*x+Oml-(Omo+Oml-1)*x*x;
  if(temp<=0.0) return -1.0e30;                   // non-physical values
  return 1.0/sqrt(temp);
}

/** 
 * \brief Same as drdz, but incorporates dark energy through ww and ww1.
 */

double drdz_dark_wrapper(double x, void *params){
  return static_cast<CosmoHndl>(params)->COSMOLOGY::drdz_dark(x);
}
double COSMOLOGY::drdz_dark(double x) const{
	return 1.0 / sqrt(Omo*x*x*x+Oml*dark_factor(x)-(Omo+Oml-1)*x*x);
}
inline double COSMOLOGY::dark_factor(double x) const{
  if(ww == -1 && ww1 == 0) return 1;
  return pow(x,3*(1+ww+ww1))*exp(-3*ww1*(x-1)/x);
}

double COSMOLOGY::ddrdzdOmo(double x) const{
  double E=drdz_dark(x);
  return -0.5*(x*x*x-dark_factor(x))*E*E*E;
}
double COSMOLOGY::ddrdzdw(double x) const{
  double E=drdz_dark(x);
  return -0.5*3*Oml*log(x)*dark_factor(x)*E*E*E;
}

double COSMOLOGY::ddrdzdw1(double x) const{
  double E=drdz_dark(x);
  return -0.5*3*Oml*(log(x) - (x-1)/x)*dark_factor(x)*E*E*E;
}


/** 
 * \brief The coordinate distance in units Mpc.  This is the radial distance found by integrating 1/H(z).  This 
 *  is NOT the comoving angular size distance if the universe is not flat.
 */
double COSMOLOGY::coorDist(double zo,double z) const{
	// interpolation
	if(zo < z_interp && z < z_interp)
		return interp(coorDist_interp, z) - interp(coorDist_interp, zo);
	
  if(zo > z) return 0;
  
  if( (ww ==-1.0) && (ww1 == 0.0) ) return Utilities::nintegrateF(drdz_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;
  return Utilities::nintegrateF(drdzdark_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;
  //if( (ww ==-1.0) && (ww1 == 0.0) ) return nintegrateDcos(&COSMOLOGY::drdz,1+zo,1+z,1.0e-9)*Hubble_length/h;
	//return nintegrateDcos(&COSMOLOGY::drdz_dark,1+zo,1+z,1.0e-9)*Hubble_length/h;
}

double COSMOLOGY::d_coorDist_dOmo(double zo,double z) const{
  return Utilities::nintegrateF(ddrdzdOmo_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;
  //return nintegrateDcos(&COSMOLOGY::ddrdzdOmo,1+zo,1+z,1.0e-9)*Hubble_length/h;
}
double COSMOLOGY::d_coorDist_dw(double zo,double z) const{
  return Utilities::nintegrateF(ddrdzdw_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;
  //return nintegrateDcos(&COSMOLOGY::ddrdzdw,1+zo,1+z,1.0e-9)*Hubble_length/h;
}
double COSMOLOGY::d_coorDist_dw1(double zo,double z) const{
  return Utilities::nintegrateF(ddrdzdw1_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;
  //return nintegrateDcos(&COSMOLOGY::ddrdzdw1,1+zo,1+z,1.0e-9)*Hubble_length/h;
}


/** 
 * \brief Non-comoving radial distance in units Mpc also known as the lookback time.  This is coorDist only integrated with the scale factor a=1/(1+z).
 */
double COSMOLOGY::radDist(double zo,double z) const {
	if(zo < z_interp && z < z_interp)
		return interp(radDist_interp, z) - interp(radDist_interp, zo);
	
  if( (ww ==-1.0) && (ww1 == 0.0) ) return Utilities::nintegrateF(adrdz_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;
  return Utilities::nintegrateF(adrdzdark_struct(*this),1+zo,1+z,1.0e-9)*Hubble_length/h;

//	if( (ww ==-1.0) && (ww1 == 0.0) ) return nintegrateDcos(&COSMOLOGY::adrdz,1+zo,1+z,1.0e-9)*Hubble_length/h;
//	return nintegrateDcos(&COSMOLOGY::adrdz_dark,1+zo,1+z,1.0e-9)*Hubble_length/h;
}

/** 
 * \brief The angular size distance in units Mpc
 *
 *  Converts angles to proper distance NOT comoving distance.
 */

double COSMOLOGY::angDist(double zo,double z) const{
  double Rcur;

  if(Omo+Oml==1) return coorDist(zo,z)/(1+z);
   /** curviture scale **/

  Rcur=Hubble_length/(h*sqrt(fabs(1-Omo-Oml)));
  if((Omo+Oml)<1.0) return Rcur*sinh(coorDist(zo,z)/Rcur)/(1+z);

  return Rcur*sin(coorDist(zo,z)/Rcur)/(1+z);
}

/** 
 * \brief The bolometric luminosity distance in units Mpc
 */
double COSMOLOGY::lumDist(double zo,double z) const{
	return pow(1+z,2)*angDist(zo,z);
}

/** 
 * \brief The inverse of the coordinate distance in units Mpc, returning redshift. It works within interpolation range.
 */
double COSMOLOGY::invCoorDist(double d) const
{
  return invert(coorDist_interp, d);
}
/** 
 * \brief The inverse of the angular size distance in units Mpc, works within interpolation range.
 */
double COSMOLOGY::invComovingDist(double d) const
{
  if(Omo+Oml==1) return invCoorDist(d);
  double Rcurve = rcurve();
  if((Omo+Oml)<1.0) return invert(coorDist_interp, Rcurve*asinh(d/Rcurve) );
  return invert(coorDist_interp,Rcurve*asin(d/Rcurve) );
}
/** 
 * \brief The inverse of the coordinate distance in units Mpc, works within interpolation range.
 */

/** 
 * \brief The inverse of the radial distance in units Mpc, returning redshift. It works within interpolation range.
 */
double COSMOLOGY::invRadDist(double d) const
{
	return invert(radDist_interp, d);
}

/** 
 * \brief Incorporates curvature for angular size distance.
 */
double COSMOLOGY::gradius(
		double R    /// the curvature radius
		,double rd  /// the coordinate
		) const{

  if( fabs(Omo+Oml-1) < 1.0e-4) return rd;
  if((Omo+Oml)<1.0) return R*sinh(rd/R);

  return R*sin(rd/R);
}

/**  
*  \brief Press-Schechter halo mass function - default \f$ \frac{1}{\overline{\rho}_{comoving}}\frac{dN}{d\ln M} \f$ i.e. fraction of mass or
*  if caseunit==1 the \f$ \frac{dN}{dlogM} \f$
*/
double COSMOLOGY::psdfdm(
		double z      /// redshift
		,double m      /// mass
		,int caseunit  /// if equal to 1 return the comoving number density per unit mass
		){
  if(!init_structure_functions) setinternals();

  double dc,Omz,sig,Dg;

  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);

  Omz=Omegam(z);

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  /*printf("dc=%e dsigdM=%e sig=%e m=%e D=%e\n",dc,dsigdM(m),sig,m,Dg);*/
  double nu = pow(dc/Dg/sig,2);
  switch (caseunit){
  case 1:
	  return -Omo*rho_crit_comoving(0)*(sig8*dsigdM(m))/sig*(sqrt(2/M_PI)*sqrt(nu)*exp(-0.5*nu))/m;
	  break;
  default:
	  return -(sig8*dsigdM(m))/sig*(sqrt(2/M_PI)*sqrt(nu)*exp(-0.5*nu));
  }
}

/** 
 * \brief Sheth-Tormen mass function
 */

double COSMOLOGY::stdfdm(
		double z      /// redshift
		,double m      /// mass
		,int caseunit  /// if equal to 1 return the number density per unit mass TODO: Carlo, why aren't the other options explained?
		){

  if(!init_structure_functions) setinternals();

  double dc,Omz,sig,Dg;

  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);

  Omz=Omegam(z);

  dc=1.68647;

  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  /*return psdfdm(sig8,z,m)*( 1+0.9009*pow(Dg*sig/dc,-0.6) )
   *exp(0.1465*pow(dc/(Dg*sig),2) );*/
  double aST = 0.707;
  double pST = 0.3;
  double AST = 0.322;
  double nu = aST*pow(dc/Dg/sig,2);
  switch (caseunit){
  	  case 1:
  		  return -Omo*rho_crit_comoving(0)*(sig8*dsigdM(m))/sig*(AST*(1+1/pow(nu,pST))*sqrt(nu/M_PI/2)*exp(-0.5*nu))/m;
  		  break;
  	  case 2:
  	          nu = m*sqrt(aST);
  	          return AST*(1.+1./pow(nu,pST))*sqrt(nu/2.)*exp(-0.5*nu)/sqrt(M_PI);
  		  break;
  	  default:
  		  return -(sig8*dsigdM(m))/sig*(AST*(1+1/pow(nu,pST))*sqrt(nu/2)*exp(-0.5*nu));
  }
}

/** 
 *  \brief Power-law mass function
 */

double COSMOLOGY::powerlawdfdm(
		double z          /// redshift
		,double m         /// mass
		,double alpha     /// slope (almost 1/6 for LCDM and cluster-size haloes)
		,int caseunit     /// if equal to 1 return the number density per unit mass
		){
      
  if(!init_structure_functions) setinternals();

	double mstar = nonlinMass(z);
	double mr = m/mstar;
	double alpha0 = 1./3.;
	switch (caseunit){
		case 1:
			return Omo*rho_crit_comoving(0)*0.6*alpha0/2.*sqrt(2/M_PI)*pow(mr,alpha)*exp(-0.5*pow(mr*pow((1+z),3./2.),alpha0*1.3))/m/m;
			break;
		default:
			return 0.6*alpha0/2.*sqrt(2/M_PI)*pow(mr,alpha)*exp(-0.5*pow(mr*pow((1+z),3./2.),alpha0*1.3))/m/m;
	}
}

/** 
 * \brief The cumulative comoving number density of haloes with mass larger than m
 * (in solar masses) at redshift z; if the dummy argument a is given, the mass function
 * times m^a is integrated. If a is omitted, a default of zero is assumed. The dummy
 * argument type specifies which type of mass function is to be used, 0 PS, 1 ST or 2 power-law
 * (in this case also the slope alpha can be set as an additional parameter)
 */
double COSMOLOGY::haloNumberDensity(
		double m      /// minimum mass of halos in Msun
		, double z    /// redshift
		, double a    /// moment of mass function that is wanted
		, int type       /// mass function type: 0 Press-Schecter, 1 Sheth-Torman, 2 power-law
		,double alpha /// exponent of power law if type==2
		){
      
  if(!init_structure_functions) setinternals();

	double n=0.0;
	double lm1 = log(m);
	double lm2 = 2.3*100.; // upper limit in natural log: I think 10^100 M_sun/h should be enough
	double lx,x,y,y1;

	for (int i=0;i<ni;i++){
		lx=(lm2-lm1)*xf[i]+lm1;
		x = exp(lx);
		switch (type){
			case 1: // Sheth-Tormen
				y1=log(stdfdm(z,x,1));
				break;
			case 2: // power-law mass function
				y1=log(powerlawdfdm(z,x,alpha,1));
				break;
			default: // Press-Schechter
				y1=log(psdfdm(z,x,1));
				break;
		}
		y = exp(y1+(a+1.0)*lx);
		n+=wf[i]*y;
	}
	return n*(lm2-lm1);
}

/** 
 * \brief The halo total surface mass density in haloes with mass larger than m_min (solar masses)
 * in the redshift bin [z1,z2] and projected to redhisft z
 * in solar mass / proper Mpc^2
 *
 */
double COSMOLOGY::totalMassDensityinHalos(
		int type	       /// choice of mass function, 0 Press-Shechter, 1 Sheth-Tormen, 2 Power-law
		,double alpha      /// slope of power law mass function if type==2
		,double m_min      /// minimum halo mass in Msun
		,double z
		,double z1
		,double z2
		){
  double n = 0.0;
  double x,d,v,c;
  
  for (int i=0;i<ni;i++){
	  x=(z2-z1)*xf[i]+z1;
	  d=coorDist(0.0,x);
	  v=4.0*PI*d*d*drdz(1+x)*Hubble_length/h;
	  c=haloNumberDensity(m_min,x,1.0,type,alpha);
	  n+=wf[i]*v*c;
  }
  
  return  n*(z2-z1)/(4*PI*pow(angDist(0,z),2));


/*
	tmp_type = type;
	tmp_alpha = alpha;
	tmp_mass = m_min;
	tmp_a = 1.0;


	return nintegrateDcos(&COSMOLOGY::dNdz,z1,z2,1.0e-3)/(4*PI*pow(angDist(0,z),2));
  */
}


/** 
 * \brief Number of halos with mass larger than m (in solar masses)
 * between redshifts z1 and z2 per square degree.
 */
double COSMOLOGY::haloNumberDensityOnSky (
		double mass                 /// minimum mass
		,double z1                  /// lower redshift limit
		,double z2                  /// higher redshift limit
		,int type                   /// The flag type specifies which type of mass function is to be used, 0 PS or 1 ST
		,double alpha               /// slope of power law mass function if type==2
		){
  if(!init_structure_functions) setinternals();

  double n = 0.0;
  double x,d,v,c;

  for (int i=0;i<ni;i++){
	  x=(z2-z1)*xf[i]+z1;
	  d=coorDist(0.0,x);
	  v=4.0*PI*d*d*drdz(1+x)*Hubble_length/h;
	  c=haloNumberDensity(mass,x,0.0,type,alpha);
	  n+=wf[i]*v*c;
  }
 
 return n*(z2-z1)/41253.;
 
/*
  tmp_type = type;
  tmp_alpha = alpha;
  tmp_mass = mass;
  tmp_a = 0.0;

	std::cout << "totalMassDensityOnSky error = " << (ans-nintegrateDcos(&COSMOLOGY::dNdz,z1,z2,1.0e-3)/41253.)/ans << std::endl;

	return ans;

  return nintegrateDcos(&COSMOLOGY::dNdz,z1,z2,1.0e-3)/41253.;
*/
}
/**
 * \brief The number of halos in a buffered cone between redshift z1 and z2.
 *
 * A buffered cone is a cone with an extra perpendicular fixed physical distance added (area(z) = pi*(\theta D + buffer*(1+z))^2).
 * This geometry is useful for reducing edge effects which can be particularly bad a low redshift for small cones.
 */
double COSMOLOGY::haloNumberInBufferedCone (
		double mass                 /// minimum mass in Msun
		,double z1                  /// lower redshift limit
		,double z2                  /// higher redshift limit
		,double fov                 /// field of view of cone in steradians
		,double buffer              /// buffer length in physical Mpc (not comoving)
		,int type                   /// The flag type specifies which type of mass function is to be used, 0 PS or 1 ST
		,double alpha               /// slope of power law mass function if type==2
		){
  if(!init_structure_functions) setinternals();

	double n = 0.0;
	double z,d,v,c;

	for (int i=0;i<ni;i++){
	  z=(z2-z1)*xf[i]+z1;
	  d=sqrt(fov/PI)*coorDist(0.0,z) + buffer*(1+z);
	  v= PI*d*d*drdz(1+z)*Hubble_length/h;
	  c=haloNumberDensity(mass,z,0.0,type,alpha);
	  n+=wf[i]*v*c;
	}
	return n*(z2-z1);
}

/**
 * \brief Total mass contained in halos in a buffered cone between redshift z1 and z2.
 *
 * A buffered cone is a cone with an extra perpendicular fixed physical distance added (area(z) = pi*(\theta D + buffer*(1+z))^2).
 * This geometry is useful for reducing edge effects which can be particularly bad a low redshift for small cones.
 */

double COSMOLOGY::haloMassInBufferedCone (
		double mass                 /// minimum mass
		,double z1                  /// lower redshift limit
		,double z2                  /// higher redshift limit
		,double fov                 /// field of view of cone in steradians
		,double buffer              /// buffer length in physical Mpc (not comoving)
		,int type                   /// The flag type specifies which type of mass function is to be used, 0 PS or 1 ST
		,double alpha               /// slope of power law mass function if type==2
		){
      
  if(!init_structure_functions) setinternals();

	double n = 0.0;
	double z,d,v,c;

	for (int i=0;i<ni;i++){
	  z=(z2-z1)*xf[i]+z1;
	  d=sqrt(fov/PI)*coorDist(0.0,z) + buffer*(1+z);
	  v= PI*d*d*drdz(1+z)*Hubble_length/h;
	  c=haloNumberDensity(mass,z,1.0,type,alpha);
	  n+=wf[i]*v*c;
	}
	return n*(z2-z1);
}

/*
 * Total halo number (tmp_a=0) or mass (tmp_a=1) in solar masses in a redshift bin for the entire sphere.
 */
double COSMOLOGY::dNdz(double z){
  return 4.0*PI*pow(coorDist(0,z),2)*haloNumberDensity(tmp_mass,z,tmp_a,tmp_type,tmp_alpha)*drdz(1+z)*Hubble_length/h;
}

/** 
 * \brief Mass variance \f$ S(m)=\sigma^2(m) \f$.  This uses a fitting
 * formula for the CDM model which might not be perfectly accurate.  See
 * TopHatVarianceR() for an alternative.
 */
double COSMOLOGY::TopHatVariance(double m) const{
  assert(init_structure_functions);

	double v = sig8*Deltao(m);
	return v*v;
}


/** 
 * \brief Derivative of \f$ \sigma(m) \f$ in CDM model
 */
double COSMOLOGY::dsigdM(double m){
#ifdef ENABLE_GSL
  double result, error;
  gsl_function F;
  F.function = &Deltao_wrapper;
  F.params = this;

  gsl_deriv_central(&F,m,0.1*m,&result,&error);

  return result;
#else
	double err;
	return dfridrDcos(&COSMOLOGY::Deltao,m,0.1*m,&err);
#endif
}

/** 
 * \brief Virial overdensity
 */
double COSMOLOGY:: DeltaVir(
		double z         /// redshift
		,int caseunit    /// by default uses the Brayan and Norman fit, if equal to 1 uses the fit by Felix Stšhr
		) const{
	double omz = Omegam(z);
	double af,bf;
	switch (caseunit){
		case 1:
		// Felix Stšhr fit
			af=0.7076;
			bf=0.4403;
		    if (Oml<1e-4){
		    	af=0.1210;
		        bf=0.6756;
		    }
		    return 0.5*18*M_PI*M_PI*(1+(omz-1)*af+pow(omz,bf));
	    	break;
	    default:
	    // Bryan and Norman fit
	    	double x = omz -1;
	    	if(Oml<1e-4){
	    		return 18*M_PI*M_PI+60*x-32*x*x;
	    	}
	    	else{
	    		return 18*M_PI*M_PI+82*x-39*x*x;
	    	}
	}
}

/** 
 * \brief Return the redshift  given \f$ \delta_c/D_+ \f$
 */
double COSMOLOGY::getZfromDeltaC(double dc){
  if(!init_structure_functions) setinternals();

	  if(dc>vDeltaCz[ni-1]) return -1+pow(10.,vlz[ni-1]);
	  if(dc<vDeltaCz[0]) return -1+pow(10.,vlz[0]);
    int i = Utilities::locate (vDeltaCz,dc);
	  i = min (max (i,0), int (ni)-2);
	  double f=(dc-vDeltaCz[i])/(vDeltaCz[i+1]-vDeltaCz[i]);
	  if(i>1 && i<n-2){
		  // cubic interpolation
		  double a0,a1,a2,a3,f2;
		  f2 = f*f;
		  a0 = vlz[i+2] - vlz[i+1] - vlz[i-1] + vlz[i];
		  a1 = vlz[i-1] - vlz[i] - a0;
		  a2 = vlz[i+1] - vlz[i-1];
		  a3 = vlz[i];
		  double lzz =(a0*f*f2+a1*f2+a2*f+a3);
		  return -1+pow(10.,lzz);
	  }
	  // linear
	  else return -1 + pow(10,(f*vlz[i+1]+(1-f)*vlz[i]));
}

/** 
 * \brief Return the time from the Big Bang in Gyr given
 * \f$ \delta_c/D_+ \f$
 */
double COSMOLOGY::getTimefromDeltaC(double dc){
  if(!init_structure_functions) setinternals();

	  if(dc>vDeltaCz[ni-1]) return vt[ni-1];
	  if(dc<vDeltaCz[0]) return vt[0];
    int i = Utilities::locate (vDeltaCz,dc);
	  i = min (max (i,0), int (ni)-2);
	  double f=(dc-vDeltaCz[i])/(vDeltaCz[i+1]-vDeltaCz[i]);
	  if(i>1 && i<n-2){
		  // cubic interpolation
		  double a0,a1,a2,a3,f2;
		  f2 = f*f;
		  a0 = vt[i+2] - vt[i+1] - vt[i-1] + vt[i];
		  a1 = vt[i-1] - vt[i] - a0;
		  a2 = vt[i+1] - vt[i-1];
		  a3 = vt[i];
		  return a0*f*f2+a1*f2+a2*f+a3;
	  }
	  // linear
	  else return f*vt[i+1]+(1-f)*vt[i];
}

double COSMOLOGY::timeEarly(double a) const{
	double aEqual=8.3696e-05/Omo; // set by hand
	double r=aEqual/a;
	return 2.0/3.0*a*sqrt(a)/sqrt(Omo)*((1-2*r)*sqrt(1+r)+2*r*sqrt(r));
}

/** 
 * \brief Return the time from the Big Bang in Gyr at a given redshift z
 */
double COSMOLOGY::time(double z) const{
	double CfactorT = 0.102271;
        CfactorT = 1/(h*CfactorT);
	double a=1/(z+1.);
	double aEqual=8.3696e-05/Omo;  // equality
	double e=5.*aEqual;
	if(a>=e){
		double n=0,x,y;
		for (int i=0;i<ni;i++){
			x = e+(a-e)*xf[i];
			y = 1./x*drdz_dark(1/x);
			n+=wf[i]*y;
		}
		return (n*(a-e)+timeEarly(e))*CfactorT;
	}
	else{
		return CfactorT*timeEarly(a);
	}
}

double COSMOLOGY::nonlinMass(double z) const{
	  double Omz=Omegam(z);
	  double dc;
	  if(Omo<1 && Oml==0) dc=1.68647*pow(Omz,0.0185);
	  if(Omo+Oml==1) dc=1.68647*pow(Omz,0.0055);
	  double g = Dgrowth(z);
	  double lm[3]={2.0,6.5,15.0};
	  double s[3];
	  for (;;){
		for (size_t i=0;i<3;i++)
	      s[i]=TopHatVariance(pow(10.,lm[i] ));
	    if (*max_element(s,s+3)<dc/g)
	      for (size_t i=0;i<3;i++) lm[i]*=0.5;
	    else if (*min_element(s,s+3)>dc/g)
	      for (size_t i=0;i<3;i++) lm[i]*=2.0;
	    else
	    {
	      int i=0;
	      if (s[1]>=dc/g) i=1;
	      if (fabs(lm[2]-lm[0])<0.0001) break;
	      double lml=lm[i];
	      double lmu=lm[i+1];
	      lm[0]=lml;
	      lm[1]=0.5*(lml+lmu);
	      lm[2]=lmu;
	    }
	  }
	  double lm0 = (lm[0]+lm[1]+lm[2])/3.0;
	  return pow(10.,lm0);
}

/** 
 * \brief \f$ \sigma(m) \f$: the rms top-hat power in standard CDM model, normalized to sig8
 *
 */

double Deltao_wrapper(double m, void *params){
  return static_cast<CosmoHndl>(params)->COSMOLOGY::Deltao(m);
}

double COSMOLOGY::delta_c() const{
  double dc;
  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omo,0.0185);
  if(Omo+Oml==1) dc*=pow(Omo,0.0055);

  return dc;
}

double COSMOLOGY::Deltao(double m) const{
  /*double dc = delta_c();
  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omo,0.0185);
  if(Omo+Oml==1) dc*=pow(Omo,0.0055);
    return dc*pow(m/1.0e10,-0.25); */
  return f4(6.005e14*pow(h*Omo,3))/f4(m*h*h*h*h*Omo*Omo);
}

double  COSMOLOGY::f4(double u) const{
  return 8.6594e-12*pow(u,0.67)*pow( 1+pow( 3.5*pow(u,-0.1) +1.628e9*pow(u,-0.63),0.255) ,3.92157);
}

/** 
 * \brief Halo bias, uses formalism by Mo-White
 *
 * by default it gives the halo bias by Mo-White
 * type=1 returns the Sheth-Tormen 99 while
 * setting type=2 the Sheth-Mo-Tormen 2001 bias
 */
double COSMOLOGY::halo_bias (
		double m       /// halo mass in solar masses
		,double z      /// redshift
		,int type         /// (0) Mo & White bias (1) Sheth-Tormen 99 (2) Sheth-Mo-Tormen 2001
		){
      
  if(!init_structure_functions) setinternals();

  double dc,Omz,sig,Dg;
  
  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);

  Omz=Omegam(z);

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  double nu2 = pow(dc/Dg/sig,2);
  double nu = dc/Dg/sig;
  
  switch (type){
  case 1: // Sheth-Tormen 99
    return 1+((0.707*nu2-1)/dc)+(2*0.3/dc)/(1 + pow(0.707*nu2,0.3));
    break;
  case 2: // Sheth-Mo-Tormen 2001
    return 1+1/sqrt(0.707)/dc*(sqrt(0.707)*(0.707*nu*nu)+sqrt(0.707)*0.5*pow(0.707*nu*nu,0.4)
			       - pow(0.707*nu*nu,0.6)/
			       (pow(0.707*nu*nu,0.6) + 0.5*(1.-0.6)*(1-0.6/2.)));
    break;
  default: // Mo-White
    return 1+(nu2-1)/dc;
  }
}

/***************************************************************/
/*** isolated integrator used for cosmological calculations ***/
/***************************************************************/
double COSMOLOGY::nintegrateDcos(pt2MemFunc func, double a,double b,double tols) const
{
  
   double ss,dss;
   double s2[JMAXP],h2[JMAXP+1];
   int j;
   double ss2=0;

   h2[1]=1.0;
   for (j=1;j<=JMAX;j++) {
     s2[j]=trapzdDcoslocal(func,a,b,j,&ss2);
     if (j>=K) {
       COSMOLOGY::polintD(&h2[j-K],&s2[j-K],K,0.0,&ss,&dss);
       if(fabs(dss) <= tols*fabs(ss)) return ss;
     }
     h2[j+1]=0.25*h2[j];
   }
   cout << "s2= "; for(j=1;j<=JMAX;j++) cout << s2[j] << " ";
   cout << "\n";
   cout << "Too many steps in routine nintegrateDcos\n";
   throw std::runtime_error("COSMOLOGY::nintegrateDcos failure");
   return 0.0;
}

double COSMOLOGY::trapzdDcoslocal(pt2MemFunc func, double a, double b, int n, double *s2) const
{
  double x,tnm,del,sum;
   int it,j;

   if (n == 1) {
	return (*s2=0.5*(b-a)*( (this->*func)(a) +(this->*func)(b) ));
   } else {

	for (it=1,j=1;j<n-1;j++) it <<= 1;
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (this->*func)(x);
	*s2=0.5*(*s2+(b-a)*sum/tnm);

	return *s2;
   }
}
void COSMOLOGY::polintD(double xa[], double ya[], int n, double x, double *y, double *dy) const
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  
  dif=fabs(x-xa[1]);
  
  double *c = new double[n];
  double *d = new double[n];
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i-1]=ya[i];
    d[i-1]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i]-d[i-1];
      if ( (den=ho-hp) == 0.0) throw runtime_error("Error in routine polint");
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns] : d[ns-1]));
    --ns;
  }
  delete [] d;
  delete [] c;
}

double COSMOLOGY::nintegrateDcos(pt2MemFuncNonConst func, double a,double b,double tols)
{
  double ss,dss;
  double s2[JMAXP],h2[JMAXP+1];
  int j;
  double ss2=0;
  
  h2[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s2[j]=trapzdDcoslocal(func,a,b,j,&ss2);
    if (j>=K) {
      polintD(&h2[j-K],&s2[j-K],K,0.0,&ss,&dss);
      if(fabs(dss) <= tols*fabs(ss)) return ss;
    }
    h2[j+1]=0.25*h2[j];
	}
  cout << "s2= "; for(j=1;j<=JMAX;j++) cout << s2[j] << " ";
  cout << "\n";
  cout << "Too many steps in routine nintegrateDcos\n";
  return 0.0;
}

double COSMOLOGY::trapzdDcoslocal(pt2MemFuncNonConst func, double a, double b, int n, double *s2)
{
  double x,tnm,del,sum;
  int it,j;
  
  if (n == 1) {
    return (*s2=0.5*(b-a)*( (this->*func)(a) +(this->*func)(b) ));
  } else {
    
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (this->*func)(x);
    *s2=0.5*(*s2+(b-a)*sum/tnm);
    
    return *s2;
  }
}

double COSMOLOGY::dfridrDcos(pt2MemFunc func, double x, double b, double *err)
{
	double CON = 1.4;
	double CON2 = (CON*CON);
	double BIG = 1.0e30;
	long NTAB = 10;
	double SAFE = 2.0;
	int i,j;
	double errt,fac,bb,**a,ans;

	if (b == 0.0) cout << "b must be nonzero in dfridr." << endl;
	a=dmatrix(1,NTAB,1,NTAB);
	bb=b;
	a[1][1]=((this->*func)(x+bb)-(this->*func)(x-bb))/(2.0*bb);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		bb /= CON;
		a[1][i]=((this->*func)(x+bb)-(this->*func)(x-bb))/(2.0*bb);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	free_dmatrix(a,1,NTAB,1,NTAB);
	return ans;
}


void COSMOLOGY::setInterpolation(double my_z_interp, std::size_t my_n_interp)
{
  z_interp = my_z_interp;
  n_interp = my_n_interp;
  calc_interp_dist();
}

void COSMOLOGY::calc_interp_dist()
{

	// prepare vectors
	redshift_interp.resize(n_interp+1);
	coorDist_interp.resize(n_interp+1);
  radDist_interp.resize(n_interp+1);
	
	// step size going like square root
	double dz = z_interp/(n_interp*n_interp);

  size_t n = n_interp;
  double ztmp = z_interp;
  n_interp = 0;
  z_interp = 0;
	for(std::size_t i = 0; i <= n; ++i)
	{
		// make sure last value is exactly z_max
		double z = (i < n) ? (i*i*dz) : ztmp;
		
		redshift_interp[i] = z;
		coorDist_interp[i] = coorDist(0, z);
		radDist_interp[i] = radDist(0, z);
	}
  n_interp = n;
  z_interp = ztmp;
}

double COSMOLOGY::interp(const std::vector<double> & table, double z) const
{
	double dz = z_interp/(n_interp*n_interp);
	double di = sqrt(z/dz);
	std::size_t i = di;
	
	// interpolation range
	if(i < n_interp)
		return table[i] + (di-i)*(table[i+1]-table[i]);
	
	// maximum value
	if (i == n_interp)
		return table.back();
	
	// out of range
	throw std::out_of_range("redshift out of interpolation range");
}

double COSMOLOGY::invert(const std::vector<double>& table, double f_z) const
{
	if(f_z < 0)
		throw std::out_of_range("cannot invert negative function values");
	
	// find value in table
	std::size_t i = std::distance(table.begin(), std::lower_bound(table.begin(), table.end(), f_z));
	
	// out of range
	if(i > n_interp)
		throw std::out_of_range("function value out of interpolation range");
	
	double f_0 = (i > 0) ? table[i-1] : 0;
	double f_1 = table[i];
	
	double z_0 = (i > 0) ? redshift_interp[i-1] : 0;
	double z_1 = redshift_interp[i];
	
	return z_0 + ((f_z-f_0)/(f_1-f_0))*(z_1-z_0);
}

double COSMOLOGY::SigmaCrit(double zlens,double zsource) const {
  return angDist(0,zsource)/angDist(0,zlens)/angDist(zlens,zsource)/(4*PI*Grav);
}

