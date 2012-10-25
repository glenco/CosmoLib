/** ******************************************

  cosmo.cpp calculates some useful comological parameters

*******************************************
*******************************************/
#include <math.h>
#include <nrD.h>
#include <nr.h>
#include <nrutil.h>
#include <cosmo.h>
#include <utilities.h>
#include <algorithm>

#define JMAX 34
#define JMAXP (JMAX+1)
#define K 6
#define NRANSI
#define CRITD 2.49783e18    /* critical density / Ho^2 solar/Mpc */
#define CRITD2 2.7752543e11 /* critical density / h^2 M_sun/Mpc^3 */
#define Hubble_length 3.0e3 /* Mpc/h */

int kmax,kount;
double *xp,**yp,dxsav;
static double alph_static;  /* DR-distance parameter */
static double Omo_static, Oml_static, h_static;

using namespace std;

COSMOLOGY::COSMOLOGY(double omegam,double omegal,double hubble, double ww) :
	Omo(omegam), Oml(omegal), h(hubble), w(ww){
	n=1.0;
	//Gamma=0.0;
	Omnu=0;
	Nnu=3.0;
	dndlnk=0.0;
	gamma=0.55;

	darkenergy=1;

	/* if 2 gamma parameterization is used for dark energy */
	/* if 1 w,w_1 parameterization is used for dark energy */

	// set parameters for Eisenstein&Hu power spectrum

	TFmdm_set_cosm();
	power_normalize(0.812);

	// allocate step and weight for gauleg integration
	ni = 64;
	xf=new float[ni];
	wf=new float[ni];
	gauleg(0.,1.,xf-1,wf-1,ni);
	// construct table of log(1+z), time, and \delta_c for interpolation
	fill_linear(vlz,ni,0.,1.7);
	double dc;
	for(int i=0;i<ni;i++){
	  double z = -1. + pow(10.,vlz[i]);
	  double Omz=Omegam(z);
	  if(Omo<1 && Oml==0) dc=1.68647*pow(Omz,0.0185);
	  if(Omo+Oml==1) dc=1.68647*pow(Omz,0.0055);
	  vDeltaCz.push_back(dc/Dgrowth(z));
	  vt.push_back(time(z));
	}

}

COSMOLOGY::COSMOLOGY(){
	SetConcordenceCosmology();

	// allocate step and weight for gauleg integration
	ni = 64;
	xf=new float[ni];
	wf=new float[ni];
	gauleg(0.,1.,xf-1,wf-1,ni);
	// construct table of log(1+z), time, and \delta_c for interpolation
	fill_linear(vlz,ni,0.,1.7);
	double dc;
	for(int i=0;i<ni;i++){
	  double z = -1. + pow(10.,vlz[i]);
	  double Omz=Omegam(z);
	  if(Omo<1 && Oml==0) dc=1.68647*pow(Omz,0.0185);
	  if(Omo+Oml==1) dc=1.68647*pow(Omz,0.0055);
	  vDeltaCz.push_back(dc/Dgrowth(z));
	  vt.push_back(time(z));
	}

}

COSMOLOGY::~COSMOLOGY(){
  // std::cout << "deleting cosmology" << std::endl;
  delete[] xf;
  delete[] wf;
}

/** \ingroup cosmolib
 * \brief Sets cosmology to WMAP 2009 model.  This is done automatically in the constructor.
 */
void COSMOLOGY::SetConcordenceCosmology(){
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

	w=-1.0;
	w1=0.0;
	n=1.0;
	//Gamma=0.0;
	Omnu=0;
	Nnu=3.0;
	dndlnk=0.0;
	gamma=0.55;

	darkenergy=1;

	/* if 2 gamma parameterization is used for dark energy */
	/* if 1 w,w_1 parameterization is used for dark energy */

	// set parameters for Eisenstein&Hu power spectrum

	TFmdm_set_cosm();
	power_normalize(0.812);
}

/** \ingroup cosmolib
 * \brief Print cosmological parameters 
 */

void COSMOLOGY::PrintCosmology(short physical){

	cout << "h: " << h << "\n";
	cout << "n: " << n << "\n";
	cout << "dndlnk: " << dndlnk << "\n";
	cout << "A: " << A << "\n";
	cout << "sig8: " << sig8 << "\n";
	//cout << "Gamma: " << Gamma << "\n";

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
  else cout << "darkenery: "<< darkenergy << " w: " << w << " w1: " << w1 << "\n";
}

/** \ingroup cosmolib
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

/** \ingroup cosmolib
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

/** \ingroup cosmolib
 * \brief Curvature radius in Mpc
 */

double COSMOLOGY::rcurve(){
  if(Omo+Oml != 1.0){
    return Hubble_length/(h*sqrt(fabs(1-Omo-Oml)));  /** curviture scale **/
  }
  return 0;
}


double COSMOLOGY::DpropDz(double z){
	double a=1.0/(1.0+z);
	return a*drdz_dark(1/a);
}


/** \ingroup cosmolib
 * \brief Matter density parameter at redshift z
 */
double COSMOLOGY::Omegam(double z){
	return Omo*pow(1+z,3)/( Omo*pow(1+z,3)+(1-Omo-Oml)*pow(1+z,2)+Oml );
}

/** \ingroup cosmolib
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

/** \ingroup cosmolib
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

/** \ingroup cosmolib
 * \brief Linear growth factor normalized to 1 at z=0
 */

double COSMOLOGY::Dgrowth(double z){
  double a,Omot,Omlt,g,go;

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
/** \ingroup cosmolib
 * \brief Critical density in M_sun/Mpc^3
 */
double COSMOLOGY::rho_crit(double z){
	return CRITD2*h*h*( Omo*pow(1+z,3)+Oml-(Omo+Oml-1)*pow(1+z,2) );
}

/** \ingroup cosmolib
 * \brief The derivative of the comoving radial distance with respect to redshift in units of Ho^-1, x=1+z
 */
double COSMOLOGY::drdz_empty(double x){
 double temp;

  temp=Oml;
  if(temp<=0.0) return -1.0e30;                   // non-physical values
  return 1.0/sqrt(temp);
}
double COSMOLOGY::drdz(double x){
 double temp;

  temp=Omo*x*x*x+Oml-(Omo+Oml-1)*x*x;
  if(temp<=0.0) return -1.0e30;                   // non-physical values
  return 1.0/sqrt(temp);
}
double COSMOLOGY::adrdz(double x){
  return drdz(x)/x;
}

/** \ingroup cosmolib
 * \brief Same as drdz, but incorporates dark energy through w and w1.
 */
double COSMOLOGY::drdz_dark(double x){
	return 1.0 / sqrt(Omo*x*x*x+Oml*pow(x,3*(1+w+w1))*exp(-3*w1*(x-1)/x)-(Omo+Oml-1)*x*x);
}
double COSMOLOGY::adrdz_dark(double x){
  return drdz_dark(x)/x;
}

/** \ingroup cosmolib
 * \brief The coordinate distance in units Mpc.  This is the radial distance found by integrating 1/H(z).
 */
double COSMOLOGY::coorDist(double zo,double z){
	if( (w ==-1.0) && (w1 == 0.0) ) return nintegrateDcos(&COSMOLOGY::drdz,1+zo,1+z,1.0e-9)*Hubble_length/h;
	return nintegrateDcos(&COSMOLOGY::drdz_dark,1+zo,1+z,1.0e-9)*Hubble_length/h;
}

double COSMOLOGY::emptyDist(double zo,double z){
  return nintegrateDcos(&COSMOLOGY::drdz_empty,1+zo,1+z,1.0e-9)*Hubble_length/h;
}

/** \ingroup cosmolib
 * \brief Non-comoving radial distance in units Mpc.  This is coorDist only integrated with the scale factor a=1/(1+z).
 */
double COSMOLOGY::radDist(double zo,double z){
	if( (w ==-1.0) && (w1 == 0.0) ) return nintegrateDcos(&COSMOLOGY::adrdz,1+zo,1+z,1.0e-9)*Hubble_length/h;
	return nintegrateDcos(&COSMOLOGY::adrdz_dark,1+zo,1+z,1.0e-9)*Hubble_length/h;
}

/** \ingroup cosmolib
 * \brief The angular size distance in units Mpc
 *
 *  Converts angles to proper distance NOT comoving distance.
 */

double COSMOLOGY::angDist(double zo,double z){
  double Rcur;

  if(Omo+Oml==1) return coorDist(zo,z)/(1+z);
   /** curviture scale **/

  Rcur=Hubble_length/(h*sqrt(fabs(1-Omo-Oml)));
  if((Omo+Oml)<1.0) return Rcur*sinh(coorDist(zo,z)/Rcur)/(1+z);

  return Rcur*sin(coorDist(zo,z)/Rcur)/(1+z);
}

/** \ingroup cosmolib
 * \brief The luminosity distance in units Mpc
 */
double COSMOLOGY::lumDist(double zo,double z){
	return pow(1+z,2)*angDist(zo,z);
}

/** \ingroup cosmolib
 * \brief Incorporates curvature for angular size distance.
 */
double COSMOLOGY::gradius(
		double R    /// the curvature radius
		,double rd  /// the coordinate
		){

  if( fabs(Omo+Oml-1) < 1.0e-4) return rd;
  if((Omo+Oml)<1.0) return R*sinh(rd/R);

  return R*sin(rd/R);
}

/**  \ingroup cosmolib
  \brief Press-Schechter mass function in unit of logarithmic mass bin and mean density
*/
double COSMOLOGY::psdfdm(
		double z      /// redshift
		,double m      /// mass
		,int caseunit  /// if equal to 1 return the number density per unit mass
		){

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
    		  return -(sig8*dsigdM(m))/sig*(sqrt(2/M_PI)*sqrt(nu)*exp(-0.5*nu))/m*Omo*CRITD2;
    		  break;
    	  default:
    		  return -(sig8*dsigdM(m))/sig*(sqrt(2/M_PI)*sqrt(nu)*exp(-0.5*nu));
    }
}

/** \ingroup cosmolib
 * \brief Sheth-Tormen mass function
 */

double COSMOLOGY::stdfdm(
		double z      /// redshift
		,double m      /// mass
		,int caseunit  /// if equal to 1 return the number density per unit mass
		){

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
  		  return -(sig8*dsigdM(m))/sig*(AST*(1+1/pow(nu,pST))*sqrt(nu/2)*exp(-0.5*nu))/m*Omo*CRITD2;
  		  break;
  	  case 2:
  	          nu = m*sqrt(aST);
  	          return AST*(1.+1./pow(nu,pST))*sqrt(nu/2.)*exp(-0.5*nu)/sqrt(M_PI);
  		  break;
  	  default:
  		  return -(sig8*dsigdM(m))/sig*(AST*(1+1/pow(nu,pST))*sqrt(nu/2)*exp(-0.5*nu));
  }
}

/** \ingroup cosmolib
 *  \brief Power-law mass function
 */

double COSMOLOGY::powerlawdfdm(
		double z          /// redshift
		,double m         /// mass
		,double alpha     /// slope (almost 1/6 for LCDM and cluster-size haloes)
		,int caseunit     /// if equal to 1 return the number density per unit mass
		){
	double mstar = nonlinMass(z);
	double mr = m/mstar;
	double alpha0 = 1./3.;
	switch (caseunit){
		case 1:
			return Omo*CRITD2*0.6*alpha0/2.*sqrt(2/M_PI)*pow(mr,alpha)*exp(-0.5*pow(mr*pow((1+z),3./2.),alpha0*1.3))/m/m;
			break;
		default:
			return 0.6*alpha0/2.*sqrt(2/M_PI)*pow(mr,alpha)*exp(-0.5*pow(mr*pow((1+z),3./2.),alpha0*1.3))/m/m;
	}
}

/** \ingroup cosmolib
 * \brief The cumulative comoving number density of haloes with mass larger than m
 * (in solar masses/h) at redshift z; if the dummy argument a is given, the mass function
 * times m^a is integrated. If a is omitted, a default of zero is assumed. The dummy
 * argument t specifies which type of mass function is to be used, 0 PS, 1 ST or 2 power-law
 * (in this case also the slope alpha can be set as an additional parameter)
 */
double COSMOLOGY::haloNumberDensity(
		double m      /// minimum mass of halos
		, double z    /// redshift
		, double a    /// moment of mass function that is wanted
		, int t       /// mass function type: 0 Press-Schecter, 1 Sheth-Torman, 2 power-law
		,double alpha /// exponent of power law if t==2
		){

	double n=0.0;
	double lm1 = log(m);
	double lm2 = 2.3*100.; // upper limit in natural log: I think 10^100 M_sun/h should be enough
	double lx,x,y,y1;

	for (int i=0;i<ni;i++){
		lx=(lm2-lm1)*xf[i]+lm1;
		x = exp(lx);
		switch (t){
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

/** \ingroup cosmolib
 * \brief The halo total surface mass density in haloes with mass larger than m_min
 * in the redshift bin [z1,z2] and projected to redhisft z
 * in solar mass / Mpc^2
 */
double COSMOLOGY::totalMassDensityinHalos(
		int t	       /// choice of mass function, 0 Press-Shechter, 1 Sheth-Tormen, 2 Power-law
		,double alpha  /// slope of power law if t==2
		,double m_min
		,double z
		,double z1
		,double z2){
  double n=0.0;
  double x,d,f,v,c;

   for (int i=0;i<ni;i++){
	   x=(z2-z1)*xf[i]+z1;
	   d=angDist(0.0,x)/Hubble_length*h;
	   f=1.0+x;
	   v=4.0*pi*d*d*DpropDz(x)*f*f*f;
	   c=haloNumberDensity(m_min,x,1,t,alpha);
	   n+=wf[i]*v*c;
  }

   double DL = angDist(0,z)*pi/180.;
   
   return n*(z2-z1)*pow(Hubble_length,3)/41253./DL/DL;

  //return haloNumberDensity(m_min,z,1,t,alpha)*h*h*pow(1+z,3)*angDist(z1,z2);
}

/** \ingroup cosmolib
 * \brief Number of haloes with mass larger than m (in solar masses/h) between
 * redshifts z1 and z2 per square degree
 * The flag t specifies which type of mass function is to be used, 0 PS or 1 ST
 */
double COSMOLOGY::haloNumberDensityOnSky (double m, double z1, double z2,int t, double alpha){
  double n=0.0;
  double x,d,f,v,c;

   for (int i=0;i<ni;i++){
	   x=(z2-z1)*xf[i]+z1;
	   d=angDist(0.0,x)/Hubble_length*h;
	   f=1.0+x;
	   v=4.0*pi*d*d*DpropDz(x)*f*f*f;
	   c=haloNumberDensity(m,x,0.0,t,alpha);
	   n+=wf[i]*v*c;
  }
   return n*(z2-z1)*pow(Hubble_length,3)/41253.;
}

/** \ingroup cosmolib
 * \brief Mass variance \f$ S(m)=\sigma^2(m) \f$.  This uses a fitting
 * formula for the CDM model which might not be perfectly accurate.  See
 * TopHatVarianceR() for an alternative.
 */
double COSMOLOGY::TopHatVariance(double m){
	double v = sig8*Deltao(m);
	return v*v;
}


/** \ingroup cosmolib
 * \brief Derivative of \f$ \sigma(m) \f$ in CDM model
 */
double COSMOLOGY::dsigdM(double m){
	double err;
	return dfridrDcos(&COSMOLOGY::Deltao,m,0.1*m,&err);
}

/** \ingroup cosmolib
 * \brief Virial overdensity
 */
double COSMOLOGY:: DeltaVir(
		double z         /// redshift
		,int caseunit    /// by default uses the Brayan and Norman fit, if equal to 1 uses the fit by Felix Stšhr
		){
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

/** \ingroup cosmolib
 * \brief Return the redshift  given \f$ \delta_c/D_+ \f$
 */
double COSMOLOGY::getZfromDeltaC(double dc){
	  if(dc>vDeltaCz[ni-1]) return -1+pow(10.,vlz[ni-1]);
	  if(dc<vDeltaCz[0]) return -1+pow(10.,vlz[0]);
	  int i = locate (vDeltaCz,dc);
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

/** \ingroup cosmolib
 * \brief Return the time from the Big Bang in Gyr given
 * \f$ \delta_c/D_+ \f$
 */
double COSMOLOGY::getTimefromDeltaC(double dc){
	  if(dc>vDeltaCz[ni-1]) return vt[ni-1];
	  if(dc<vDeltaCz[0]) return vt[0];
	  int i = locate (vDeltaCz,dc);
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

double COSMOLOGY::timeEarly(double a){
	double aEqual=8.3696e-05/Omo; // set by hand
	double r=aEqual/a;
	return 2.0/3.0*a*sqrt(a)/sqrt(Omo)*((1-2*r)*sqrt(1+r)+2*r*sqrt(r));
}

/** \ingroup cosmolib
 * \brief Return the time from the Big Bang in Gyr at a given redshift z
 */
double COSMOLOGY::time(double z){
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

double COSMOLOGY::nonlinMass(double z){
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


/** \ingroup cosmolib
 * \brief \f$ \sigma(m) \f$: the rms top-hat power in standard CDM model, normalized to sig8
 *
 * Warning! Uses global variables Omo_static,Oml_static and h_static. Not perfectly accurate.
 *
double Deltao(double m){
  double dc;
  dc=1.68647;
  if(Omo_static<1 && Oml_static==0) dc*=pow(Omo_static,0.0185);
  if(Omo_static+Oml_static==1) dc*=pow(Omo_static,0.0055);
  //return dc*pow(m/1.0e10,-0.25);
  return f4(6.005e14*pow(h_static*Omo_static,3))/f4(m*h_static*h_static*h_static*h_static*Omo_static*Omo_static);
}

*/
/** \ingroup cosmolib
 * \brief \f$ \sigma(m) \f$: the rms top-hat power in standard CDM model, normalized to sig8
 *
 */
double COSMOLOGY::Deltao(double m){
  double dc;
  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omo,0.0185);
  if(Omo+Oml==1) dc*=pow(Omo,0.0055);
  /*  return dc*pow(m/1.0e10,-0.25); */
  return f4(6.005e14*pow(h*Omo,3))/f4(m*h*h*h*h*Omo*Omo);
}

double f4(double u){
  return 8.6594e-12*pow(u,0.67)*pow( 1+pow( 3.5*pow(u,-0.1) +1.628e9*pow(u,-0.63),0.255) ,3.92157);
}

/** \ingroup cosmolib
 * \brief Halo bias, uses formalism by Mo-White
 *
 * by default it gives the halo bias by Mo-White
 * t=1 returns the Sheth-Tormen 99 while
 * setting t=2 the Sheth-Mo-Tormen 2001 bias
 */
double COSMOLOGY::halo_bias (
		double m       /// halo mass in solar masses
		,double z      /// redshift
		,int t         /// (0) Mo & White bias (1) Sheth-Tormen 99 (2) Sheth-Mo-Tormen 2001
		){
  double dc,Omz,sig,Dg;
  
  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);

  Omz=Omegam(z);

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  double nu2 = pow(dc/Dg/sig,2);
  double nu = dc/Dg/sig;
  
  switch (t){
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

/// The radius to which a halo must shrink to be 200 times as dense as the average density of the universe.
double COSMOLOGY::R200(double z,double mass){
	return pow(3*mass/800/pi/rho_crit(z),1.0/3.);
}

/***************************************************************/
/*** isolated integrator used for cosmological calculations ***/
/***************************************************************/
double COSMOLOGY::nintegrateDcos(pt2MemFunc func, double a,double b,double tols)
{
   double ss,dss;
   double s2[JMAXP],h2[JMAXP+1];
   int j;
   double ss2=0, ssum2=0;

   h2[1]=1.0;
   for (j=1;j<=JMAX;j++) {
     s2[j]=trapzdDcoslocal(func,a,b,j,&ss2,&ssum2);
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

double COSMOLOGY::trapzdDcoslocal(pt2MemFunc func, double a, double b, int n, double *s2, double *sum2)
{
   double x,tnm,del;
   int it,j;

   if (n == 1) {
	return (*s2=0.5*(b-a)*( (this->*func)(a) +(this->*func)(b) ));
   } else {
	for (it=1,j=1;j<n-1;j++) it <<= 1;
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (*sum2=0.0,j=1;j<=it;j++,x+=del) *sum2 += (this->*func)(x);
	*s2=0.5*(*s2+(b-a)*(*sum2)/tnm);
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


