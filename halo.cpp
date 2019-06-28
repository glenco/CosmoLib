/*
 * halo.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: cgiocoli
 */
#include <halo.h>

 /** 
  * \brief Constructor initializing a NFW-halo
  */
HALOCalculator::HALOCalculator (
		COSMOLOGY *cos     /// pointer to a COSMOLOGY
		,double mass       /// halo mass in M_sun/h   TODO Carlo:  Is this really M_sun/h
		,double redshift   /// redshift
		)
: co(cos),m(mass),z(redshift){
	Set_Parameters();
}

void HALOCalculator::Set_Parameters(){
	Omz=co->Omegam(z);
	Omo=co->getOmega_matter();
	Oml=co->getOmega_lambda();
	sigma2M=co->TopHatVariance(m);
	deltac0=1.68647;
	if(Omo<1 && Oml==0) deltac0*=pow(Omz,0.0185);
	if(Omo+Oml==1) deltac0*=pow(Omz,0.0055);
}

HALOCalculator::~HALOCalculator () {}

/** 
 * \brief Reset halo mass and redshift
 */
void HALOCalculator::reset(double mr,double zr){
	m = mr;
	z = zr;
	Set_Parameters();
}

/** 
 * \brief Virial radius of the halo in physical Mpc
 */
double HALOCalculator:: getRvir(
		      int caseunit    /// by default uses the Brayan and Norman fit, if equal to 1 uses the fit by Felix and Stoer
		      ){
  double d=co->DeltaVir(z,caseunit)*co->rho_crit_comoving(z);
  return pow( 3*m/(4*M_PI*d), 0.3333 )/(1+z);
}

/** 
 * \brief Radius in physical Mpc at which the enclosed density reach 200 times the critical value at that redshift
 */
double HALOCalculator:: getR200(){
  double d=200.0*co->rho_crit_comoving(z);
  return pow( 3*m/(4*M_PI*d), 0.3333 )/(1+z);
}

/** 
 * \brief The median redshift at which the main halo progenitor assembles
 * 50% of the halo mass, if 0<f<1 is given, is returned the redshift at which
 * this fraction is assembled
 */
double HALOCalculator:: getFormationRedshift(double f){
	double sigma2fM=co->TopHatVariance(f*m);
	double alphaf = 0.815*exp(-2*f*f*f)/pow(f,0.707);
	double wmed = sqrt(2*log(alphaf+1));
	double deltacz = deltac0/co->Dgrowth(z)+wmed*sqrt(sigma2fM-sigma2M);
	return co->getZfromDeltaC(deltacz);
}

/** 
 * \brief The median time at which the main halo progenitor assembles
 * 50% of the halo mass, if 0<f<1 is given, is returned the redshift at which
 * this fraction is assembled
 */
double HALOCalculator:: getFormationTime(double f){
	double sigma2fM=co->TopHatVariance(f*m);
	double alphaf = 0.815*exp(-2*f*f*f)/pow(f,0.707);
	double wmed = sqrt(2*log(alphaf+1));
	double deltacz = deltac0/co->Dgrowth(z)+wmed*sqrt(sigma2fM-sigma2M);
	return co->getTimefromDeltaC(deltacz);
}

/** 
 * \brief The halo concentration according to Zhao et al. 2009 is returned,
 * if the integer 1 is given the concentration according to Munoz-Cuartas et al. 2011 is returned
 * if 2 the concentration according to Giocoli et al. 2012 is returned
 * is returned
 */
double HALOCalculator:: getConcentration(
		int caseunit     /// set relation used - (0) Zhao et al. 2009, (1) Munoz-Cuartas et al. 2011 (2) Giocoli et al. 2012 (3) power-law c-m relation
		,double alpha0   /// power-law slope if caseunit == 3
		){
		double w=0.029;
		double mu=0.097;
		double alpha=-110.001;
		double beta=2469.720;
		double gamma=16.885;
		double a,b,logc;
		double t004,t05,t0;
		double h0,hz;
		switch (caseunit){
	    	case (1): // Munoz-Cuartas et al. 2011
	    		a = w*z-mu;
	    		b=alpha/(z+gamma)+beta/(z+gamma)/(z+gamma);
	    		logc=a*log10(m*co->gethubble())+b;
	    		return pow(10.,logc);
	    		break;
	    	case (-1): // Munoz-Cuartas et al. 2011, as the previous but M is assumet in M/h
	    		a = w*z-mu;
	    		b=alpha/(z+gamma)+beta/(z+gamma)/(z+gamma);
	    		logc=a*log10(m)+b;
	    		return pow(10.,logc);
	    		break;
	    	case (2): // Giocoli et al. 2012
			  t004=getFormationTime(0.04);
			  t05=getFormationTime(0.5);
			  t0=co->time(z);
			  return 0.45*(4.23+pow(t0/t004,1.15)+pow(t0/t05,2.3));
			  break;
	    	case (3): // power-law c-m relation
		        alpha0 = alpha0/(1+z);
	    		h0 = co->gethubble();
	    		hz = h0/co->drdz(1+z);
	    		return 10.*pow(m*h0/1.e+12,alpha0)*pow(h0/hz,2./3.);
	    		break;
	    	case (-3): // power-law c-m relation, as the previous but M is assumed in M/h
		        alpha0 = alpha0/(1+z);
	    		h0 = co->gethubble();
	    		hz = h0/co->drdz(1+z);
	    		return 10.*pow(m/1.e+12,alpha0)*pow(h0/hz,2./3.);
	    		break;
 	    	default: // Zhao et al. 2009
			  t004=getFormationTime(0.04);
			  t0=co->time(z);
			  return 4*pow(1+pow(t0/3.75/t004,8.4),1./8.);
		}

		return 0;
}

/** \brief returns the fraction of halo+galaxy in stars according to from Moster et al. 2010ApJ...710..903M
 */
double HALOCalculator::MosterStellarMassFraction(
                                       double Mtotal /// total mass of halo+galaxy in solar masses
){
  static const double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
  // from Moster et al. 2010ApJ...710..903M
  
  return mo*pow(Mtotal/M1,gam1)
  /pow(1+pow(Mtotal/M1,be),(gam1-gam2)/be)/Mtotal;
}

