/*
 * halo.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: cgiocoli
 */
#include <halo.h>

#define CRITD 2.49783e18    /* critical density / Ho^2 solar/Mpc */
#define CRITD2 2.7752543e11 /* critical density / h^2 M_sun/Mpc^3 */

 /** \ingroup cosmolib
  * \brief Constructor initializing a NFW-halo
  */
HALO::HALO (
		COSMOLOGY *cos     /// pointer to a COSMOLOGY
		,double mass       /// halo mass in M_sun/h
		,double redshift   /// redshift
		)
: co(cos),m(mass),z(redshift){
	Set_Parameters();
}

void HALO::Set_Parameters(){
	Omz=co->Omegam(z);
	Omo=co->getOmega_matter();
	Oml=co->getOmega_lambda();
	sigma2M=co->TopHatVariance(m);
	deltac0=1.68647;
	if(Omo<1 && Oml==0) deltac0*=pow(Omz,0.0185);
	if(Omo+Oml==1) deltac0*=pow(Omz,0.0055);
}

HALO::~HALO () {}

/** \ingroup cosmolib
 * \brief Reset halo mass and redshift
 */
void HALO::reset(double mr,double zr){
	m = mr;
	z = zr;
	Set_Parameters();
}

/** \ingroup cosmolib
 * \brief Virial radius of the halo in Mpc
 */
double HALO:: getRvir(int caseunit){
	double d=co->DeltaVir(z,caseunit)*Omo*CRITD2/Omz;
	return co->gethubble()*pow( 3*m/(4*M_PI*d), 0.3333 )/(1+z);
}

/** \ingroup cosmolib
 * \brief Radius at which the enclosed density reach 200 times the critical value
 */
double HALO:: getR200(){
	return 1.63e-5*pow( m*Omz/Omo, 0.3333 )/( 1.0+z );
}

/** \ingroup cosmolib
 * \brief The median redshift at which the main halo progenitor assembles
 * 50% of the halo mass, if 0<f<1 is given, is returned the redshift at which
 * this fraction is assembled
 */
double HALO:: getFormationRedshift(double f){
	double sigma2fM=co->TopHatVariance(f*m);
	double alphaf = 0.815*exp(-2*f*f*f)/pow(f,0.707);
	double wmed = sqrt(2*log(alphaf+1));
	double deltacz = deltac0/co->Dgrowth(z)+wmed*sqrt(sigma2fM-sigma2M);
	return co->getZfromDeltaC(deltacz);
}

/** \ingroup cosmolib
 * \brief The median time at which the main halo progenitor assembles
 * 50% of the halo mass, if 0<f<1 is given, is returned the redshift at which
 * this fraction is assembled
 */
double HALO:: getFormationTime(double f){
	double sigma2fM=co->TopHatVariance(f*m);
	double alphaf = 0.815*exp(-2*f*f*f)/pow(f,0.707);
	double wmed = sqrt(2*log(alphaf+1));
	double deltacz = deltac0/co->Dgrowth(z)+wmed*sqrt(sigma2fM-sigma2M);
	return co->getTimefromDeltaC(deltacz);
}

/** \ingroup cosmolib
 * \brief The halo concentration according to Zhao et al. 2009 is returned,
 * if the integer 1 is given the concentration according to Munoz-Cuartas et al. 2011 is returned
 * if 2 the concentration according to Giocoli et al. 2012 is returned
 * is returned
 */
double HALO:: getConcentration(int caseunit,double alpha0){
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
	    		return 10.*pow(m/1.e+12,alpha0)*pow(h0/hz,2./3.);
	    		break;
 	    	default: // Zhao et al. 2009
			  t004=getFormationTime(0.04);
			  t0=co->time(z);
			  return 4*pow(1+pow(t0/3.75/t004,8.4),1./8.);
		}

		return 0;
}

/** \ingroup cosmolib
 * \brief The halo total mass density in haloes with mass larger than m is returned,
 * if the integer 1 is given total mass density in haloes rescaled to the background is returned
 * two more parameter can be set,
 *   the first that define the halo mass function to use
 *       0 Press-Shechter
 *       1 Sheth-Tormen
 *       2 Power-law
 *  the second set the slope of the power-law mass function
 */
double HALO :: totalMassDensityinHalos(int caseunit,int t,double alpha){
  switch (caseunit){
	    	case (1):
	    		return co->haloNumberDensity(m,z,1,t,alpha)/(CRITD2*co->Omegam(0.));
	    		break;
	    	default:
	    		return co->haloNumberDensity(m,z,1,t,alpha);
  }

}

