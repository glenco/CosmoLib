/*
 * nfw.cpp
 *
 *  These are routines related to the Navarro, Frenk & White profile.
 *  These routines are limited to properties of the profile and do not
 *  include properties that depend on the context.
 *
 *  These formuli are mostly taken from:
 *  Navarro, Frenk & White, 1997 ApJ 490, 493
 *
 *  Created on: Mar 28, 2012
 *      Author: bmetcalf
 */

#include <math.h>
#include <nrD.h>
#include <cosmo.h>

using namespace std;

double vg;

/// Circular velocity at R200 in km/s
double COSMOLOGY::NFW_V200(
		double M200   /// Mass
		,double R200   /// Radius
		){
	return lightspeed*sqrt(Grav*M200/R200);
}
/// Maximum circular velocity in km/s
double COSMOLOGY::NFW_Vmax(
		double cons    /// concentration = R_200/R_s
		,double M200   /// Mass
		,double R200   /// Radius
		){
	double f = log(1+cons) - cons/(1+cons);

	return sqrt(0.216*cons/f)*NFW_V200(M200,R200);
}
/// Circular velocity in km/s
double COSMOLOGY::NFW_Vr(
		double x       /// radius , r/R_200
		,double cons    /// concentration = R_200/R_s
		,double M200    /// Mass
		,double R200    /// Radius
		){
	double f = log(1+cons) - cons/(1+cons);

	return sqrt( (log(1+cons*x) -cons*x/(1+cons*x) ) /f/x)*NFW_V200(M200,R200);
}
/// central over-density of nfw halo
double COSMOLOGY::NFW_deltac(
		double cons    /// concentration = R_200/R_s
		){
	return 200*cons*cons*cons/3/(log(1+cons)-cons/(1+cons));
}
/// Concentration  of NFW given Vmax in km/s
double COSMOLOGY::NFW_Concentration(
		double Vmax    /// Maximum circular velocity
		,double M200   /// Mass
		,double R200   /// Radius
		){
	double funcforconcentration(double cons);

	double tmp=NFW_Vmax(1.0,M200,R200);

	vg = Vmax/NFW_V200(M200,R200);
	if( vg < 1){
		cout << "ERROR: Vmax is too small in COSMOLOGY::Concentration_nfw()" << endl;
		return 1.0;  //// !!!!!!! must change this !!!!!!
		exit(0);
	}
	if(Vmax > NFW_Vmax(1000,M200,R200) ){
	    cout << "ERROR: Vmax is too large in COSMOLOGY::Concentration_nfw() concentration will be over 1000" << endl;
	    exit(0);
	}
	if(Vmax < NFW_Vmax(2.175,M200,R200) ){
	    cout << "ERROR: Vmax is too small in COSMOLOGY::Concentration_nfw() concentration will be under 2.175" << endl;
	    exit(0);
	}

	return zbrentD(funcforconcentration,2.175,1000,1.0e-8);
}

double funcforconcentration(double cons){
	double f = log(1+cons) - cons/(1+cons);
	return vg*vg - 0.216*cons/f;
}
/// The density of an NFW profile in units of the critical density
double COSMOLOGY::NFW_rho(
		double cons    /// concentration = R_200/R_s
		,double x      /// radius , r/R_200
		){
	if(x<=0) return 0;
	return NFW_deltac(cons)/x/pow(1+x,2);
}
