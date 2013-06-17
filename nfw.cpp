/*
 * nfw.cpp
 *
 *  These are routines related to the Navarro, Frenk & White profile.
 *  These routines are limited to properties of the profile and do not
 *  include properties that depend on the context.
 *
 *  These formula are mostly taken from:
 *  Navarro, Frenk & White, 1997 ApJ 490, 493
 *
 *  Created on: Mar 28, 2012
 *      Author: bmetcalf
 */

#include <math.h>
#include <nrD.h>
#include <cosmo.h>
#include <assert.h>

using namespace std;

double vg;

/// Circular velocity at R200 in km/s
double NFW_Utility::NFW_V200(
		double M200   /// Mass
		,double R200   /// Radius
		){
	return lightspeed*sqrt(Grav*M200/R200);
}
/// Maximum circular velocity in km/s
double NFW_Utility::NFW_Vmax(
		double cons    /// concentration = R_200/R_s
		,double M200   /// Mass
		,double R200   /// Radius
		){
	double f = log(1+cons) - cons/(1+cons);

	return sqrt(0.216*cons/f)*NFW_V200(M200,R200);
}
/// Circular velocity in km/s
double NFW_Utility::NFW_Vr(
		double x       /// radius , r/R_200
		,double cons    /// concentration = R_200/R_s
		,double M200    /// Mass
		,double R200    /// Radius
		){
	double f = log(1+cons) - cons/(1+cons);

	return sqrt( (log(1+cons*x) -cons*x/(1+cons*x) ) /f/x)*NFW_V200(M200,R200);
}
/// central over-density of nfw halo
double NFW_Utility::NFW_deltac(
		double cons    /// concentration = R_200/R_s
		){
	return 200*cons*cons*cons/3/(log(1+cons)-cons/(1+cons));
}
/// Concentration  of NFW given Vmax in km/s
double NFW_Utility::NFW_Concentration(
		double Vmax    /// Maximum circular velocity
		,double M200   /// Mass
		,double R200   /// Radius
		){
	double funcforconcentration(double cons);

	//double tmp=NFW_Vmax(1.0,M200,R200);

	vg = Vmax/NFW_V200(M200,R200);
	if( vg < 1){
		std::cout << "ERROR: Vmax is too small in NFW_Utility::Concentration_nfw()" << std::endl;
		exit(0);
		return 1.0;  //// !!!!!!! must change this !!!!!!
	}
	if(Vmax > NFW_Vmax(1000,M200,R200) ){
	    std::cout << "ERROR: Vmax is too large in NFW_Utility::Concentration_nfw() concentration will be over 1000" << std::endl;
	    exit(0);
	}
	if(Vmax < NFW_Vmax(2.175,M200,R200) ){
	    std::cout << "ERROR: Vmax is too small in NFW_Utility::Concentration_nfw() concentration will be under 2.175" << std::endl;
	    exit(0);
	}

	return zbrentD(&NFW_Utility::funcforconcentration,2.175,1000,1.0e-8);
}

float NFW_Utility::funcforconcentration(float cons){
	float f = log(1+cons) - cons/(1+cons);
	return vg*vg - 0.216*cons/f;
}
/// The density of an NFW profile in units of the critical density
double NFW_Utility::NFW_rho(
		double cons    /// concentration = R_200/R_s
		,double x      /// radius , r/R_200
		){
	if(x<=0) return 0;
	return NFW_deltac(cons)/x/pow(1+x,2);
}

/// Returns the concentration and radius of an NFW halo with the mass, half mass radius and Vmax provided
void NFW_Utility::match_nfw(
		float my_Vmax        /// Maximum circular velocity (km/s)
		,float my_R_half     /// Half mass radius (Mpc)
		,float my_mass       /// Mass (solar masses)
		,float *my_cons      /// output concentration
		,float *my_Rmax      /// Radius of halo,  Not necessarily R200 or Rvir.
		){

	//std::cout << "NFW_Utility Test: " << " mass: " << my_mass << " R_half: " << my_R_half << " Vmax: " << my_Vmax << std::endl;
	if(my_mass <= 0.0){
		*my_cons = 0;
		*my_Rmax = 0;
		return;
	}
	assert(my_Vmax > 0.0);
	assert(my_R_half > 0.0);

	mass = my_mass;
	R_half = my_R_half;
	Vmax = my_Vmax;

	if(nfwfunc(1.0e-4)*nfwfunc(1.0e4) > 0.0){
		ERROR_MESSAGE();
		std::cout << "ERROR: Vmax, R_half & mass are inconsistent!" << std::endl;
		exit(1);
	}

	*my_cons = zbrentD(&NFW_Utility::nfwfunc,1.0e-5,1.0e4,1.0e-8);
	//std::cout << "NFW_Utility Test: " << nfwfunc(1.0e-4)<< nfwfunc(*my_cons)<< nfwfunc(1.0e4) << " mass: " << mass << " R_half: " << R_half << " Vmax: " << Vmax << std::endl;
	*my_Rmax = Rmax(*my_cons,Vmax,mass);
	assert(*my_Rmax > my_R_half);
}

float NFW_Utility::nfwfunc(float cons){
	return 2*g_func(R_half*cons/Rmax(cons,Vmax,mass) ) - g_func(cons);
}

float NFW_Utility::Rmax(float cons,float Vmax,float mass){
	return Grav*mass*pow(0.216*lightspeed*cons/Vmax,2);
}
float NFW_Utility::g_func(float x){
	return log(1+x) - x/(1+x);
}

float NFW_Utility::zbrentD(MemFunc func, float x1, float x2, float tol){
	int iter,ITMAX = 100;
	float EPS = 3.0e-8;
	float a=x1,b=x2,c=x2,d,e,min1,min2;
	float fa=(this->*func)(a),fb=(this->*func)(b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
	        printf("fa=%e  fb=%e  x1=%e  x2=%e\n",fa,fb,x1,x2);
	        ERROR_MESSAGE();
	        std::cout << "Root must be bracketed in NFW_Utility::zbrentD" << std::endl;
	        exit(1);
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += ((xm) >= 0.0 ? fabs(tol1) : -fabs(tol1));
		fb=(this->*func)(b);
	}
	ERROR_MESSAGE();
	std::cout << "Maximum number of iterations exceeded in NFW_Utility::zbrentD" << std::endl;
	return 0.0;
}
