#include <math.h>
#include <nrD.h>
#include <cosmo.h>
#include <assert.h>
#include "utilities.h"

static double omo, oml, hh;


double COSMOLOGY::powerCDMz(
		double k    /// scale in the Fourier space
		,double z   /// redshift
		){

  if(!init_structure_functions) setinternals();

  double kn=0.0,kl,powL,powNL,nin;
  double knl,knh,kll,klh,g,a,Omot,Omlt,AA,B,aa,b,V,go;
  int m=0;
  double omo,oml;

  omo=Omo;
  oml=Oml;

  go=2.5*omo/( pow(omo,4.0/7.0)-oml+(1+0.5*omo)*(1+oml/70) );

  a=1.0/(1+z);
  if(omo==1.0){g=go;
  }else{
    Omot=omo/(omo+oml*a*a*a-a*(omo+oml-1));
    Omlt=a*a*a*oml*Omot/omo;
    g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
  }

  kl=k;
  nin=1+npow(0.75*kl)/3.0;
  AA=0.482*pow(nin,-0.947);
  B=0.226*pow(nin,-1.778);
  aa=3.31*pow(nin,-0.244);
  b=0.862*pow(nin,-0.287);
  V=11.55*pow(nin,-0.423);

  powL=5.066e-2*a*a*kl*kl*kl*powerloc(kl,z);
  powL*=g*g/(go*go);
  powNL=powL*pow( (1+B*b*powL+pow(AA*powL,aa*b))
		    /(1+pow( g*g*g*pow(AA*powL,aa)/(V*sqrt(powL)),b) ), 1.0/b );
  kn=kl*pow(1.0+powNL,0.333333);

  kll=0.0;  klh=k;                                     /* bracket */
  knl=0.0;  knh=kn;

  while( k<=1.0e4*fabs(kn-k) && m < 500 ){

    kl=kll-(knl-k)*(kll-klh)/(knl-knh);
    nin=1+npow(0.75*kl)/3.0;
    AA=0.482*pow(nin,-0.947);
    B=0.226*pow(nin,-1.778);
    aa=3.31*pow(nin,-0.244);
    b=0.862*pow(nin,-0.287);
    V=11.55*pow(nin,-0.423);

    powL=5.066e-2*a*a*kl*kl*kl*powerloc(kl,z)*g*g/(go*go);
    powNL=powL*pow( (1+B*b*powL+pow(AA*powL,aa*b))
		    /(1+pow( g*g*g*pow(AA*powL,aa)/(V*sqrt(powL)),b) ), 1.0/b );
    kn=kl*pow(1.0+powNL,0.33333);

    if( (kn-k) < 0.0){
      kll=kl;
      knl=kn;
    }else{
      klh=kl;
      knh=kn;
    }
   ++m;
  }
  //printf("%i %e %e %e %e\n",m,kl,k,powL,powNL);    /*test line */
  return 19.739*powNL/(k*k*k*a*a);
}

/**  
 * \brief The scale factor, a = 1/(1+z), as a function of radius in Mpc
 */

double COSMOLOGY::scalefactor(double rad) const{
  double a[1];
  int nok,nbad;
  void dir(double,double [],double []);
  double omo,oml;

  assert(ww == -1);
  assert(ww1 == 0.0);
  omo=Omo;
  oml=Oml;

  hh = h;

  if(omo==1.0) return pow(1-0.5*h*rad/3.0e3,2);

  a[0]=1.0;

  if(rad<1.0e1*h){ return 1.0 - h*rad/3.0e3;
  }else{
	  odeintD(a-1,1,0.0,rad,1.0e-6,rad/5,rad/1000,&nok,&nbad,dir,bsstepD);
  }
  if(a[0] < 100){
	  std::cerr << "COSMOLOGY::scalefactor() should be updated to do high redshift!" << std::endl;
	  exit(1);
  }
  return a[0];
}

void dir(double r,double a[],double dadr[]){
  dadr[1] = -hh*sqrt( a[1]*(omo+oml*pow(a[1],3)+(omo+oml-1.0)*a[1]) )/3.0e3;
}


//void COSMOLOGY::dzdangDist(double D,double z[],double dzdD[]){
//	dzdD[1] = (1+z[1])/( drdz(1+z[1]) - angDist(0,z[1]) );
//}

/** 
 * \brief Logorithmic slope of the power spectrum
 */
double COSMOLOGY::npow(double k){
  double qt;

  qt = k*exp(2*Omb)/(Omo*h*h);

  return n-2.0+4.68*qt/(log(1+2.34*qt)*(1+2.34*qt))
    -0.5*qt*(3.89+qt*(5.1842e2+qt*(4.8831e2+8.1088e3*qt)))/( 1+qt*(3.89+qt*(2.5921e2+qt*(1.6277e2+2.0272e3*qt))) );
}

double COSMOLOGY::power_linear(double k,double z){
  if(!init_structure_functions) setinternals();

  if(z==0.0) return powerloc(k,0);
  return pow(Dgrowth(z)*(1+z),2)*powerloc(k,z);
}

/** 
 * \brief The linear power spectrum without growth factor
 * growth factor should be normalized to 1 at z=0
 */
double COSMOLOGY::powerloc(double k,double z){

  if(Omnu == 0.0) return powerEHv2(k);
  return powerEH(k,z);
}
/** 
 * \brief Set the linear normalization for the power spectrum.
 *
 * This function keeps the internal normalization parameters in sync.  The normalization
 * should not be changed in any other way.
 */
double COSMOLOGY::power_normalize(double sigma8){
  if(!init_structure_functions) setinternals();
    
  double powfactor;
  sig8=sigma8;
  A=1.0;
  Rtophat = 8;
  ztmp = 0;
  
  powfactor=9*Utilities::nintegrateF(normL_struct(*this),log(1.0e-3),log(1.0e4),1.0e-9)/(2*PI*PI);

//  powfactor=9*nintegrateDcos(&COSMOLOGY::normL,log(1.0e-3),log(1.0e4),1.0e-9)/(2*PI*PI);
  A=sigma8*sigma8/powfactor;
  
  cosmo_set = CosmoParamSet::none;
  return powfactor;
}

double COSMOLOGY::normL(double lgk){
  double R,win,k;

  k=exp(lgk);
  R=k*Rtophat/h;
  if(R<=1.0e-4){ win = (1-R*R/10)/3; 
  }else{win=(sin(R)/R - cos(R))/(R*R);}

  return k*k*k*powerloc(k,ztmp)*win*win;
}
/** 
 * \brief Variance within a spherical top-hat filter of size R (Mpc), \f$ S(R)=\sigma^2(R) \f$, at redshift z.
 *
 * The variance is found through directly integrating linear power spectrum.
 */
double COSMOLOGY::TopHatVarianceR(double R,double z){
  
  if(!init_structure_functions) setinternals();

	double ans;

	Rtophat = R;
	ztmp = z;

  ans=9*Utilities::nintegrateF(normL_struct(*this),log(1.0e-3),log(1.0e4),1.0e-9)/(2*PI*PI);
  
	//ans = 9*nintegrateDcos(&COSMOLOGY::normL,log(1.0e-3),log(1.0e4),1.0e-9)/(2*PI*PI);
	return pow(Dgrowth(z)*(1+z),2)*ans;
}
/** 
 * \brief Variance within a spherical top-hat filter at mass scale M at redshift z.
 *
 * The variance is found through directly integrating linear power spectrum.
 */
double COSMOLOGY::TopHatVarianceM(double M,double z){
  if(!init_structure_functions) setinternals();

	double R = pow(M/rho_crit_comoving(0)/Omo,1./3.);

	return TopHatVarianceR(R,z);
}

/** 
 * \brief The power spectrum from Eisinstein & Hu with neutrinos but no BAO
 */
double COSMOLOGY::powerEH(double k,double z){
  //CosmoHndl cosmo_old;
  //static double zloc=-100;
  double Trans;

  TFmdm_set_cosm_change_z(z);
  Trans=TFmdm_onek_mpc(k);

  return A*pow(k/h,n+dndlnk*log(k))*Trans*Trans/pow(h/3.0e3,3);
}


/** 
 * \brief This is the power spectrum from Eisinstein & Hu 
 * with BAO but no neutrinos  
 */

double COSMOLOGY::powerEHv2(double k){

  double Trans;
  double baryon_piece,cdm_piece;
 
  Trans=TFfit_onek(k, &baryon_piece, &cdm_piece);
  //printf("trans=%e A=%e h=%e n=%e\n",Trans,A,h,n);
  return A*pow(k/h,n+dndlnk*log(k))*Trans*Trans/pow(h/3.0e3,3);
}


double COSMOLOGY::CorrelationFunction(double radius,double redshift
                                      ,double k_max,double k_min){
  
  CorrFunctorType func(this,1.0,redshift);
  
  if(k_max < k_min) std::swap(k_min,k_max);

  
  double a = PI/radius/2;
    
  func.r = radius;
    
  double kmin = k_min;
  double kmax = std::min(a,k_max);
  double tmp,ans=0;
  do{
    tmp = Utilities::nintegrate<CorrFunctorType,double>(func,kmin,kmax,1.0e-4);
    ans += tmp;
    kmin = kmax;
    kmax = std::min(kmax + a,k_max);
  }while(fabs(tmp/ans) > 1.0e-4);
    
  return ans;
  
}


