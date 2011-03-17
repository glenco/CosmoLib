/********************************************
  powerCDM.c calculates the nonlinear P(k,r)/a(r)^2
    example of how to use is in testpower.c
  *******************************************/
//#include "RecipesD/nrD.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <nrD.h>
#include <cosmo.h>
//#include "Cosmo/powerEH.c"
//#include "Cosmo/powerEHv2.c"

//extern struct cosmology cosmo;
static CosmoHndl cosmog;
/*** These cosmological parameters are redundant. if Gammma=0 
     it is not used else it is used and Omb appears nowhere in 
     this calculation  ***/

double powerCDMz(double k,double z,const CosmoHndl cosmo){

  double kn=0.0,kl,powL,powNL,nin;
  double knl,knh,kll,klh,g,a,Omot,Omlt,AA,B,aa,b,V,go;
  int m=0;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) );

  a=1.0/(1+z);
  if(Omo==1.0){g=go;
  }else{
    Omot=Omo/(Omo+Oml*a*a*a-a*(Omo+Oml-1));
    Omlt=a*a*a*Oml*Omot/Omo;
    g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
  }

  kl=k;
  nin=1+npow(0.75*kl,cosmo)/3.0;
  AA=0.482*pow(nin,-0.947);
  B=0.226*pow(nin,-1.778);
  aa=3.31*pow(nin,-0.244);
  b=0.862*pow(nin,-0.287);
  V=11.55*pow(nin,-0.423);

  powL=5.066e-2*a*a*kl*kl*kl*powerloc(kl,z,cosmo);
  powL*=g*g/(go*go);
  powNL=powL*pow( (1+B*b*powL+pow(AA*powL,aa*b))
		    /(1+pow( g*g*g*pow(AA*powL,aa)/(V*sqrt(powL)),b) ), 1.0/b );
  kn=kl*pow(1.0+powNL,0.333333);

  kll=0.0;  klh=k;                                     /* bracket */
  knl=0.0;  knh=kn;
  
  while( k<=1.0e4*fabs(kn-k) && m < 500 ){

    kl=kll-(knl-k)*(kll-klh)/(knl-knh);
    nin=1+npow(0.75*kl,cosmo)/3.0;
    AA=0.482*pow(nin,-0.947);
    B=0.226*pow(nin,-1.778);
    aa=3.31*pow(nin,-0.244);
    b=0.862*pow(nin,-0.287);
    V=11.55*pow(nin,-0.423);

    powL=5.066e-2*a*a*kl*kl*kl*powerloc(kl,z,cosmo)*g*g/(go*go);
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

/* double powerCDM(double k,double rt){ */

/*   double kn=0.0,kl,powL,powNL,nin; */
/*   double knl,knh,kll,klh,g,a,Omot,Omlt,AA,B,aa,b,V,go; */
/*   double powerloc(double,double),De(double),npow(double); */
/*   int m=0; */
/*   double Omo,Oml; */

/*   if(cosmo->physical){ */
/*     Omo=cosmo->Omo/cosmo->h/cosmo->h; */
/*     Oml=cosmo->Oml/cosmo->h/cosmo->h; */
/*   }else{ */
/*     Omo=cosmo->Omo; */
/*     Oml=cosmo->Oml; */
/*   } */

/*   go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) ); */

/*   a=De(rt); */
/*   if(Omo==1.0){g=go; */
/*   }else{ */
/*     Omot=Omo/(Omo+Oml*a*a*a-a*(Omo+Oml-1)); */
/*     Omlt=a*a*a*Oml*Omot/Omo; */
/*     g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) ); */
/*   } */

/*   kl=k; */
/*   nin=1+npow(0.75*kl)/3.0; */
/*   AA=0.482*pow(nin,-0.947); */
/*   B=0.226*pow(nin,-1.778); */
/*   aa=3.31*pow(nin,-0.244); */
/*   b=0.862*pow(nin,-0.287); */
/*   V=11.55*pow(nin,-0.423); */

/*   powL=5.066e-2*a*a*kl*kl*kl*powerloc(kl,1./a-1.); */
/*   powL*=g*g/(go*go); */
/*   powNL=powL*pow( (1+B*b*powL+pow(AA*powL,aa*b)) */
/* 		    /(1+pow( g*g*g*pow(AA*powL,aa)/(V*sqrt(powL)),b) ), 1.0/b ); */
/*   kn=kl*pow(1.0+powNL,0.333333); */

/*   kll=0.0;  klh=k;                                     /\* bracket *\/ */
/*   knl=0.0;  knh=kn; */
  
/*   while( k<=1.0e4*fabs(kn-k) && m < 500 ){ */

/*     kl=kll-(knl-k)*(kll-klh)/(knl-knh); */
/*     nin=1+npow(0.75*kl)/3.0; */
/*     AA=0.482*pow(nin,-0.947); */
/*     B=0.226*pow(nin,-1.778); */
/*     aa=3.31*pow(nin,-0.244); */
/*     b=0.862*pow(nin,-0.287); */
/*     V=11.55*pow(nin,-0.423); */

/*     powL=5.066e-2*a*a*kl*kl*kl*powerloc(kl,1/a-1.0)*g*g/(go*go); */
/*     powNL=powL*pow( (1+B*b*powL+pow(AA*powL,aa*b)) */
/* 		    /(1+pow( g*g*g*pow(AA*powL,aa)/(V*sqrt(powL)),b) ), 1.0/b ); */
/*     kn=kl*pow(1.0+powNL,0.33333); */

/*     if( (kn-k) < 0.0){ */
/*       kll=kl; */
/*       knl=kn; */
/*     }else{ */
/*       klh=kl; */
/*       knh=kn; */
/*     } */
/*    ++m; */
/*   } */
/*   /\*printf("%i %e %e %e %e\n",m,kl,k,powL,powNL);    /\*test line *\/ */
/*   return 19.739*powNL/(k*k*k*a*a); */
/* } */


/** the scale factor, a, as a function of radius in Mpc **/

double De(double rad,const CosmoHndl cosmo){
  double a[1];
  int nok,nbad;
  void dir(double,double [],double []);
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  if(Omo==1.0) return pow(1-0.5*cosmo->h*rad/3.0e3,2);

  a[0]=1.0;

  if(rad<1.0e1*cosmo->h){ return 1.0 - cosmo->h*rad/3.0e3;
  }else{
	  cosmog = cosmo;
	  odeintD(a-1,1,0.0,rad,1.0e-6,rad/5,rad/1000,&nok,&nbad,dir,bsstepD);
  }
  return a[0];
}

void dir(double r,double a[],double dadr[]){
  double Omo,Oml;

  if(cosmog->physical){
    Omo=cosmog->Omo/cosmog->h/cosmog->h;
    Oml=cosmog->Oml/cosmog->h/cosmog->h;
  }else{
    Omo=cosmog->Omo;
    Oml=cosmog->Oml;
  }

  dadr[1] = -cosmog->h*sqrt( a[1]*(Omo+Oml*pow(a[1],3)+(Omo+Oml-1.0)*a[1]) )/3.0e3;
}

double npow(double k,const CosmoHndl cosmo){
  double qt;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  /* return -1;*/
  /*  qt = k*exp(cosmo->Omb+cosmo->Omb/Omo)/(Omo*cosmo->h*cosmo->h);  */

  if( cosmo->Gamma==0) qt = k*exp(2*cosmo->Omb)/(Omo*cosmo->h*cosmo->h);
  else qt= k/cosmo->Gamma/cosmo->h;

  return cosmo->n-2.0+4.68*qt/(log(1+2.34*qt)*(1+2.34*qt))
    -0.5*qt*(3.89+qt*(5.1842e2+qt*(4.8831e2+8.1088e3*qt)))/( 1+qt*(3.89+qt*(2.5921e2+qt*(1.6277e2+2.0272e3*qt))) );
}

/* linear power spectrum P(k,z)/a^2 */
double power_linear(double k,double z,const CosmoHndl cosmo){
  double powerloc(double k,double z,const CosmoHndl cosmo);

  if(z==0.0) return powerloc(k,0,cosmo);

/*   go=2.5*Omo/( pow(Omo,4.0/7.0)-cosmo->Oml+(1+0.5*Omo)*(1+cosmo->Oml/70) ); */
/*   aa=1/(1+z); */

/*   if(Omo==1.0){g=go; */
/*   }else{ */
/*     Omot=Omo/(Omo+cosmo->Oml*aa*aa*aa-aa*(Omo+cosmo->Oml-1)); */
/*     Omlt=aa*aa*aa*cosmo->Oml*Omot/Omo; */
/*     g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) ); */
/*   } */

  return pow(Dgrowth(z,cosmo)*(1+z),2)*powerloc(k,z,cosmo);
  /*return g*g*powerloc(k,z)/go/go;*/
}

/** linear power spectrum without growth factor **/
/**   growth factor should be normalized to 1 at z=0 **/
double powerloc(double k,double z,const CosmoHndl cosmo){
  double qt,ans;
  //double powerEH(double,double);
  //double powerEHv2(double);

//  return powerEHv2(k);
  return powerEH(k,z,cosmo);
  if(cosmo->Omnu>0.0) return powerEH(k,z,cosmo);  /* with neutrinos but no BAO */
  else return powerEHv2(k,cosmo);/* with BAO but no neutrinos */

  if( cosmo->Gamma==0){
    if(cosmo->physical) qt = k*exp(2*cosmo->Omb/cosmo->h/cosmo->h)/(cosmo->Omo);
    else qt = k*exp(2*cosmo->Omb)/(cosmo->Omo*cosmo->h*cosmo->h);
  }else{
    qt= k/cosmo->Gamma/cosmo->h;
  }

  ans=cosmo->A*pow(k/cosmo->h,cosmo->n)*pow( log(1+2.34*qt)/(2.34*qt) ,2)
    /sqrt(1+qt*(3.89+qt*(2.5921e2+qt*(1.6277e2+2.0272e3*qt))) );

  ans/=pow(cosmo->h/3.0e3,3);

  return ans;
}

double power_normalize(double sigma8,CosmoHndl cosmo){
  double normL(double lgk),powfactor;

  cosmo->sig8=sigma8;
  cosmo->A=1.0;
  cosmog = cosmo;
  powfactor=9*nintegrateDcos(normL,log(1.0e-3),log(1.0e4),1.0e-9)/(2*pi*pi);
  cosmo->A=sigma8*sigma8/powfactor;                        /* linear normalization */

  return powfactor;
}

double normL(double lgk){
  double R,win,k;

  k=exp(lgk);
  R=k*8/cosmog->h;
  if(R<=1.0e-4){ win = (1-R*R/10)/3; 
  }else{win=(sin(R)/R - cos(R))/(R*R);}

  return k*k*k*powerloc(k,0,cosmog)*win*win;
}

/** this is the power spectrum from Eisinstein & Hu **/
/** with neutrinos but no BAO  **/
double powerEH(double k,double z,const CosmoHndl cosmo){
  static struct cosmology cosmo_old;
  static double zloc=-100;
  double Trans;

  /*PrintCosmology(cosmo_old);*/
  /*  PrintCosmology(cosmo_old);*/
  /*printf("compare = %i\n",cosmo_compare(cosmo_old,cosmo));*/

  /** remove this after tests **/
  /** if( zin < 10) z=7;  **/

  if( z != zloc){
    if(cosmo->physical) {
      TFmdm_set_cosm((cosmo->Omo/cosmo->h/cosmo->h)
		     ,(cosmo->Omb/cosmo->h/cosmo->h)
		     ,(cosmo->Omnu/cosmo->h/cosmo->h)
		     ,cosmo->Nnu,(cosmo->Oml/cosmo->h/cosmo->h)
		     ,(cosmo->h),(z));
    }else{
      //PrintCosmology(cosmo);
      //printf("%e %e %e %e %e %e %e\n",(cosmo->Omo),(cosmo->Omb),(cosmo->Omnu)
    	//,(cosmo->Nnu),(cosmo->Oml),(cosmo->h),(z));
      TFmdm_set_cosm((cosmo->Omo),(cosmo->Omb),(cosmo->Omnu)
		     ,cosmo->Nnu,(cosmo->Oml),(cosmo->h),(z));
 
    }
    cosmo_copy(&cosmo_old,cosmo);
    zloc=z;
  }else if(cosmo_compare(&cosmo_old,cosmo)){
    if(cosmo->physical) {
      TFmdm_set_cosm((cosmo->Omo/cosmo->h/cosmo->h)
		     ,(cosmo->Omb/cosmo->h/cosmo->h)
		     ,(cosmo->Omnu/cosmo->h/cosmo->h)
		     ,cosmo->Nnu,(cosmo->Oml/cosmo->h/cosmo->h)
		     ,(cosmo->h),(z));
    }else{
    TFmdm_set_cosm((cosmo->Omo),(cosmo->Omb),(cosmo->Omnu)
		   ,cosmo->Nnu,(cosmo->Oml),(cosmo->h),(z));
     }

    cosmo_copy(&cosmo_old,cosmo);
  }
 
  Trans=TFmdm_onek_mpc(k);
  //printf("trans=%e A=%e h=%e n=%e\n",Trans,cosmo->A,cosmo->h,cosmo->n);
  return cosmo->A*pow(k/cosmo->h,cosmo->n+cosmo->dndlnk*log(k))*Trans*Trans/pow(cosmo->h/3.0e3,3);
}

/** this is the power spectrum from Eisinstein & Hu **/
/** with BAO bcosmo->A*ut no neutrinos  **/
double powerEHv2(double k,const CosmoHndl cosmo){
  static struct cosmology cosmo_old;
  double Trans;
  double baryon_piece,cdm_piece;

  /*PrintCosmology(cosmo_old);*/
  /*  PrintCosmology(cosmo_old);*/
  /*printf("compare = %i\n",cosmo_compare(cosmo_old,cosmo));*/
  
  if(cosmo_compare(&cosmo_old,cosmo)){
    if(cosmo->physical) {
      TFset_parameters((cosmo->Omo),(cosmo->Omb/cosmo->Omo),2.728);
    }else{
      TFset_parameters((cosmo->Omo*cosmo->h*cosmo->h),(cosmo->Omb/cosmo->Omo),2.728);
    }
    cosmo_copy(&cosmo_old,cosmo);
  }
 
  Trans=TFfit_onek(k, &baryon_piece, &cdm_piece);
  //printf("trans=%e A=%e h=%e n=%e\n",Trans,cosmo->A,cosmo->h,cosmo->n);
  return cosmo->A*pow(k/cosmo->h,cosmo->n)*Trans*Trans/pow(cosmo->h/3.0e3,3);
}


