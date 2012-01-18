/** ******************************************
  powerCDM.c calculates the nonlinear P(k,z)/a(r)^2
    example of how to use is in testpower.c
  *******************************************/
#include <math.h>
#include <nrD.h>
#include <cosmo.h>

double COSMOLOGY::powerCDMz(double k,double z){

  double kn=0.0,kl,powL,powNL,nin;
  double knl,knh,kll,klh,g,a,Omot,Omlt,AA,B,aa,b,V,go;
  int m=0;
  double omo,oml;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }

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

/* double powerCDM(double k,double rt){ */

/*   double kn=0.0,kl,powL,powNL,nin; */
/*   double knl,knh,kll,klh,g,a,Omot,Omlt,AA,B,aa,b,V,go; */
/*   double powerloc(double,double),De(double),npow(double); */
/*   int m=0; */
/*   double Omo,Oml; */

/*   if(physical){ */
/*     Omo=Omo/h/h; */
/*     Oml=Oml/h/h; */
/*   }else{ */
/*     Omo=Omo; */
/*     Oml=Oml; */
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

double COSMOLOGY::De(double rad){
  double a[1];
  int nok,nbad;
  void dir(double,double [],double []);
  double omo,oml;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }

  if(omo==1.0) return pow(1-0.5*h*rad/3.0e3,2);

  a[0]=1.0;

  if(rad<1.0e1*h){ return 1.0 - h*rad/3.0e3;
  }else{
	  odeintD(a-1,1,0.0,rad,1.0e-6,rad/5,rad/1000,&nok,&nbad,dir,bsstepD);
  }
  return a[0];
}

void COSMOLOGY::dir(double r,double a[],double dadr[]){
  double omo,oml;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }

  dadr[1] = -h*sqrt( a[1]*(omo+oml*pow(a[1],3)+(omo+oml-1.0)*a[1]) )/3.0e3;
}

double COSMOLOGY::npow(double k){
  double qt, omo, oml;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }

  /* return -1;*/
  /*  qt = k*exp(Omb+Omb/Omo)/(Omo*h*h);  */

  if( Gamma==0) qt = k*exp(2*Omb)/(omo*h*h);
  else qt= k/Gamma/h;

  return n-2.0+4.68*qt/(log(1+2.34*qt)*(1+2.34*qt))
    -0.5*qt*(3.89+qt*(5.1842e2+qt*(4.8831e2+8.1088e3*qt)))/( 1+qt*(3.89+qt*(2.5921e2+qt*(1.6277e2+2.0272e3*qt))) );
}

/** linear power spectrum P(k,z)/a^2 */
double COSMOLOGY::power_linear(double k,double z){

  if(z==0.0) return powerloc(k,0);

/*   go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) ); */
/*   aa=1/(1+z); */

/*   if(Omo==1.0){g=go; */
/*   }else{ */
/*     Omot=Omo/(Omo+Oml*aa*aa*aa-aa*(Omo+Oml-1)); */
/*     Omlt=aa*aa*aa*Oml*Omot/Omo; */
/*     g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) ); */
/*   } */

  return pow(Dgrowth(z)*(1+z),2)*powerloc(k,z);
  /*return g*g*powerloc(k,z)/go/go;*/
}

/** linear power spectrum without growth factor
**   growth factor should be normalized to 1 at z=0
**   **/
double COSMOLOGY::powerloc(double k,double z){
  //double qt,ans;
  //double powerEH(double,double);
  //double powerEHv2(double);

//  return powerEHv2(k);
  return powerEH(k,z);
}

double COSMOLOGY::power_normalize(double sigma8){
  double powfactor;

  sig8=sigma8;
  A=1.0;

  powfactor=9*nintegrateDcos(&COSMOLOGY::normL,log(1.0e-3),log(1.0e4),1.0e-9)/(2*pi*pi);
  A=sigma8*sigma8/powfactor;                        /* linear normalization */

  return powfactor;
}

double COSMOLOGY::normL(double lgk){
  double R,win,k;

  k=exp(lgk);
  R=k*8/h;
  if(R<=1.0e-4){ win = (1-R*R/10)/3; 
  }else{win=(sin(R)/R - cos(R))/(R*R);}

  return k*k*k*powerEH(k,0)*win*win;
}

/** this is the power spectrum from Eisinstein & Hu
 * with neutrinos but no BAO
 * **/
double COSMOLOGY::powerEH(double k,double z){
  COSMOLOGY cosmo_old;
  static double zloc=-100;
  double Trans;

  /*PrintCosmology(cosmo_old);*/
  /*  PrintCosmology(cosmo_old);*/
  /*printf("compare = %i\n",cosmo_compare(cosmo_old,cosmo));*/

  /** remove this after tests **/
  /** if( zin < 10) z=7;  **/

  if( z != zloc){
    if(physical) {
      TFmdm_set_cosm((Omo/h/h)
		     ,(Omb/h/h)
		     ,(Omnu/h/h)
		     ,Nnu,(Oml/h/h)
		     ,(h),(z));
    }else{
      //PrintCosmology(cosmo);
      //printf("%e %e %e %e %e %e %e\n",(Omo),(Omb),(Omnu)
    	//,(Nnu),(Oml),(h),(z));
      TFmdm_set_cosm((Omo),(Omb),(Omnu)
		     ,Nnu,(Oml),(h),(z));
 
    }
    cosmo_copy(&cosmo_old,this);
    zloc=z;
  }else if(cosmo_compare(&cosmo_old,this)){
    if(physical) {
      TFmdm_set_cosm((Omo/h/h)
		     ,(Omb/h/h)
		     ,(Omnu/h/h)
		     ,Nnu,(Oml/h/h)
		     ,(h),(z));
    }else{
    TFmdm_set_cosm((Omo),(Omb),(Omnu)
		   ,Nnu,(Oml),(h),(z));
     }

    cosmo_copy(&cosmo_old,this);
  }
 
  Trans=TFmdm_onek_mpc(k);
  //printf("trans=%e A=%e h=%e n=%e\n",Trans,A,h,n);
  return A*pow(k/h,n+dndlnk*log(k))*Trans*Trans/pow(h/3.0e3,3);
}

/** this is the power spectrum from Eisinstein & Hu **/
/** with BAO bA*ut no neutrinos  **/
double COSMOLOGY::powerEHv2(double k){
  COSMOLOGY cosmo_old;
  double Trans;
  double baryon_piece,cdm_piece;

  /*PrintCosmology(cosmo_old);*/
  /*  PrintCosmology(cosmo_old);*/
  /*printf("compare = %i\n",cosmo_compare(cosmo_old,cosmo));*/
  
  if(cosmo_compare(&cosmo_old,this)){
    if(physical) {
      TFset_parameters((Omo),(Omb/Omo),2.728);
    }else{
      TFset_parameters((Omo*h*h),(Omb/Omo),2.728);
    }
    cosmo_copy(&cosmo_old,this);
  }
 
  Trans=TFfit_onek(k, &baryon_piece, &cdm_piece);
  //printf("trans=%e A=%e h=%e n=%e\n",Trans,A,h,n);
  return A*pow(k/h,n)*Trans*Trans/pow(h/3.0e3,3);
}


