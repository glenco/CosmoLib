/********************************************
  powerCDM.c calculates the nonlinear P(k,r)/a(r)^2
    example of how to use is in powerCDMt.c
  *******************************************/
#include <math.h>
//#include "Recipes/nrutil.h"
//#include "RecipesD/nrD.h"
#include <nrutil.h>
#include <nrD.h>
#include <cosmo.h>

//extern struct cosmology cosmo;
static CosmoHndl cosmog;

double powerCDMdark(double k,double z,const CosmoHndl cosmo){

  double kn=0.0,kl,powL,powNL,nin;
  double knl,knh,kll,klh,a,AA,B,aa,b,V;
  //double powerloc(double,double),npow(double);
  //double growth_dark(double a),powerCDMz(double,double);
  int m=0;

  double g,go;

  if(cosmo->w == -1.0 && cosmo->w1== 0.0) return powerCDMz(k,z,cosmo);

  a=1.0/(1+z);

  if(cosmo->growth == NULL){
      go=growth_dark(1.0,cosmo);
      g=growth_dark(a,cosmo);
  }else{
      go=growth_table(1.0,cosmo->a,cosmo->growth,cosmo->Ntable);
      g=growth_table(a,cosmo->a,cosmo->growth,cosmo->Ntable);
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
  //printf("%i %e %e %e %e\n",m,kl,k,powL,powNL);
  return 19.739*powNL/(k*k*k*a*a);
}


double dEx_darkda(double x,const CosmoHndl cosmo){
  double ao=1.0;

  if(cosmo->physical) return  -(1.5*cosmo->Omo*pow(x,4)
	  +1.5*cosmo->Oml*((1+cosmo->w+ao*cosmo->w1)*x-cosmo->w1)
            *pow(x,3*(1+cosmo->w+ao*cosmo->w1))*exp(-3*cosmo->w1*(x-1)/x)
	   +(cosmo->h*cosmo->h-cosmo->Omo-cosmo->Oml)*pow(x,3))/Ex_dark(x,cosmo)/cosmo->h/cosmo->h;

  return -(1.5*cosmo->Omo*pow(x,4)
	  +1.5*cosmo->Oml*((1+cosmo->w+ao*cosmo->w1)*x-cosmo->w1)
            *pow(x,3*(1+cosmo->w+ao*cosmo->w1))*exp(-3*cosmo->w1*(x-1)/x)
	   +(1-cosmo->Omo-cosmo->Oml)*pow(x,3))/Ex_dark(x,cosmo);
}

/******* the growth factor relative to matter dominated *******/
double growth_dark1(double a,const CosmoHndl cosmo){
  double Omot,Omlt,alpha,aa;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  if(fabs(Omo-1.0)<1.0e-5) return 1.0;

  Omot=Omo/pow(Ex_dark(1/a,cosmo),2)/pow(a,3);
  Omlt=Oml*pow(a,-3*(1+cosmo->w))/pow(Ex_dark(1/a,cosmo),2);
  alpha= 3/(5-cosmo->w/(1-cosmo->w))
    + (1-Omot)*3*(1-cosmo->w)*(1-3*cosmo->w/2)/pow(1-6*cosmo->w/5,3)/125;
  aa=-0.28/(cosmo->w+0.08) -0.3;

  return 5*Omot/( pow(Omot,alpha) - Omlt + (1+0.5*Omot)*(1+aa*Omlt) )/2;
}

/** fitting formula of E. Linder 2005 **/
double growth_dark2(double a,const CosmoHndl cosmo){
  double dg2dla(double lna);

  if(a == 1.0) return 1.0;
  if(1.0-a < 1.0e-6){
    if(cosmo->physical)  return 1+( pow(cosmo->Omo/cosmo->h/cosmo->h,0.55+0.05*(1+cosmo->w+0.5*cosmo->w1)) -1)*(a-1)/a;
    return 1+( pow(cosmo->Omo,0.55+0.05*(1+cosmo->w+0.5*cosmo->w1)) -1)*(a-1)/a;
  }

  cosmog = cosmo;
  return exp(nintegrateDcos(dg2dla,0,log(a),1.0e-6));
}

double dg2dla(double lna){
  double a,Omot;

  a=exp(lna);
  if(cosmog->physical) Omot=cosmog->Omo/pow(Ex_dark(1/a,cosmog),2)/pow(a,3)/cosmog->h/cosmog->h;
  else Omot=cosmog->Omo/pow(Ex_dark(1/a,cosmog),2)/pow(a,3);

  if(cosmog->darkenergy == 2) return pow(Omot,cosmog->gamma)-1;
  return pow(Omot,0.55+0.05*(1+cosmog->w+0.5*cosmog->w1))-1;
}

double growth_dark(double a,const CosmoHndl cosmo){
  double y[2];
  int nok,nbad;
  void dgrowth_dark(double lna,double y[],double dydlna[]);
  double Omot,Omlt,go;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) );

  if( a==1.0) return go;
  if(fabs(Omo-1.0)<1.0e-5) return go;

  /* cosmological constant case */
  if( (cosmo->w==-1) && (cosmo->w1==0) ){
    Omot=Omo/(Omo+Oml*a*a*a-a*(Omo+Oml-1));
    Omlt=a*a*a*Oml*Omot/Omo;
    return 2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
  }

  return go*growth_dark2(a,cosmo);

  /*  if(cosmo->w1==0) return growth(a); */

  y[0]=1.0;
  y[1]=0.0;

  cosmog = cosmo;
  odeintD(y-1,2,log(1.0e-3),log(a),1.0e-6,log(1.0e3/a)/5.,log(1.0e3/a)/1000.,&nok,&nbad,dgrowth_dark,bsstepD);
  return go*y[0]/a;
}


/* linearly interpolates the growth function from a table */
double growth_table(double a,double *a_table,double *g_table,int Ntable){
    unsigned long k,j,m=3;
    double g,dg;

    locateD(a_table-1,Ntable,a,&j);
    k=IMIN(IMAX(j-(m-1)/2,1),Ntable+1-m)-1;
    polintD(&a_table[k-1],&g_table[k-1],m,a,&g,&dg);
    return g;
}
/* makes lookup table for growth function */
void make_growth_table(double z,double *a_table,double *g_table,int Ntable,const CosmoHndl cosmo){

  double y[2];
  int nok,nbad,i;
  void dgrowth_dark(double lna,double y[],double dydlna[]);
  double Omot,Omlt,go;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) );

  for(i=0;i<Ntable;++i){

      a_table[i]=1-i*(1-1./(1+z))/(Ntable-1);
/*      a_table[i]=1./(1+z)+i*(1-1./(1+z))/(Ntable-1);*/

/*       if( a_table[i]==1.0) g_table[i]= go; */
/*       else if(fabs(Omo-1.0)<1.0e-5) g_table[i]= go; */
/*       /\* cosmological constant case *\/ */
/*       else  */
      if( (cosmo->w==-1) * (cosmo->w1==0) ){
	  Omot=Omo/(Omo+Oml*a_table[i]*a_table[i]*a_table[i]-a_table[i]*(Omo+Oml-1));
	  Omlt=a_table[i]*a_table[i]*a_table[i]*Oml*Omot/Omo;
	  g_table[i]=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
      }else{
	  /* calculate by brute force */
	  /* this will be normalized to 1 at a=1.0e-3 */
	  y[0]=1.0;
	  y[1]=0.0;

	  cosmog = cosmo;
	  odeintD(y-1,2,log(1.0e-3),log(a_table[i]),1.0e-6,log(1.0e3/a_table[i])/5.,log(1.0e3/a_table[i])/1000.,&nok,&nbad,dgrowth_dark,bsstepD);
	  g_table[i]=y[0]/a_table[i];
      }
  }

  /* renormalize to z=0 relative to Omega=1.0 */
  for(i=Ntable-1;i>=0;--i) g_table[i]=go*g_table[i]/g_table[0];
}

void dgrowth_dark(double lna,double y[],double dydlna[]){
  double a,Et,dEdat;

  a=exp(lna);
  Et=Ex_dark(1/a,cosmog);
  dEdat=dEx_darkda(1/a,cosmog);

  dydlna[1]=y[2];
  if(cosmog->physical) dydlna[2]=-(2+a*dEdat/Et)*y[2]
    +1.5*cosmog->Omo*pow(a,-3)/Et/Et*y[1]/cosmog->h/cosmog->h;
  else dydlna[2]=-(2+a*dEdat/Et)*y[2]
    +1.5*cosmog->Omo*pow(a,-3)/Et/Et*y[1];
}

/* the supression factor to the linear power spectrum caused by massive */
/* neutrinos according to Kiakotou, Elgaroy & Lahav, astro-ph/0709.0253v2 */

/*  this has not been tested !!!!!*/

double growthKEL(double z,double k,const CosmoHndl cosmo){
  double dlnddlnaKEL(double lna),muKEL(double,const CosmoHndl cosmo),go;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) );

  cosmog = cosmo;
  return go*exp(-2*muKEL(k,cosmo)*nintegrateDcos(dlnddlnaKEL,0,log(1.0/(1+z)),1.0e-8));
}

double dlnddlnaKEL(double lna){
  float a,Omot,alpha;
 
  a=exp(lna);
  if(cosmog->physical) Omot = cosmog->Omo
		  /(cosmog->Omo + cosmog->Oml*a*a*a - a*(cosmog->Omo+cosmog->Oml-cosmog->h*cosmog->h));
  else Omot=cosmog->Omo/(cosmog->Omo+cosmog->Oml*a*a*a-a*(cosmog->Omo+cosmog->Oml-1));
  alpha=3./(5.-cosmog->w/(1-cosmog->w)) + (1-Omot)*3*(1-cosmog->w)*(1-3*cosmog->w/2)*pow(1-6*cosmog->w/5,-3)/125;
  return pow(Omot,alpha);
}

double muKEL(double k,const CosmoHndl cosmo){
  double fnu;
  double At[6],Bt[6],Ct[6],kh[6],A,B,C;
  unsigned long j;

  if( k <=  0.001*cosmo->h) return 1.0;
  if( k > 0.5*cosmo->h) return pow(1-cosmo->Omnu/cosmo->Omo,3/(5-cosmo->w/(1-cosmo->w)) );

  kh[0]=0.001; kh[1]=0.01; kh[2]=0.05; kh[3]=0.07; kh[4]=0.1; kh[5]=0.5;
  At[0] =0;     At[1]=0.132;  At[2]=0.613; At[3]=0.733; At[4]=0.786; At[5]=0.813;
  Bt[0] =0;     Bt[1]=1.62;  Bt[2]=5.59; Bt[3]=6.0; Bt[4]=5.09; Bt[5]=0.803;
  Ct[0] =0;     Ct[1]=7.13;  Ct[2]=21.13; Ct[3]=21.45; Ct[4]=15.5; Ct[5]=-0.844;

  locateD(kh-1,6,k/cosmo->h,&j);
  A=At[j] + (k/cosmo->h - kh[j])*(At[j+1]-At[j])/(kh[j+1]-kh[j]);
  B=Bt[j] + (k/cosmo->h - kh[j])*(Bt[j+1]-Bt[j])/(kh[j+1]-kh[j]);
  C=Ct[j] + (k/cosmo->h - kh[j])*(Ct[j+1]-Ct[j])/(kh[j+1]-kh[j]);

  fnu=cosmo->Omnu/cosmo->Omo;
  return 1 - A*cosmo->Oml*fnu + B*fnu*fnu - C*fnu*fnu*fnu;
}
