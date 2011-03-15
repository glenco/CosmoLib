/********************************************

  cosmo.c calculates some useful comological parameters

to run:

struct cosmology{
  double h;
  double n;
  double A;
  double Omo;
  double Oml;
  double Omb;
  double w;
  double w1;
  double Gamma;
};

#include "../Library/cosmo.c"
#include "../Library/Recipes/polint.c"


int kmax,kount;
double *xp,**yp,dxsav;

in main:
  xp=dvector(1,1);                  
  yp=dmatrix(1,1,1,1);
  kmax=0.0;
  dxsav=1.0;

*******************************************
*******************************************/
#define CRITD 2.49783e18 /* critical density / Ho^2 ; solar/Mpc */ 

#include <math.h>
//#include <nrD.h>
//#include "../Library/Recipes/nrutil.h"
#include <nrD.h>
#include <nrutil.h>
#include <cosmo.h>

/*#include "../Library/bsstepD.c"*/
/*#include "../Library/odeintD.c"*/
/*#include "../Library/mmidD.c"*/
/*#include "../Library/pzextrD.c"*/
/*#include "../Library/dfridrD.c"*/
/** things for integrator **/
#define JMAX 34
#define JMAXP (JMAX+1)
#define K 6
#define NRANSI

#define FUNC(x) ((*func)(x))

int kmax,kount;
double *xp,**yp,dxsav;

//extern struct cosmology cosmo;
static CosmoHndl cosmog;
static double alph;  /* DR-distance parameter */

void SetConcordenceCosmology(CosmoHndl cosmo){
	// set cosmological parameters to standard WMAP 5r values
	// Komatsu et al. 2009, ApJ 180, 330
	// does not set power spectrum normalization
	// which needs to be done separately

	/* default parameterization */
	if( (cosmo->physical != 0)*(cosmo->physical != 1) ) cosmo->physical = 0;

	cosmo->Omo=0.1358;
	cosmo->Omb=0.02267;
	cosmo->h=0.705;

	if(cosmo->physical==0){ /* use standard parameter set */
		cosmo->Omo/=cosmo->h*cosmo->h;
		cosmo->Omb/=cosmo->h*cosmo->h;
		cosmo->Oml=1.0-cosmo->Omo;
	}else if(cosmo->physical==1){
		cosmo->Oml=(1-cosmo->Omo/cosmo->h/cosmo->h)*cosmo->h*cosmo->h;
	}

	cosmo->w=-1.0;
	cosmo->w1=0.0;
	cosmo->n=1.0;
	cosmo->Gamma=0.0;
	cosmo->Omnu=0;
	cosmo->Nnu=3.0;
	cosmo->dndlnk=0.0;
	cosmo->gamma=0.55;

	cosmo->darkenergy=1;
	/* if 2 gamma parameterization is used for dark energy */
	/* if 1 w,w_1 parameterization is used for dark energy */
	power_normalize(0.812,cosmo);
}

void PrintCosmology(COSMOLOGY cos1){
	  printf("h: %e n: %e dndlnk: %e\nA: %e sig8: %e Gamma: %e \n"
			  ,cos1.h,cos1.n,cos1.dndlnk,cos1.A,cos1.sig8,cos1.Gamma);
	if(cos1.physical==0){
	  printf("Omo: %e Oml: %e Omb: %e\nOmnu : %e Nnu : %e\n"
			  ,cos1.Omo,cos1.Oml,cos1.Omb,cos1.Omnu,cos1.Nnu);
  }
  if(cos1.physical==1){
	  printf("Omo hh: %e Oml hh: %e Omb hh: %e\nOmnu hh: %e Nnu : %e\n"
				  ,cos1.Omo,cos1.Oml,cos1.Omb,cos1.Omnu,cos1.Nnu);
  }
  if(cos1.darkenergy==2) printf("darkenery=%i gamma=%e\n",cos1.darkenergy,cos1.gamma);
  else printf("darkenery: %i w: %f w1: %f\n",cos1.darkenergy,cos1.w,cos1.w1);
}

/** see if cosmologies are identical **/
int cosmo_compare(COSMOLOGY *cos1,COSMOLOGY *cos2){

  return 1-(cos1->h == cos2->h)*(cos1->n == cos2->n)*(cos1->A == cos2->A)*(cos1->Omo == cos2->Omo)*(cos1->Oml == cos2->Oml)*(cos1->Omb == cos2->Omb)*(cos1->Omnu == cos2->Omnu)*(cos1->w == cos2->w)*(cos1->w1 == cos2->w1)*(cos1->Gamma == cos2->Gamma)*(cos1->Nnu == cos2->Nnu);
}
void cosmo_copy(struct cosmology *cos1,struct cosmology *cos2){
	cos1->physical=cos2->physical;
	cos1->Omo=cos2->Omo;
	cos1->Oml=cos2->Oml;
	cos1->Omb=cos2->Omb;
	cos1->h=cos2->h;
	cos1->w=cos2->w;
	cos1->w1=cos2->w1;
	cos1->n=cos2->n;
	cos1->Gamma=cos2->Gamma;
	cos1->Omnu=cos2->Omnu;
	cos1->Nnu=cos2->Nnu;
	cos1->dndlnk=cos2->dndlnk;
	cos1->gamma=cos2->gamma;
	cos1->darkenergy=cos2->darkenergy;
	cos1->A=cos2->A;
	cos1->sig8=cos2->sig8;
}

/*** The curvature radius in Mpc ***/

double rcurve(const CosmoHndl cosmo){
  if(cosmo->Omo+cosmo->Oml != 1.0){
    if(cosmo->physical) return 3.0e3/(sqrt(fabs(cosmo->h*cosmo->h-cosmo->Omo-cosmo->Oml)));  /** curviture scale **/
    return 3.0e3/(cosmo->h*sqrt(fabs(1-cosmo->Omo-cosmo->Oml)));  /** curviture scale **/
  }
  return 0;
}

/****** derivative of the coordinate distance in units of Ho^-1, x=1+z ****/

double radiusm(double x){
 double temp;
  /*printf("cosmog->->Omo=%f cosmog->->Oml=%e x=%e\n",cosmog->->Omo,cosmog->->Oml,x);*/
  /*  if( cosmog->->Omo==1.0) return 2.0*(1-1/sqrt(x)); */

  if(cosmog->physical) temp=(cosmog->Omo*x*x*x+cosmog->Oml-(cosmog->Omo+cosmog->Oml-cosmog->h*cosmog->h)*x*x)/cosmog->h/cosmog->h;
  temp=cosmog->Omo*x*x*x+cosmog->Oml-(cosmog->Omo+cosmog->Oml-1)*x*x;
  if(temp<=0.0) return -1.0e30;                   // nonphysical values
  return 1.0/sqrt(temp);
}

double radiusm_dark(double x){
  double Ex_dark(double,const CosmoHndl cosmo),Ex(double,const CosmoHndl cosmo);

  if(cosmog->w == -1.0 && cosmog->w1 == 0.0) return 1.0/Ex(x,cosmog);
  return 1.0/Ex_dark(x,cosmog);
}

double Ex_dark(double x,const CosmoHndl cosmo){
  double ao=1.0;
  if(cosmo->physical) return sqrt(cosmo->Omo*x*x*x+cosmo->Oml*pow(x,3*(1+cosmo->w+ao*cosmo->w1))*exp(-3*cosmo->w1*(x-1)/x)-(cosmo->Omo+cosmo->Oml-cosmo->h*cosmo->h)*x*x)/cosmo->h;

  return sqrt(cosmo->Omo*x*x*x+cosmo->Oml*pow(x,3*(1+cosmo->w+ao*cosmo->w1))*exp(-3*cosmo->w1*(x-1)/x)-(cosmo->Omo+cosmo->Oml-1)*x*x);
}

/***************************************************************/
/***** comoving Dyer-Roeder angular size distance for lambda=0 ****/
/***************************************************************/

double DRradius(double zo,double z,double pfrac,CosmoHndl cosmo){
  double zl,Omr,Omz,Da[2],Ho;
  void ders(double,double *,double *);
  int nok,nbad;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  Ho=cosmo->h/3.0e3;
  alph=1.5*(1-pfrac);

    /*    if(Oml !=0.0){
        printf("ERROR in DRradius: Oml != 0\n");
        return 0.0;
        }*/
    zl=1+z;

    if(zo==0 && (Oml==0 && alph==0)){
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
      cosmog = cosmo;
      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
    return 0;
}

double DRradius2(double zo,double z,CosmoHndl cosmo){
  double zl,Omr,Omz,Da[2],ds,dl,Ho;
    void ders(double,double *,double *);
    int nok,nbad;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

    Ho=cosmo->h/3.0e3;
    /*    if(Oml !=0.0){
        printf("ERROR in DRradius: Oml != 0\n");
        return 0.0;
        }*/
    alph=0;
    zl=1+z;

    if(Omo==1.0){
      ds=2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
      if(zo==0.0){ return ds;
      }else{
	zl=zo+1;
	dl=2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
 	return   ds-dl*(1+z)/(1+zo);/**/
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
      cosmog = cosmo;
      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
}

void ders(double z,double Da[],double dDdz[]){
  double Omo,Oml;

  if(cosmog->physical){
    Omo=cosmog->Omo/cosmog->h/cosmog->h;
    Oml=cosmog->Oml/cosmog->h/cosmog->h;
  }else{
    Omo=cosmog->Omo;
    Oml=cosmog->Oml;
  }

  dDdz[1]=Da[2];
  dDdz[2]=-3*Da[2]/(1+z) - (0.5*(Omo-2*Oml*pow(1+z,-3))*Da[2]+alph*Omo*Da[1]/(1+z) )/(1+z*Omo+(pow(1+z,-2)-1)*Oml);
}

double arctanh(double x){
    return 0.5*log((1+x)/(1-x));
}

double fmini(double a,double b){
  if(a<b)return a;
  return b;
}
double fmaxi(double a,double b){
  if(a>b)return a;
  return b;
}

/** linear growth factor normalized to 1 at z=0 **/

double Dgrowth(double z,const CosmoHndl cosmo){
  double a,Omot,Omlt,g,go;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  a=1./(1+z);
  if(Omo==1 && Oml==0){ return a;
  }else{
    go=2.5*Omo/( pow(Omo,4.0/7.0)-Oml+(1+0.5*Omo)*(1+Oml/70) );
    Omot=Omo/(Omo+Oml*a*a*a-a*(Omo+Oml-1));
    Omlt=a*a*a*Oml*Omot/Omo;
    g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
    return a*g/go;
  }
}

/****** the coordinate distance in units Mpc  ****/

double coorDist(double zo,double z,const CosmoHndl cosmo){
  //double radiusm_dark(double x);

  //printf("In coorDist\n");
  cosmog = cosmo;
  //PrintCosmology(*cosmog);
  if( (cosmo->w ==-1.0) && (cosmo->w1 == 0.0) ) return nintegrateDcos(radiusm,1+zo,1+z,1.0e-9)*3.0e3/cosmo->h;
  return nintegrateDcos(radiusm_dark,1+zo,1+z,1.0e-9)*3.0e3/cosmo->h;
}

double angDist(double zo,double z,const CosmoHndl cosmo){
  double Omo,Oml,Rcur;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  if(Omo+Oml==1) return coorDist(zo,z,cosmo)/(1+z);
   /** curviture scale **/

  Rcur=3.0e3/(cosmo->h*sqrt(fabs(1-Omo-Oml)));
  if((Omo+Oml)<1.0) return Rcur*sinh(coorDist(zo,z,cosmo)/Rcur)/(1+z);
  return Rcur*sin(coorDist(zo,z,cosmo)/Rcur)/(1+z);
}

double lumDist(double zo,double z,const CosmoHndl cosmo){
	return pow(1+z,2)*angDist(zo,z,cosmo);
}

/****** angular size distence - R is the curviture radius, rd the coordinate ***/

double gradius(double R,double rd,const CosmoHndl cosmo){
  /*  printf("cosmo->Omo=%e cosmo->Oml=%e (cosmo->Omo+cosmo->Oml)-1=*/
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  if( fabs(Omo+Oml-1) < 1.0e-4) return rd;
  if((Omo+Oml)<1.0) return R*sinh(rd/R);
  return R*sin(rd/R);
}


/***************************************************************/
/** Press-Schechter mass function in fraction of mean density **/
/***************************************************************/

double psdfdm(double sig8,double z,double m,const CosmoHndl cosmo){
  double dc,Omz,sig,Dg;
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  Dg=Dgrowth(z,cosmo);
  sig=sig8*Deltao(m,cosmo);
  Omz=Omo*pow(1+z,3)/( Omo*pow(1+z,3)+(1-Omo-Oml)*pow(1+z,2)+Oml );

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  /*printf("dc=%e dsigdM=%e sig=%e m=%e D=%e\n",dc,dsigdM(m),sig,m,Dg);*/
  return -0.797885*dc*sig8*dsigdM(m,cosmo)*exp(-0.5*pow(dc/(Dg*sig),2) )
    /(Dg*sig*sig);
}

/***************************************************************/
/** Sheth-Tormen, most be normalized to 1 by calling code **/
/***************************************************************/

double stdfdm(double sig8,double z,double m,const CosmoHndl cosmo){
  double dc,Omz,sig,Dg
    ,Dgrowth(double,const CosmoHndl),Deltao(double,const CosmoHndl)
    ,psdfdm(double,double,double,const CosmoHndl cosmo),dsigdM(double,const CosmoHndl cosmo);
  double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  Dg=Dgrowth(z,cosmo);
  sig=sig8*Deltao(m,cosmo);
  Omz=Omo*pow(1+z,3)/( Omo*pow(1+z,3)+(1-Omo-Oml)*pow(1+z,2)+Oml );

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  /*return psdfdm(sig8,z,m)*( 1+0.9009*pow(Dg*sig/dc,-0.6) )
   *exp(0.1465*pow(dc/(Dg*sig),2) );*/

  return -0.797885*dc*sig8*dsigdM(m,cosmo)*( 1+1.11*pow(Dg*sig/dc,0.6) )
    *exp(-0.3535*pow(dc/(Dg*sig),2) )/(Dg*sig*sig);
}

/*** derivative of Deltao in CDM model ***/
double dsigdM(double m,const CosmoHndl cosmo){
	double Delta_tmp(double m),err;

	cosmog = cosmo;
	return dfridrD(Delta_tmp,m,0.1*m,&err);
}
double Delta_tmp(double m){
	return Deltao(m,cosmog);
}

/** rms top-hat power in CDM model, normalized to sig8 **/ 
double Deltao(double m,const CosmoHndl cosmo){
  double f4(double),dc;
   double Omo,Oml;

  if(cosmo->physical){
    Omo=cosmo->Omo/cosmo->h/cosmo->h;
    Oml=cosmo->Oml/cosmo->h/cosmo->h;
  }else{
    Omo=cosmo->Omo;
    Oml=cosmo->Oml;
  }

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omo,0.0185);
  if(Omo+Oml==1) dc*=pow(Omo,0.0055);
  /*  return dc*pow(m/1.0e10,-0.25); */
  return f4(6.005e14*pow(cosmo->h*Omo,3))/f4(m*cosmo->h*cosmo->h*cosmo->h*cosmo->h*Omo*Omo);
}

double f4(double u){

  return 8.6594e-12*pow(u,0.67)*pow( 1+pow( 3.5*pow(u,-0.1) +1.628e9*pow(u,-0.63),0.255) ,3.92157);
}

double Ex(double x,const CosmoHndl cosmo){
  if(cosmo->physical) return sqrt(cosmo->Omo*x*x*x+cosmo->Oml+(cosmo->h*cosmo->h-cosmo->Omo+cosmo->Oml)*x*x)/cosmo->h;
  return sqrt(cosmo->Omo*x*x*x+cosmo->Oml+(1-cosmo->Omo+cosmo->Oml)*x*x);
}

/***************************************************************/
/*** isolated integrator used for cosmological calculations ***/
/***************************************************************/

double nintegrateDcos(double (*func)(double),double a,double b,double tols)
{
   void polintD(double xa[], double ya[], int n, double x, double *y, double *dy);
   double trapzdDcoslocal(double (*func)(double), double a, double b, int n);
   void nrerror(char error_text[]);
   double ss,dss;
   double s2[JMAXP],h2[JMAXP+1];
   int j;

   h2[1]=1.0;
   for (j=1;j<=JMAX;j++) {
	s2[j]=trapzdDcoslocal(func,a,b,j);
	if (j>=K) {
	   polintD(&h2[j-K],&s2[j-K],K,0.0,&ss,&dss);
	   if(fabs(dss) <= tols*fabs(ss)) return ss;
	}
	h2[j+1]=0.25*h2[j];
   }
   printf("s2= "); for(j=1;j<=JMAX;j++) printf(" %e ",s2[j]);
   printf("\n");
   nrerror("Too many steps in routine nintegrateDcos");
   return 0.0;
}

double trapzdDcoslocal(double (*func)(double), double a, double b, int n)
{
   double x,tnm,del;
   static double s2,sum2;
   int it,j;

   if (n == 1) {
	return (s2=0.5*(b-a)*(FUNC(a)+FUNC(b)));
   } else {
	for (it=1,j=1;j<n-1;j++) it <<= 1;
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (sum2=0.0,j=1;j<=it;j++,x+=del) sum2 += FUNC(x);
	s2=0.5*(s2+(b-a)*sum2/tnm);
	return s2;
   }
}
