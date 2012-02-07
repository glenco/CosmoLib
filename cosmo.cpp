/** ******************************************

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


*******************************************
*******************************************/
#include <math.h>
#include <nrD.h>
#include <nrutil.h>
#include <cosmo.h>


#define JMAX 34
#define JMAXP (JMAX+1)
#define K 6
#define NRANSI
#define CRITD 2.49783e18 /* critical density / Ho^2 ; solar/Mpc */


int kmax,kount;
double *xp,**yp,dxsav;

static double alph;  /* DR-distance parameter */
static double omo, oml, hh;

COSMOLOGY::COSMOLOGY(){
	SetConcordenceCosmology();
	// set parameters for Eisenstein&Hu power spectrum
	TFmdm_set_cosm();
}

COSMOLOGY::~COSMOLOGY(){
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
	power_normalize(0.812);
}

/** \ingroup cosmolib
 * \brief
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
 * \brief
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
 * \brief The curvature radius in Mpc
 */
double COSMOLOGY::rcurve(){
  if(Omo+Oml != 1.0){
    return 3.0e3/(h*sqrt(fabs(1-Omo-Oml)));  /** curviture scale **/
  }
  return 0;
}


/** \ingroup cosmolib
 *
 * \brief Comoving Dyer-Roeder angular size distance for lambda=0
***************************************************************/

double COSMOLOGY::DRradius(double zo,double z,double pfrac){
  double zl,Omr,Omz,Da[2],Ho;
  int nok,nbad;

  omo=Omo;
  oml=Oml;

  Ho=h/3.0e3;
  alph=1.5*(1-pfrac);

    /*    if(oml !=0.0){
        printf("ERROR in DRradius: oml != 0\n");
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

  omo=Omo;
  oml=Oml;


    Ho=h/3.0e3;
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

      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
}

void ders(double z,double Da[],double dDdz[]){
  dDdz[1]=Da[2];
  dDdz[2]=-3*Da[2]/(1+z) - (0.5*(omo-2*oml*pow(1+z,-3))*Da[2]+alph*omo*Da[1]/(1+z) )/(1+z*omo+(pow(1+z,-2)-1)*oml);
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

/** \ingroup cosmolib
 * \brief linear growth factor normalized to 1 at z=0
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
 * derivative of the coordinate distance in units of Ho^-1, x=1+z
 * ****/

double COSMOLOGY::radiusm(double x){
 double temp;
  /*printf("->Omo=%f ->Oml=%e x=%e\n",->Omo,->Oml,x);*/
  /*  if( ->Omo==1.0) return 2.0*(1-1/sqrt(x)); */

  temp=Omo*x*x*x+Oml-(Omo+Oml-1)*x*x;
  if(temp<=0.0) return -1.0e30;                   // nonphysical values
  return 1.0/sqrt(temp);
}
/** \ingroup cosmolib
 *  Same as radiusm, but incorporates dark energy through w and w1.
 */
double COSMOLOGY::radiusm_dark(double x){
	return 1.0 / sqrt(Omo*x*x*x+Oml*pow(x,3*(1+w+w1))*exp(-3*w1*(x-1)/x)-(Omo+Oml-1)*x*x);
}

/** \ingroup cosmolib
 * \brief The coordinate distance in units Mpc
 */
double COSMOLOGY::coorDist(double zo,double z){

	if( (w ==-1.0) && (w1 == 0.0) ) return nintegrateDcos(&COSMOLOGY::radiusm,1+zo,1+z,1.0e-9)*3.0e3/h;
	return nintegrateDcos(&COSMOLOGY::radiusm_dark,1+zo,1+z,1.0e-9)*3.0e3/h;
}

/** \ingroup cosmolib
 * \brief The angular size distance in units Mpc
 *
 * Converts angles to proper distance NOT comoving distance.
 */

double COSMOLOGY::angDist(double zo,double z){
  double Rcur;

  if(Omo+Oml==1) return coorDist(zo,z)/(1+z);
   /** curviture scale **/

  Rcur=3.0e3/(h*sqrt(fabs(1-Omo-Oml)));
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
 * \brief Incorporates curvature for angular size distance
 */
double COSMOLOGY::gradius(double R,double rd){
  /*  printf("Omo=%e Oml=%e (Omo+Oml)-1=*/

  if( fabs(Omo+Oml-1) < 1.0e-4) return rd;
  if((Omo+Oml)<1.0) return R*sinh(rd/R);
  return R*sin(rd/R);
}


/** \ingroup cosmolib
* \brief Press-Schechter mass function in fraction of mean density
**/

double COSMOLOGY::psdfdm(double z,double m){
  double dc,Omz,sig,Dg;

  omo=Omo;
  oml=Oml;
  hh = h;

  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);
  Omz=Omo*pow(1+z,3)/( Omo*pow(1+z,3)+(1-Omo-Oml)*pow(1+z,2)+Oml );

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  /*printf("dc=%e dsigdM=%e sig=%e m=%e D=%e\n",dc,dsigdM(m),sig,m,Dg);*/
  return -0.797885*dc*sig8*dsigdM(m)*exp(-0.5*pow(dc/(Dg*sig),2) )
    /(Dg*sig*sig);
}

/** \ingroup cosmolib
 * *
 * \brief Sheth-Tormen mass function, most be normalized to 1 by calling code
 */

double COSMOLOGY::stdfdm(double z,double m){
  double dc,Omz,sig,Dg;

  omo=Omo;
  oml=Oml;
  hh = h;

  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);
  Omz=Omo*pow(1+z,3)/( Omo*pow(1+z,3)+(1-Omo-Oml)*pow(1+z,2)+Oml );

  dc=1.68647;
  if(Omo<1 && Oml==0) dc*=pow(Omz,0.0185);
  if(Omo+Oml==1) dc*=pow(Omz,0.0055);

  /*return psdfdm(sig8,z,m)*( 1+0.9009*pow(Dg*sig/dc,-0.6) )
   *exp(0.1465*pow(dc/(Dg*sig),2) );*/

  return -0.797885*dc*sig8*dsigdM(m)*( 1+1.11*pow(Dg*sig/dc,0.6) )
    *exp(-0.3535*pow(dc/(Dg*sig),2) )/(Dg*sig*sig);
}

/*** derivative of Deltao in CDM model ***/
double dsigdM(double m){
	double err;

	return dfridrD(Deltao,m,0.1*m,&err);
}

/*
**  rms top-hat power in CDM model, normalized to sig8
** **/
double Deltao(double m){
  double dc;

  dc=1.68647;
  if(omo<1 && oml==0) dc*=pow(omo,0.0185);
  if(omo+oml==1) dc*=pow(omo,0.0055);
  /*  return dc*pow(m/1.0e10,-0.25); */
  return f4(6.005e14*pow(hh*omo,3))/f4(m*hh*hh*hh*hh*omo*omo);
}

double f4(double u){
  return 8.6594e-12*pow(u,0.67)*pow( 1+pow( 3.5*pow(u,-0.1) +1.628e9*pow(u,-0.63),0.255) ,3.92157);
}


/***************************************************************/
/*** isolated integrator used for cosmological calculations ***/
/***************************************************************/

double COSMOLOGY::nintegrateDcos(pt2MemFunc func, double a,double b,double tols)
{
   void polintD(double xa[], double ya[], int n, double x, double *y, double *dy);
   //double trapzdDcoslocal(double (*func)(double), double a, double b, int n);
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
   cout << "s2= "; for(j=1;j<=JMAX;j++) cout << s2[j] << " ";
   cout << "\n";
   cout << "Too many steps in routine nintegrateDcos\n";
   return 0.0;
}

double COSMOLOGY::trapzdDcoslocal(pt2MemFunc func, double a, double b, int n)
{
   double x,tnm,del;
   static double s2,sum2;
   int it,j;

   if (n == 1) {
	return (s2=0.5*(b-a)*( (this->*func)(a) +(this->*func)(b) ));
   } else {
	for (it=1,j=1;j<n-1;j++) it <<= 1;
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (sum2=0.0,j=1;j<=it;j++,x+=del) sum2 += (this->*func)(x);
	s2=0.5*(s2+(b-a)*sum2/tnm);
	return s2;
   }
}
