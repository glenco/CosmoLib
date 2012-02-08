/********************************************

  cosmo.cpp calculates some useful comological parameters

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
#include <nr.h>
#include <nrutil.h>
#include <cosmo.h>
#include <utilities.h>

#define JMAX 34
#define JMAXP (JMAX+1)
#define K 6
#define NRANSI
#define CRITD 2.49783e18    /* critical density / Ho^2 solar/Mpc */
#define CRITD2 2.7752543e11 /* critical density / h^2 M_sun/Mpc^3 */

int kmax,kount;
double *xp,**yp,dxsav;
float *xf,*wf;
int ni=513;
double CfactorT;
static double alph;  /* DR-distance parameter */
static double omo, oml, hh;

COSMOLOGY::COSMOLOGY(){
	SetConcordenceCosmology();
	xf=new float[ni];
	wf=new float[ni];
	gauleg(0.,1.,xf,wf,ni);
	fill_linear(vlz,ni,0.,1.7);
	// CfactorT=100.*1.e9*365.25*24.0*3600.0/(3.0856775806e21*1.e3/1.e5);
	CfactorT = 0.102271;
    CfactorT = 1/(h*CfactorT);
	double dc;
	for(int i=0;i<ni;i++){
		double z = -1. + pow(10.,vlz[i]);
		double Omz=Omegam(z);
		if(omo<1 && oml==0) dc=  dc=1.68647*pow(Omz,0.0185);
		if(omo+oml==1) dc=1.68647*pow(Omz,0.0055);
		vDeltaCz.push_back(dc/Dgrowth(z));
		vt.push_back(time(z));
	}
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

	/* default parameterization */
	//if( (physical != 0)*(physical != 1) ) physical = 0;

	physical = 0;

	Omo=0.1358;
	Omb=0.02267;
	h=0.705;

	if(physical==0){ /* use standard parameter set */
		Omo/=h*h;
		Omb/=h*h;
		Oml=1.0-Omo;
	}else if(physical==1){
		Oml=(1-Omo/h/h)*h*h;
	}

	w=-1.0;
	w1=0.0;
	n=1.0;
	Gamma=0.0;
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
 * \brief Print cosmological parameters 
 */
void COSMOLOGY::PrintCosmology(){
	cout << "h: " << h << "\n";
	cout << "n: " << n << "\n";
	cout << "dndlnk: " << dndlnk << "\n";
	cout << "A: " << A << "\n";
	cout << "sig8: " << sig8 << "\n";
	cout << "Gamma: " << Gamma << "\n";

	if(physical==0){
		cout << "Omo: " << Omo << "\n";
		cout << "Oml: " << Oml << "\n";
		cout << "Omb: " << Omb << "\n";
		cout << "Omnu: " << Omnu << "\n";
		cout << "Nnu: " << Nnu << "\n";
  }
  if(physical==1){
		cout << "Omo hh: " << Omo << "\n";
		cout << "Oml hh: " << Oml << "\n";
		cout << "Omb hh: " << Omb << "\n";
		cout << "Omnu hh: " << Omnu << "\n";
		cout << "Nnu: " << Nnu << "\n";
  }
  if(darkenergy==2) cout << "darkenery=" << darkenergy << " gamma=" << gamma << "\n";
  else cout << "darkenery: "<< darkenergy << " w: " << w << " w1: " << w1 << "\n";
}

/** \ingroup cosmolib
 * See if cosmologies are identical
 */
int cosmo_compare(COSMOLOGY *cos1, COSMOLOGY *cos2){

  return 1-(cos1->h == cos2->h)*(cos1->n == cos2->n)*(cos1->getSigma8() == cos2->getSigma8())
		  *(cos1->Omo == cos2->Omo)*(cos1->Oml == cos2->Oml)*(cos1->Omb == cos2->Omb)
		  *(cos1->Omnu == cos2->Omnu)*(cos1->w == cos2->w)*(cos1->w1 == cos2->w1)
		  *(cos1->Gamma == cos2->Gamma)*(cos1->Nnu == cos2->Nnu);
}

/** \ingroup cosmolib
 * \brief Copy a cosmology
 */
void cosmo_copy(CosmoHndl cos1, CosmoHndl cos2){
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

	cos1->power_normalize(cos2->getSigma8());
}

/** \ingroup cosmolib
 * \brief Curvature radius in Mpc
 */
double COSMOLOGY::rcurve(){
  if(Omo+Oml != 1.0){
    if(physical) return 3.0e3/(sqrt(fabs(h*h-Omo-Oml)));  // curvature scale //
    return 3.0e3/(h*sqrt(fabs(1-Omo-Oml)));               // curvature scale // 
  }
  return 0;
}

/** \ingroup cosmolib
 * \brief The expansion function, i.e. the right-hand side of Friedmann's equation
 * without H_0. // radius_dark
 */
double COSMOLOGY::eH(double a){
	if(physical){
		omo=Omo/h/h;
	    oml=Oml/h/h;
	}else{
		omo=Omo;
		oml=Oml;
	}
	return sqrt(omo/pow(a,3)+oml+(1.0-omo-oml)/pow(a,2));
}

double COSMOLOGY::DpropDz(double z){
	double a=1.0/(1.0+z);
	return a/eH(a);
}

/** \ingroup cosmolib
 * \brief Dark matter density parameter at redshift z
 */
double COSMOLOGY::Omegam(double z){
	if(physical){
		omo=Omo/h/h;
	    oml=Oml/h/h;
	}else{
		omo=Omo;
		oml=Oml;
	}
	return omo*pow(1+z,3)/( omo*pow(1+z,3)+(1-omo-oml)*pow(1+z,2)+oml );
}

/** \ingroup cosmolib
 * \brief Comoving Dyer-Roeder angular size distance for lambda=0
 */
double COSMOLOGY::DRradius(double zo,double z,double pfrac){
  double zl,Omr,Omz,Da[2],Ho;
  int nok,nbad;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }

  Ho=h/3.0e3;
  alph=1.5*(1-pfrac);

    /*    if(oml !=0.0){
        printf("ERROR in DRradius: oml != 0\n");
        return 0.0;
        }*/
    zl=1+z;

    if(zo==0 && (oml==0 && alph==0)){
      if(omo==1.0){
	return 2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
      }else if(oml==0){
	Omr=sqrt(1-omo);
	Omz=sqrt(1+omo*(zl-1));
	return zl*( (2-5*omo) + Omz*( omo*(3*zl+2)-2 )/(zl*zl) + 1.5*omo*omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-omo,2));
      }
    }else{
      if((z-zo)<0.001) return (z-zo)/(Ho*(1+zo)*(1+zo)*sqrt(1+omo*zo+oml*(pow(1+zo,-2)-1)) );
      Da[0]=0.0;
      Da[1]=-(1+zo)/(zl*zl*sqrt(1+omo*z+oml*(pow(zl,-2)-1)));

      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
    return 0;
}

/** \ingroup cosmolib
 * \brief Comoving Dyer-Roeder angular size distance
 */
double COSMOLOGY::DRradius2(double zo,double z){
  double zl,Omr,Omz,Da[2],ds,dl,Ho;
  int nok,nbad;

  if(physical){
    omo=omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }
    Ho=h/3.0e3;
    /*    if(oml !=0.0){
        printf("ERROR in DRradius: oml != 0\n");
        return 0.0;
        }*/
    alph=0;
    zl=1+z;

    if(omo==1.0){
      ds=2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
      if(zo==0.0){ return ds;
      }else{
	zl=zo+1;
	dl=2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho);
 	return   ds-dl*(1+z)/(1+zo);
 	/*return   ds*(1+zo)-dl*(1+z);*/
      }
      /*        return 2*( zl*zl-pow(zl,-0.5) )/(zl*5*Ho); */
    }else if(oml==0){
      Omr=sqrt(1-omo);
      Omz=sqrt(1+omo*(zl-1));
      ds= zl*( (2-5*omo) + Omz*( omo*(3*zl+2)-2 )/(zl*zl) + 1.5*omo*omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-omo,2));
      if(zo==0.0){ return ds;
      }else{
	zl=zo+1;
	dl=zl*( (2-5*omo) + Omz*( omo*(3*zl+2)-2 )/(zl*zl) + 1.5*omo*omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-omo,2));
 	return   ds-dl*(1+z)/(1+zo);
	/* 	return   ds*(1+zo)-dl*(1+z); */
	/*return  (1+zo)*( ds-dl );*/
      }
	/*        return zl*( (2-5*omo) + Omz*( omo*(3*zl+2)-2 )/(zl*zl) + 1.5*omo*omo*log( (1+Omr)*(Omz-Omr)/( (1-Omr)*(Omz+Omr) ) )/Omr)/(Ho*4*pow(1-omo,2));*/

    }else{
      if((z-zo)<0.001) return (z-zo)/(Ho*(1+zo)*(1+zo)*sqrt(1+omo*zo+oml*(pow(1+zo,-2)-1)) );
      Da[0]=0.0;
      Da[1]=-1./(zl*zl*sqrt(1+omo*z+oml*(pow(zl,-2)-1)));

      odeintD(Da-1,2,z,zo,1.0e-6,(z-zo)*0.1,0,&nok,&nbad,ders,bsstepD); 
      /*return (1+z)*Da[0]/Ho;*/
      return Da[0]/Ho;
    }
}

void ders(double z,double Da[],double dDdz[]){
  dDdz[1]=Da[2];
  dDdz[2]=-3*Da[2]/(1+z) - (0.5*(omo-2*oml*pow(1+z,-3))*Da[2]+alph*omo*Da[1]/(1+z) )/(1+z*omo+(pow(1+z,-2)-1)*oml);
}

/** \ingroup cosmolib
 * \brief Linear growth factor normalized to 1 at z=0
 */
double COSMOLOGY::Dgrowth(double z){
  double a,Omot,Omlt,g,go;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }
  a=1./(1+z);
  if(omo==1 && oml==0){
	  return a;
  }else{
	  go=2.5*omo/( pow(omo,4.0/7.0)-oml+(1+0.5*omo)*(1+oml/70) );
	  Omot=omo/(omo+oml*a*a*a-a*(omo+oml-1));
	  Omlt=a*a*a*oml*Omot/omo;
	  g=2.5*Omot/( pow(Omot,4.0/7.0)-Omlt+(1+0.5*Omot)*(1+Omlt/70) );
	  return a*g/go;
  }
}

/** \ingroup cosmolib
 * \brief Derivative of the coordinate distance in units of Ho^-1, x=1+z
 */
double COSMOLOGY::radiusm(double x){
 double temp;
  /*printf("->Omo=%f ->Oml=%e x=%e\n",->Omo,->Oml,x);*/
  /*  if( ->Omo==1.0) return 2.0*(1-1/sqrt(x)); */

  temp=omo*x*x*x+oml-(omo+oml-1)*x*x;
  if(temp<=0.0) return -1.0e30;                   // nonphysical values
  return 1.0/sqrt(temp);
}

/** \ingroup cosmolib
 * \brief
 */
double COSMOLOGY::radiusm_dark(double x){
	return 1.0 / sqrt(omo*x*x*x+oml*pow(x,3*(1+w+w1))*exp(-3*w1*(x-1)/x)-(omo+oml-1)*x*x);
}

/** \ingroup cosmolib
 * \brief The coordinate distance in units Mpc
 */
double COSMOLOGY::coorDist(double zo,double z){
	if(physical){
		omo=Omo/h/h;
		oml=Oml/h/h;
	}else{
		omo=Omo;
		oml=Oml;
	}
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
  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }
  if(omo+oml==1) return coorDist(zo,z)/(1+z);
  // curvature scale 
  Rcur=3.0e3/(h*sqrt(fabs(1-omo-oml)));
  if((omo+oml)<1.0) return Rcur*sinh(coorDist(zo,z)/Rcur)/(1+z);
  return Rcur*sin(coorDist(zo,z)/Rcur)/(1+z);
}

/** \ingroup cosmolib
 * \brief The luminosity distance in units Mpc
 */
double COSMOLOGY::lumDist(double zo,double z){
	return pow(1+z,2)*angDist(zo,z);
}

/** \ingroup cosmolib
 * \brief The angular size distance
 */
double COSMOLOGY::gradius(
		double R    /// the curvature radius
		,double rd  /// the coordinate
		){
  /*  printf("Omo=%e Oml=%e (Omo+Oml)-1=*/
  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }

  if( fabs(omo+oml-1) < 1.0e-4) return rd;
  if((omo+oml)<1.0) return R*sinh(rd/R);
  return R*sin(rd/R);
}


/**
  \ingroup cosmo
  \brief Press-Schechter mass function in unit of logarithmic mass bin and mean density
*/
double COSMOLOGY::psdfdm(
		double sig8    /// amplitude of the linear power spectrum on the scale of \f$ 8 h^{-1} Mpc \f$
		,double z      /// redshift
		,double m      /// mass
		,int caseunit  /// if equal to 1 return the number density per unit mass
		){
  double dc,Omz,sig,Dg;

  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }
  hh = h;
  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);
  Omz=Omegam(z);
  dc=1.68647;
  if(omo<1 && oml==0) dc*=pow(Omz,0.0185);
  if(omo+oml==1) dc*=pow(Omz,0.0055);
  /*printf("dc=%e dsigdM=%e sig=%e m=%e D=%e\n",dc,dsigdM(m),sig,m,Dg);*/
  double nu = pow(dc/Dg/sig,2);
  switch (caseunit){
    	  case 1:
    		  return -(sig8*dsigdM(m))/sig*(sqrt(2/M_PI)*sqrt(nu)*exp(-0.5*nu))/m*omo*CRITD2;
    		  break;
    	  default:
    		  return -(sig8*dsigdM(m))/sig*(sqrt(2/M_PI)*sqrt(nu)*exp(-0.5*nu));
    }
}

/** \ingroup cosmolib
  \brief Sheth-Tormen mass function in unit of logarithmic mass bin and mean density
*/
double COSMOLOGY::stdfdm(
		double sig8    /// amplitude of the linear power spectrum on the scale of \f$ 8 h^{-1} Mpc \f$
		,double z      /// redshift
		,double m      /// mass
		,int caseunit  /// if equal to 1 return the number density per unit mass
		){
  double dc,Omz,sig,Dg;
  if(physical){
    omo=Omo/h/h;
    oml=Oml/h/h;
  }else{
    omo=Omo;
    oml=Oml;
  }
  hh = h;
  Dg=Dgrowth(z);
  sig=sig8*Deltao(m);
  Omz=Omegam(z);
  dc=1.68647;
  if(omo<1 && oml==0) dc*=pow(Omz,0.0185);
  if(omo+oml==1) dc*=pow(Omz,0.0055);
  /*return psdfdm(sig8,z,m)*( 1+0.9009*pow(Dg*sig/dc,-0.6) )
   *exp(0.1465*pow(dc/(Dg*sig),2) );*/
  double aST = 0.707;
  double pST = 0.3;
  double AST = 0.322;
  double nu = aST*pow(dc/Dg/sig,2);
  switch (caseunit){
  	  case 1:
  		  return -(sig8*dsigdM(m))/sig*(AST*(1+1/pow(nu,pST))*sqrt(nu/2)*exp(-0.5*nu))/m*omo*CRITD2;
  		  break;
  	  default:
  		  return -(sig8*dsigdM(m))/sig*(AST*(1+1/pow(nu,pST))*sqrt(nu/2)*exp(-0.5*nu));
  }
}

/** \ingroup cosmolib
 * \brief The cumulative comoving number density of haloes with mass larger than m
 * (in solar masses/h) at redshift z; if the dummy argument a is given, the mass function
 * times m^a is integrated. If a is omitted, a default of zero is assumed. The dummy
 * argument t specifies which type of mass function is to be used, PS or ST.
 */
double COSMOLOGY::numberDensity(double sig8,double m, double z, double a, int t){
	double n=0.0;
	for (int i=1;i<ni;i++){
		double y,y1;
		switch (t){
		    case 1: // Sheth-Tormen
		    	y1=log(stdfdm(sig8,z,m/xf[i],1));
		    	break;
		    default: // Press-Schechter
		    	y1=log(psdfdm(sig8,z,m/xf[i],1));
		    	break;
	    }
		y=exp(y1+(a+1.0)*log(m)-(a+2.0)*log(xf[i]));
		n+=wf[i]*y;
	}
	return n;
}

/** \ingroup cosmolib
 * \brief Number of haloes with mass larger than m (in solar masses/h) between
 * redshifts z1 and z2 on the full sky. The dummy argument t specifies which type of
 * mass function is to be used, PS or ST
 */
double COSMOLOGY::number (double sig8, double m, double z1, double z2, int t){
  double n=0.0;
  for (int i=1;i<ni;i++){
    double x=(z2-z1)*xf[i]+z1;
    double d=angDist(0.0,x)/3.0e3*h;
    double f=1.0+x;
    double v=4.0*pi*d*d*DpropDz(x)*f*f*f;
    double c=numberDensity(sig8,m,x,0.0,t);
    n+=wf[i]*v*c;
  }
  return n*(z2-z1)*2.7e10; // Hubble Volume
}

/** \ingroup cosmolib
 * \brief Mass variance \f$ S(m)=\sigma^2(m) \f$
 */
double COSMOLOGY::Variance(double m){
	hh=h ;
	double v = sig8*Deltao(m);
	return v*v;
}


/** \ingroup cosmolib
 * \brief Derivative of \f$ \sigma(m) \f$ in CDM model
 */
double COSMOLOGY::dsigdM(double m){
	double err;
	return dfridrD(Deltao,m,0.1*m,&err);
}

/** \ingroup cosmolib
 * \brief Virial overdensity
 */
double COSMOLOGY:: DeltaV(
		double z         /// redshift
		,int caseunit    /// by default uses the Brayan and Norman fit, if equal to 1 uses the fit by Felix Stšhr
		){
	double omz = Omegam(z);
	double af,bf;
	switch (caseunit){
		case 1:
		// Felix Stšhr fit
			af=0.7076;
			bf=0.4403;
		    if (Oml<1e-4){
		    	af=0.1210;
		        bf=0.6756;
		    }
		    return 0.5*18*M_PI*M_PI*(1+(omz-1)*af+pow(omz,bf));
	    	break;
	    default:
	    // Bryan and Norman fit
	    	double x = omz -1;
	    	if(Oml<1e-4){
	    		return 18*M_PI*M_PI+60*x-32*x*x;
	    	}
	    	else{
	    		return 18*M_PI*M_PI+82*x-39*x*x;
	    	}
	}
}

/** \ingroup cosmolib
 * \brief Return the redshift  given \f$ \delta_c/D_+ \f$
 */
double COSMOLOGY::getZfromDeltaC(double dc){
	  if(dc>vDeltaCz[ni-1]) return -1+pow(10.,vlz[ni-1]);
	  if(dc<vDeltaCz[0]) return -1+pow(10.,vlz[0]);
	  int i = locate (vDeltaCz,dc);
	  i = min (max (i,0), int (ni)-2);
	  // cubic interpolation
	  double f=(dc-vDeltaCz[i])/(vDeltaCz[i+1]-vDeltaCz[i]);
	  double a0,a1,a2,a3,f2;
	  f2 = f*f;
	  a0 = vlz[i+2] - vlz[i+1] - vlz[i-1] + vlz[i];
	  a1 = vlz[i-1] - vlz[i] - a0;
	  a2 = vlz[i+1] - vlz[i-1];
	  a3 = vlz[i];
	  double lzz =(a0*f*f2+a1*f2+a2*f+a3);
	  return -1+pow(10.,lzz);
}

/** \ingroup cosmolib
 * \brief Return the time from the Big Bang in Gyr given
 * \f$ \delta_c/D_+ \f$
 */
double COSMOLOGY::getTimefromDeltaC(double dc){
	  if(dc>vDeltaCz[ni-1]) return vt[ni-1];
	  if(dc<vDeltaCz[0]) return vt[0];
	  int i = locate (vDeltaCz,dc);
	  i = min (max (i,0), int (ni)-2);
	  // cubic interpolation
	  double f=(dc-vDeltaCz[i])/(vDeltaCz[i+1]-vDeltaCz[i]);
	  double a0,a1,a2,a3,f2;
	  f2 = f*f;
	  a0 = vt[i+2] - vt[i+1] - vt[i-1] + vt[i];
	  a1 = vt[i-1] - vt[i] - a0;
	  a2 = vt[i+1] - vt[i-1];
	  a3 = vt[i];
	  return a0*f*f2+a1*f2+a2*f+a3;
}

double COSMOLOGY::timeEarly(double a){
	if(physical){
		omo=Omo/h/h;
		oml=Oml/h/h;
	}else{
		omo=Omo;
	    oml=Oml;
	}
	double aEqual=8.3696e-05/omo;
	double r=aEqual/a;
	return 2.0/3.0*a*sqrt(a)/sqrt(omo)*((1-2*r)*sqrt(1+r)+2*r*sqrt(r));
}

/** \ingroup cosmolib
 * \brief Return the time from the Big Bang in Gyr at a given redshift z
 */
double COSMOLOGY::time(double z){
	double a=1/(z+1.);
	double aEqual=8.3696e-05/omo;
	double e=5.*aEqual;
	if(a>=e){
		double n=0;
		for (int i=1;i<ni;i++){
			double x = e+(a-e)*xf[i];
			double y = 1./x/eH(x);
			n+=wf[i]*y;
		}
		return (n*(a-e)+timeEarly(e))*CfactorT;
	}
	else{
		return CfactorT*timeEarly(a);
	}
}

/** \ingroup cosmolib
 * \brief \f$ \sigma(m) \f$: the rms top-hat power in CDM model, normalized to sig8
 */
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
