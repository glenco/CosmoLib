/**********************************
nfwdensity.c calculates the central 
density of a NFW halo of mass M

*** rc is in Megaparsecs!!!
***
***
***
***   NOT YET MADE COMPATIBLE WITH OTHER CosmoLib ROUTINES
**********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "cosmo.h"
#include <cosmo.h>

//extern struct cosmology cosmo;

void nfwdensity(double sig8,double M,double zo,double *dc,double *rc,double *cons,CosmoHndl cosmo){
  int n;
  double dcrito,doverd,f,Omz,r200,c;
  double zb,yb,zt,yt,zm,ym;
  double dcofc(double c),dcrit(double z,const CosmoHndl cosmo);

  //double Dgrowth(double),Deltao(double),dcrit(double),dcofc(double);

  if(cosmo->physical == 1){
	printf("ERROR: nfwdensity not updated for physical parameters");
	exit(1);
  }
  f=0.01;

  if(cosmo->Omo+cosmo->Oml != 1 && cosmo->Oml != 0){
    printf("\nERROR: cosmology must be flat or Lambda=0 in nfwdensity"); exit(0);}

  /*** critical density at zo***/
  dcrito=dcrit(zo,cosmo);
  /*** Omega at zo***/
  Omz=cosmo->Omo*pow(1+zo,3)/( cosmo->Omo*pow(1+zo,3)+(1-cosmo->Omo-cosmo->Oml)*pow(1+zo,2)+cosmo->Oml );

  doverd=1+0.477*sig8*sqrt( 2*(pow(Deltao(f*M,cosmo),2)-pow(Deltao(M,cosmo),2)) )/dcrito;
  if(doverd<1){ *dc=0; *rc=0; exit(0);} /*** halo has not yet formed **/

  /*printf("doverd=%e\n",doverd); */

  /* find zcoll which gives the fraction doverd */

  /*initial values*/
  zb=zo;  yb=1-doverd;
  zt=100; yt=dcrit(zt,cosmo)/dcrito-doverd;

  n=0;
  while(fabs(ym)>1.0e-10 || n==0 ){
    zm=zt-yt*(zt-zb)/(yt-yb);
    ym=dcrit(zm,cosmo)/dcrito-doverd;
    if(ym>0){zt=zm; yt=ym;}
    if(ym<0){zb=zm; yb=ym;}
    ++n;
    /*printf("zb=%f yb=%e, zt=%f yt=%e, zm=%f ym=%e\n",zb,yb,zt,yt,zm,ym);*/
  }
  /*printf("zb=%f yb=%e, zt=%f yt=%e, zm=%f ym=%e n=%i\n",zb,yb,zt,yt,zm,ym,n);*/
  /*printf("zm=%f n=%i\n",zm,n);*/

  *dc=3.41e3*Omz*pow((1+zm)/(1+zo),3);

  /**** find scale length ****/
  r200=1.63e-5*pow(M*cosmo->h*Omz/cosmo->Omo,0.33333)/(cosmo->h*(1+zo));
  /*initial values, z is realy c */
  zb=1; yb=dcofc(zb)-*dc;
  zt=100; yt=dcofc(zt)-*dc;

  /*printf("zb=%f yb=%e, zt=%f yt=%e\n",zb,yb,zt,yt);*/

  n=0;
  while(fabs(ym)>1.0e-9 || n==0 ){
    c=zt-yt*(zt-zb)/(yt-yb);
    /*printf("c=%e\n",c); */
    ym=dcofc(c)-*dc;
    /*    if(ym>0){zt=c; yt=ym;}   // false position method *
    if(ym<0){zb=c; yb=ym;} */
    if( fabs(yb)>fabs(yt) ){zb=c; yb=ym;}  /** scant method **/
    else{zt=c; yt=ym;}
    ++n;
  /*printf("zb=%f yb=%e, zt=%f yt=%e, c=%f ym=%e n=%i\n",zb,yb,zt,yt,c,ym,n);*/
  }
  /*printf("zb=%f yb=%e, zt=%f yt=%e\n",zb,yb,zt,yt);*/
  /*printf("c=%f n=%i\n",c,n);*/
  *rc=r200/c;
  *cons=c;
}
double dcofc(double c){
  return 200*c*c*c*(1+c)/(3*((1+c)*log(1+c)-c));
}

/** redshift dependent critical density **/

double dcrit(double z,const CosmoHndl cosmo){
  double ans,Omz;

  Omz=cosmo->Omo*pow(1+z,3)/( cosmo->Omo*pow(1+z,3)+(1-cosmo->Omo-cosmo->Oml)*pow(1+z,2)+cosmo->Oml );

  ans=1.68647/Dgrowth(z,cosmo);
  if(cosmo->Omo<1 && cosmo->Oml==0) ans*=pow(Omz,0.0185);
  if(cosmo->Omo+cosmo->Oml==1) ans*=pow(Omz,0.0055);

  return ans;
}

/*** the mass in a NFW profile within R ***/
double nfwmass(double dc,double rs,double R,const CosmoHndl cosmo){
  double rhocrit,x;

  rhocrit=2.775e11*cosmo->h*cosmo->h;
  x=R/rs;

  return 4*pi*dc*rhocrit*pow(rs,3)*( log(1+x)-x/(1+x) );
}
