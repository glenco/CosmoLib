/*
 * halo.h
 *
 *  Created on: Oct 23, 2012
 *      Author: cgiocoli
 *
 *      This class implement the non-linear power spectrum using the Halo Model
 */

#ifndef POWERCDMHM_H_
#define POWERCDMHM_H_
#include <math.h>
#include <cosmo.h>
#include <halo.h>
#include <utilities.h>

#ifdef ENABLE_GSL

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>

/** 
 * \brief Class for calculating the non-linear power spectrum using ****** Halo Model.
 *  To use this class you must enable GSL.
 */
class POWERCDMHM{
  public:
  POWERCDMHM(COSMOLOGY *co, double redshift, double minHalomass, double sigma=0.25,int cmRelation=0,double slope=-0.1);
  
  double nonlinpowerCDMHM(double k);
  double nonlinpowerCDMHM1Halo(double k);
  double nonlinpowerCDMHM2Halo(double k);
  double nonlinKAPPApowerCDMHM(double l, double zs);
  double nonlinKAPPApowerCDMHM1Halo(double l, double zs);
  double nonlinKAPPApowerCDMHM2Halo(double l, double zs);
  double linKAPPApowerCDMHM(double l, double zs);
  double nonlinfitKAPPApowerCDMHM(double l, double zs);
  virtual ~POWERCDMHM ();
  
  protected:
  gsl_function intPk1,intPk2;
  double Pklin,Pk1,Pk2;
  double Pk20;
  double weight (double z1, double z2);
  double weight (double z0);

  struct weightfnc{
    std::vector<double> ai;
    std::vector<double> wi;
  };
  weightfnc wgf; // weight function
  void Initweight(double zs);
  private:
  float *xf,*wf;
};
#endif


#endif /* POWERCDMHM_H_ */




