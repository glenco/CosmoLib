/*
 * halo.h
 *
 *  Created on: Feb 2, 2012
 *      Author: cgiocoli
 *
 *      This class implement the properties of an NFW-halo
 */

#ifndef HALO_H_
#define HALO_H_
#include <math.h>
#include <cosmo.h>
#include <utilities.h>

class HALO{
public:
    HALO (COSMOLOGY *co, double mass, double redshift);
    double getRvir(int caseunit=0);
    double getR200();
    double getConcentration(int caseunit=0);
    double getFormationTime(double f=0.5);
    double getFormationRedshift(double f=0.5);
    void reset(double,double);
    virtual ~HALO ();
protected:
    COSMOLOGY *co;   // cosmological model
    double m;        // halo mass
    double z;        // redshift
    double R200;     // radius enclosing 200 times the critical density
    double Rvir;     // radius enclosing the virial overdensity
    double Omz,Omo,Oml;  // matter density parameter at redshift z and at z=0
    double sigma2M;   // the mass variance S(m)
    double deltac0;  // spherical collapse overdensity
private:
    void Set_Parameters();
};
typedef HALO *inCosmoHaloHndl;

#endif /* HALO_H_ */
