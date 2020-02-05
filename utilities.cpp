/*
 * utilities.cpp
 *
 *  Created on: Jan 27, 2012
 *      Author: cgiocoli
 */

#include <utilities.h>

namespace Utilities{
    

PosType arctanh(PosType x){
    return 0.5*log((1+x)/(1-x));
}

PosType fmini(PosType a,PosType b){
  if(a<b)return a;
  return b;
}
PosType fmaxi(PosType a,PosType b){
  if(a>b)return a;
  return b;
}

/**
 * Interpolate (cubic interpolation) the value of a function
 * \f$ y=y(x) \f$ given xi
 */
PosType InterpolateYvec(std:: vector<PosType> &x, std:: vector<PosType> &y,PosType xi){
	int n=x.size();
	if(x[n-1]>x[0]){
		if(xi>x[n-1]) return y[n-1];
		if(xi<x[0]) return y[0];
	}
	if(x[n-1]<x[0]){
		if(xi<x[n-1]) return y[n-1];
		if(xi>x[0]) return y[0];
	}
	int i = locate (x,xi);
	i = std::min (std::max (i,0), int (n)-2);
	PosType f=(xi-x[i])/(x[i+1]-x[i]);
	if(i>1 && i<n-2){
		PosType a0,a1,a2,a3,f2;
	    f2 = f*f;
	    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
	    a1 = y[i-1] - y[i] - a0;
	    a2 = y[i+1] - y[i-1];
	    a3 = y[i];
	    return a0*f*f2+a1*f2+a2*f+a3;
	}
	else return f*y[i+1]+(1-f)*y[i];
}

}
