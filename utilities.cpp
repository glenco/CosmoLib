/*
 * utilities.cpp
 *
 *  Created on: Jan 27, 2012
 *      Author: cgiocoli
 */

#include <utilities.h>

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
