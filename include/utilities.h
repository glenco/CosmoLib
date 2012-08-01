/*
 * utilities.h
 *
 *  Created on: Jan 27, 2012
 *      Author: cgiocoli
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <math.h>
#include <algorithm>

/** \ingroup Utill
 * Fills a vector with equidistant points from [min, max].
 */
template <class T>
void fill_linear ( std::vector<T> &v, size_t n, T min, T max ){
  v.resize ( n );
  for ( size_t i = 0; i < n; i++ )
    v[i] = min + (max - min) * T (i)/T (n-1);
}

/** \ingroup Utill
 * Fills a vector with logarithmically equidistant points from [min, max].
 */
template <class T>
void fill_logarithmic ( std::vector<T> &v, size_t n, T min, T max ){
  v.resize (n);
  for ( size_t i = 0; i < n; i++ )
    v[i] = exp ( log (min) + ( log (max) - log (min) ) * T (i)/T (n-1) );
}

/** \ingroup Utill
 * Locates the element of the given vector which, together with the following
 * element, brackets the given number. If x is smaller than the smallest entry or
 * larger than the largest, the result is either -1 or n.
 */
template <class T>
int locate (const std::vector<T> &v, const T x)
{
  size_t n = v.size ();
  int jl = -1;
  int ju = n;
  bool as = (v[n-1] >= v[0]);
  while (ju-jl > 1)
  {
    int jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}

/** \ingroup Utill
 * \brief Template function that maps values of inputs
 */
template <class T>
void swap (T a,T b){
	T tmp;
	tmp = a;
	a = b;
	b = tmp;
}

double InterpolateYvec(std:: vector<double> x, std:: vector<double> y,double xi);

double arctanh(double x);

double fmini(double a,double b);

double fmaxi(double a,double b);

double median(std:: vector<double> vec); 

#endif /* UTILITIES_H_ */
