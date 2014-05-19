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

#ifndef PosType_declare
#define PosType_declare
typedef double PosType;
#endif

namespace Utilities{
const double nXbin=64.;

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
 * larger than the largest, the result is either -1 or n-1.
 */
template <class T>
int locate (const std::vector<T> &v, const T x)
{
  size_t n = v.size ();

  if (x == v[0])
    return 0;
  if (x == v[n-1])
    return n-1;

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

  return jl;
}

template <class T>
void locate (T *v, unsigned long n, T x, unsigned long *index)
{
  if (x == v[0]){
    *index = 0;
    return;
  }
  if (x == v[n-1]){
    *index = n-1;
    return;
  }
  
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

  *index = jl;
  return;
}

/// Returns the index of the element of v that is closest to x.  v must be sorted.
template <class T>
size_t closest (const std::vector<T> &v, const T x)
{
  
  int index = locate(v,x);
  if(index == -1) return 0;
  if(index == v.size()-1) return index;
  
  return (x - v[index]) < (v[index+1] - x) ? index : index+1;
}


PosType InterpolateYvec(std:: vector<PosType> x, std:: vector<PosType> y,PosType xi);

PosType arctanh(PosType x);

PosType fmini(PosType a,PosType b);

PosType fmaxi(PosType a,PosType b);

PosType median(std:: vector<PosType> vec); 

}
#endif /* UTILITIES_H_ */
