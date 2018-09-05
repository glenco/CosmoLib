/*
 * utilities.h
 *
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <math.h>
#include <algorithm>
#include <functional>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>


#ifndef PosType_declare
#define PosType_declare
typedef double PosType;
typedef unsigned long ULONG;
#endif

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


namespace Utilities{
  
  inline void print_date(){
    time_t now = time(0);
    struct tm date = *localtime(&now);
    std::cout << date.tm_hour << ":" << date.tm_min << "   " << date.tm_mday << "/" << date.tm_mon << "/" << date.tm_year + 1900 << std::endl;
  }

  
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
  int long locate (const std::vector<T> &v, const T x)
  {
    size_t n = v.size ();
    if(n==0) return -1;
    
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
  /** \ingroup Utill
   * Locates the element of the given vector which, together with the following
   * element, brackets the given number. If x is smaller than the smallest entry or
   * larger than the largest, the result is either -1 or n-1.
   */
  template <class T,class F>
  int locate (const std::vector<T> &v,F x,std::function<bool (F ,const T &)> less_than)
  {
    size_t n = v.size ();
    
    if (less_than(x,v[0]))
      return -1;
    if (!less_than(x,v[n-1]))
      return n-1;
    
    int jl = -1;
    int ju = n;
    
    
    while (ju-jl > 1)
    {
      int jm = (ju+jl)/2;
      if ( less_than(x,v[jm]) )
        ju=jm;
      else
        jl=jm;
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
  
  
  PosType InterpolateYvec(std:: vector<PosType> &x, std:: vector<PosType> &y,PosType xi);
  
  PosType arctanh(PosType x);
  
  PosType fmini(PosType a,PosType b);
  
  PosType fmaxi(PosType a,PosType b);
  
  PosType median(std:: vector<PosType> vec); 
  
  /// interpolation
  template <typename T>
  void polintT(T xa[], T ya[], int n, T x, T *y, T *dy)
  {
    int i,m,ns=1;
    T den,dif,dift,ho,hp,w;
    T *c,*d;
    
    dif=fabs(x-xa[1]);
    
    c = new T[n+1];
    d = new T[n+1];
    
    for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {
        ns=i;
        dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
        ho=xa[i]-x;
        hp=xa[i+m]-x;
        w=c[i+1]-d[i];
        if ( (den=ho-hp) == 0.0) throw std::runtime_error("Error in routine polint");
        den=w/den;
        d[i]=hp*den;
        c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    delete[] d;
    delete[] c;
  }

  /// used in trapizoidal integral
  template <typename FunctorType,typename T = double>
  T trapz(FunctorType &func, T a, T b, int n, T *s2)
  {
    T x,tnm,del,sum;
    int it,j;
    
    if (n == 1) {
      return (*s2=0.5*(b-a)*( func(a) + func(b) ));
    } else {
      
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += func(x);
      *s2=0.5*(*s2+(b-a)*sum/tnm);
      
      return *s2;
    }
  }

  template <typename FunctorType,typename T = double>
  T nintegrateF(
                FunctorType func        /// struct or class to be integrated
                ,T a      /// limit of integrations
                ,T b      /// limit of integrations
                ,T tols   /// target fractional error
  )
  {
    const int JMAX = 34,K=6;
    const int JMAXP = JMAX+1;
    
    if(a==b) return 0;
    
    T ss,dss;
    T s2[JMAXP],h2[JMAXP+1];
    T ss2=0;
    
    h2[1]=1.0;
    for (int j=1;j<=JMAX;j++) {
      s2[j] = Utilities::trapz<FunctorType,T>(func,a,b,j,&ss2);
       
      if (j>=K) {
        Utilities::polintT<T>(&h2[j-K],&s2[j-K],K,0.0,&ss,&dss);
        if(fabs(dss) <= tols*fabs(ss)) return ss;
      }
      h2[j+1]=0.25*h2[j];
    }
    std::cout << "s2= "; for(int j=1;j<=JMAX;j++) std::cout << s2[j] << "  ";
    std::cout << std::endl << "Too many steps in routine nintegrate<>\n";
    return 0.0;
  }

  /** \brief Numerical integrator.  The class or structure FunctorType must have a () operator that returns
   a number of type T.  Other access functions or public variable can be used instead of global variables
   to change variables within the integrand.
   
   <pre>
   example:
   
   struct FunctorType{
   
   FunctionType(double my_a,double my_b): a(my_a),b(my_b){};
   double a;
   double b;
   
   
   double operator () (double x) { return a*x + b*x*x;}
   }
   
   ............ using it ..................
   
   FunctorType func(21,23.1);
   
   // change some internal variables
   func.a = 3;
   func.b = 10.3;
   // use as a function
   double tmp = func(10);
   
   result = Utilities::nintegrate<FunctorType,double>(func,1.0e-2,1.0e4,1.0e-4);
   
   
   ***  now a double integral ****
   
   \int dy \int dx_y^{10} a*x + b*x*cos(x*y)
   
   struct FunctorType2{
   FunctorType2(FunctorType1 *pfunc1,y): func1(pfunc1),{
   
   FunctorType1 *func1;
   double y;
   
   double operator () (double x) { return func1->a*x + func1->b*x*cos(x*y);}
   }
   
   
   struct FunctorType1{
   FunctionType1(double my_a,double my_b,double xrange1, double xrange2)
   : a(my_a),b(my_b)
   {
   xrange[0] = xrange1;
   xrange[1] = xrange2;
   };
   double a;
   double b;
   double xrange[2];
   
   double operator () (double y) {
   FunctorType2 func2(this);
   func2.a = a;
   func2.b = b;
   func2.y = y;
   xrange[0] = y;
   
   return Utilities::nintegrate<FunctorType2,double>(func2,xrange[0],xrange[1],1.0e-4);
   }
   
   }
   
   ............ using it ..................
   
   FunctorType1 func(21,23.1,0,10.0);
   
   result = Utilities::nintegrate<FunctorType1,double>(func,1.0e-3,5.0e-3,1.0e-4);
   
   <\pre>
   */
  
  template <typename FunctorType,typename T = double>
  T nintegrate(
               FunctorType &func        /// struct or class to be integrated
               ,T a      /// limit of integrations
               ,T b      /// limit of integrations
               ,T tols   /// target fractional error
  )
  {
    const int JMAX = 34,K=6;
    const int JMAXP = JMAX+1;
    
    T ss,dss;
    T s2[JMAXP],h2[JMAXP+1];
    T ss2=0;
    
    h2[1]=1.0;
    for (int j=1;j<=JMAX;j++) {
      s2[j] = Utilities::trapz<FunctorType,T>(func,a,b,j,&ss2);
      if (j>=K) {
        Utilities::polintT<T>(&h2[j-K],&s2[j-K],K,0.0,&ss,&dss);
        if(fabs(dss) <= tols*fabs(ss)) return ss;
      }
      h2[j+1]=0.25*h2[j];
    }
    std::cout << "s2= "; for(int j=1;j<=JMAX;j++) std::cout << s2[j] << "  ";
    std::cout << std::endl << "Too many steps in routine nintegrate<>\n";
    return 0.0;
  }

  
}
#endif /* UTILITIES_H_ */
