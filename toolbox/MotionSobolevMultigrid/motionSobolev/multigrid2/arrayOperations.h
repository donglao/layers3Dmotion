#ifndef ARRAYOPERATIONS_H
#define ARRAYOPERATIONS_H

#include <math.h>


namespace arrayOps {
  
  template<class T> inline static  void copy( T *lhs, const T *rhs, int GridSize ) {
    for (int p=0; p<GridSize; p++) lhs[p]=rhs[p];
    return;
  }
  
  template<class T> inline static void set( T *lhs, T scalar, int GridSize ) {
    for (int p=0; p<GridSize; p++) lhs[p]=scalar;
    return;    
  }

  template<class T> inline static void set( T *lhs, double scalar, int GridSize ) {
    for (int p=0; p<GridSize; p++) lhs[p]=scalar;
    return;    
  }

  template<class T> inline static void add( T *sum, 
					    const T *a, double scale_a, 
					    const T *b, double scale_b, 
					    int GridSize ) {
    for (int p=0; p<GridSize; p++) sum[p]=a[p]*scale_a+b[p]*scale_b;
    return;
  }

  template<class T> inline static void add( T *sum, 
					    const T *a, double scale_a, T b, 
					    int GridSize ) {
    for (int p=0; p<GridSize; p++) sum[p]=a[p]*scale_a+b;
    return;
  }

  template<class T, class S> inline static void 
    mult( T *prod, const T *a, const S *b, int inv, int GridSize ) {
    if (inv) {
      for (int p=0; p<GridSize; p++) prod[p]=a[p]/b[p];
    }
    else {
      for (int p=0; p<GridSize; p++) prod[p]=a[p]*b[p];      
    }
    return;
  }


  template<class T> inline static T mean(const T *a, int GridSize) {
    T sum=0;
    
    for (int p=0; p<GridSize; p++) {
      sum+=a[p];
    }

    return sum / (double)GridSize;
  }


  template<class T> inline static double maxNorm(const T *a, int GridSize) {
    double max_val=-1, val;
    
    for (int p=0; p<GridSize; p++) {
      val=a[p]*a[p];
      max_val=max_val<val ? val : max_val;
    }

    return sqrt(max_val);
  }
  

  template<class T> inline static double l2Norm(const T *a, int GridSize) {
    double sum=0;
    
    for (int p=0; p<GridSize; p++) sum+=a[p]*a[p];
    
    return sqrt(sum/GridSize);
  }


  template<class T> inline static T max(const T *a, int GridSize) {
    T max_val=a[0];
    
    for (int p=1; p<GridSize; p++) {
      max_val=max_val<a[p] ? a[p] : max_val;
    }

    return max_val;
  }


  template<class T> inline static T min(const T *a, int GridSize) {
    T min_val=a[0];
    
    for (int p=1; p<GridSize; p++) {
      min_val=min_val>a[p] ? a[p] : min_val;
    }
    
    return min_val;
  }

  inline static int isdiff(const int *a, const int *b, int GridSize) {
    for (int p=0; p<GridSize; p++) {
      if (a[p]!=b[p]) return 1;
    }
    return 0;
  }

  inline static void getcoords(int p, int dim, const int *increments, int *coords) {
    for (int d=dim-1; d>=0; d--) {
      coords[d]=p/increments[d];
      p%=increments[d];
    }
    return;
  }

  inline static void increment_coords(int dim, const int *sizes, int *coords) {
    for (int d=0; d<dim; d++) {
      if (coords[d]<sizes[d]-1) { coords[d]++; break; }
      else { coords[d]=0; }
    }
    return;
  }

  inline static void init_coords(int dim, int *coords) {
    for (int d=0; d<dim; d++) coords[d]=0;
    return;
  }

}

#endif
