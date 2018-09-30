#ifndef POISSON_MG_JACOBI_H
#define POISSON_MG_JACOBI_H

#include "MG_Jacobi.h"

template <class T>
class Poisson_MG_Jacobi : public MG_Jacobi<T> {
 public:
  double **region_prob;
  double *beta;

 public:
  Poisson_MG_Jacobi( int _dim, int **_sizes, int **_increments, 
		      int *_GridSize, double **_region_prob, double *_beta ) :
  MG_Jacobi<T>( _dim, _sizes, _increments, _GridSize ), 
    region_prob(_region_prob), beta(_beta) { };
  
  void jacobiStep( const T* x, const T* b, int level );
  void computeResidual( const T* x, const T* b, T* residual, int level );
};

#include "Poisson_MG_Jacobi.C"

#endif