#ifndef MG_JACOBI_H
#define MG_JACOBI_H

#include "arrayOperations.h"

template <class T>
class MG_Jacobi 
{
 public:
  int dim;
  int **sizes;
  int **increments;
  int *GridSize;

  T *x_update;
  
 public:
  MG_Jacobi( int _dim, int **_sizes, int **_increments, int *_GridSize ) :
  dim(_dim), sizes(_sizes), increments(_increments), GridSize(_GridSize) { 
    allocate( GridSize[0] );
  };
  
  int allocate(int N);
  void deallocate();

  virtual void jacobiStep( const T* x, const T* b, int level ) = 0;
  virtual void computeResidual( const T* x, const T* b, T* residual, 
				int level ) = 0;
  
  void jacobiIterate( T *x, const T*b, int iters, int level );
};

#include "MG_Jacobi.C"

#endif
