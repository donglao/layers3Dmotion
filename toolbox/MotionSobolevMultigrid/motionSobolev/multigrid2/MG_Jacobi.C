#ifndef MG_JACOBI_C
#define MG_JACOBI_C

#include "MG_Jacobi.h"

template <class T>
int MG_Jacobi<T>::allocate(int N)
{
  if ( !(x_update = new T[N]) ) return 0;
  return 1;
}


template <class T>
void MG_Jacobi<T>::deallocate()
{
  delete[] x_update; x_update=0;
}


template <class T>
void MG_Jacobi<T>::jacobiIterate( T *x, const T *b, int iters, int level )
{
  while (iters--) {
    jacobiStep( x, b, level );
    arrayOps::copy( x, x_update, GridSize[level] );
  }
  
  return;
}

#endif
