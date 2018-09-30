#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include <math.h>
#include <stdio.h>
#include "CGLinearOperator.h"

#ifdef PARALLELIZE_CODE
#include <omp.h>
#endif

#define BLOCKS(block) int block=0; block<Nblocks; block++
#define PIXELS(p) int p=start; p<end; p++ 
#define BLOCK_ENDS(N)                                         \
  int start = block*(N)/Nblocks;                              \
  int end   = block==Nblocks-1 ? (N) : (block+1)*(N)/Nblocks


/*
 * Solves the system Ax = b using the conjugate gradient method.
 * Assumes A is linear, symmetric, positive definite.
 *
 * The matrix A is not stored, instead the class uses a function pointer
 * such that when the function pointed to is called, Ax is computed.
 */
template <class T>
class conjugateGradient
{
 public:
  int N; //size of matrix

  T *x, *r, *p, *Ap, *b;
  
  char *active;

  double (*E)(double a);
  CGLinearOperator<T> *Aptr;

  // for parallel
  int Nblocks;

 public:
  conjugateGradient() {
    N=0;
    x=r=p=Ap=b=0;
    active=0;
    Aptr=0;
    Nblocks=1;
  };

  int allocate(int Size, int Nblocks);
  void deallocate();
  void initialize() {
    for (int i=0; i<N; i++) if (active[i]) x[i]=0;
  }

  void computeSolution(double errorTol);

 protected:
  // c = scalea * a + scaleb * b
  void  add(T *a, double scalea, T *b, double scaleb, T* c) {
#ifdef PARALLELIZE_CODE
  #pragma omp parallel for
#endif

    for (BLOCKS(block)) {
      BLOCK_ENDS( N );
      for (PIXELS(i)) if (active[i]) c[i]=a[i]*scalea + b[i]*scaleb;
    }
  }

  // a := b
  void  copy(T *a, const T* b) {
#ifdef PARALLELIZE_CODE
    #pragma omp parallel for 
#endif

    for (BLOCKS(block)) {
      BLOCK_ENDS( N );
      for (PIXELS(i)) if (active[i]) a[i] = b[i];
    }
  }

  // sum_i a[i]*b[i]
  double inner(const T* a, const T* b) {
    double ret=0;
    double sum[Nblocks];

#ifdef PARALLELIZE_CODE
  #pragma omp parallel for
#endif

    for (BLOCKS(block)) {
      BLOCK_ENDS(N);

      double sum_tmp=0;
      for (PIXELS(i)) if (active[i]) sum_tmp+=a[i]*b[i];

      sum[block]=sum_tmp;
    }

    for (BLOCKS(block)) ret += sum[block];

    return ret;
  }

};

#include "conjugateGradient.C"

#endif