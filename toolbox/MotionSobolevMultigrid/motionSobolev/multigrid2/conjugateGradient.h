#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include <math.h>
#include <stdio.h>
#include "arrayOperations.h"
#include "CGLinearOperator.h"

//using namespace std;

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

 public:
  conjugateGradient() {
    N=0;
    x=r=p=Ap=b=0;
    active=0;
    Aptr=0;
  };

  int allocate(int Size);
  void deallocate();
  void initialize() {
    for (int i=0; i<N; i++) if (active[i]) x[i]=0;
  }

  void computeSolution(double errorTol);

 protected:
  // c = scalea * a + scaleb * b
  void  add(T *a, double scalea, T *b, double scaleb, T* c) {
    for (int i=0; i<N; i++) if (active[i]) c[i]=a[i]*scalea + b[i]*scaleb;
  }
  // a := b
  void  copy(T *a, const T* b) {
    for (int i=0; i<N; i++) if (active[i]) a[i]=b[i];
  }
  // sum_i a[i]*b[i]
  double inner(const T* a, const T* b) {
    double ret=0;
 
    for (int i=0; i<N; i++) if (active[i]) ret+=a[i]*b[i];
    return ret;
  }

  //double error();
  //void print(double *x);
};

#include "conjugateGradient.C"

#endif
