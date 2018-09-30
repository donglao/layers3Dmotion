#ifndef CONJUGATEGRADIENT_C
#define CONJUGATEGRADIENT_C

#include "conjugateGradient.h"


template <class T>
int conjugateGradient<T>::allocate(int Size)
{
  if ( Size <= 0 ) return 0;
  N=Size;

  if ( !(x  = new T[N]) ||
       !(r  = new T[N]) ||
       !(p  = new T[N]) ||
       !(Ap = new T[N]) ||
       !(b  = new T[N]) ||
       !(active = new char[N]) ) {
    deallocate();
    return 0;
  }

  return 1;
}


template <class T>
void conjugateGradient<T>::deallocate()
{
  delete[]  x;  x=0;
  delete[]  r;  r=0;
  delete[]  p;  p=0;
  delete[] Ap; Ap=0;
  delete[]  b;  b=0;
  delete[] active; active=0;

  return;
}


template <class T>
void conjugateGradient<T>::computeSolution(double errorTol)
{
  double alpha, beta;
  double err;
  T *Ax=r;
  CGLinearOperator<T> &A=*Aptr;
  double tol;

  double bb=inner(b,b);
  tol=bb>0 ? errorTol*errorTol*bb : errorTol;

  //initialize();                  // x = 0
  int numbActive=0;
  for (int i=0; i<N; i++) if (active[i]) numbActive++;

  A(x, Ax);                      // Ax = A(x)
  add(Ax, -1, b, 1, r);          // r = b - Ax
  copy(p, r);                    // p = r

  double rr=inner(r,r), rrnew;
  int count=0;

  double rr_start = rr;

  while (rr > tol && count < 2*numbActive) {
    
    A(p, Ap);                    // Ap = A(p)
    alpha=rr/inner(p,Ap);        // alpha = r.*r / p .* Ap
    //printf("iter.=%d/%d, err=%f, errorTol=%f, alpha=%f\n", 
    //	   count, numbActive, rr, tol, alpha);
    add( x, 1,  p,  alpha, x);   // x <= x + alpha *  p
    add( r, 1, Ap, -alpha, r);   // r <= r - alpha * Ap

    //err=error();

    rrnew=inner(r,r);
    beta=rrnew/rr;
    add( r, 1, p, beta, p);      // p <= r + beta*p
    rr=rrnew;

    count++;
  }

#ifdef DEBUG
  printf("iter.=%d/%d, err_start=%f, err=%f, errorTol=%f\n", count, numbActive, rr_start, rr, tol);
#endif

  return;
}

#endif
