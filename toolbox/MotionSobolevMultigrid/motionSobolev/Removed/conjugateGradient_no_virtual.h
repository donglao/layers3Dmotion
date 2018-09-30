#ifndef CONJUGATEGRADIENT_NO_VIRTUAL_H
#define CONJUGATEGRADIENT_NO_VIRTUAL_H

#include <math.h>
#include <stdio.h>

#ifdef PARALLELIZE_CODE
#include <omp.h>
#endif

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

  // for parallel
  int Nblocks;


  // for Ax operations for Poisson eqn
  int XSize, YSize, GridSize;
  char *region;
  unsigned char *edgebits;
  int poisson;
  enum {XPOS=1,XNEG=2,YPOS=4,YNEG=8,XEDGES=3,YEDGES=12,EDGES=15};

 public:
  conjugateGradient() {
    N=0;
    x=r=p=Ap=b=0;
    active=0;
    Nblocks=1;
    poisson=1;
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
    int start, end;

#ifdef PARALLELIZE_CODE
  #pragma omp parallel for private(start,end)
#endif

    for (int block=0; block<Nblocks; block++) {
      start = block*N/Nblocks;
      end   = block==Nblocks-1 ? N : (block+1)*N/Nblocks;

      for (int i=start; i<end; i++) if (active[i]) c[i]=a[i]*scalea + b[i]*scaleb;
    }
  }

  // a := b
  void  copy(T *a, const T* b) {
    int start, end;

    #pragma omp parallel for private(start,end)

    for (int block=0; block<Nblocks; block++) {
      start = block*N/Nblocks;
      end   = block==Nblocks-1 ? N : (block+1)*N/Nblocks;
      for (int i=start; i<end; i++) if (active[i]) a[i]=b[i];
    }
  }

  // sum_i a[i]*b[i]
  double inner(const T* a, const T* b) {
    double ret=0;
    int start, end;
    double sum[Nblocks];

#ifdef PARALLELIZE_CODE
  #pragma omp parallel for private(start,end)
#endif

    for (int block=0; block<Nblocks; block++) {
      start = block*N/Nblocks;
      end= block==Nblocks-1 ? N : (block+1)*N/Nblocks;

      double sum_tmp=0;
      for (int i=start; i<end; i++) if (active[i]) sum_tmp+=a[i]*b[i];

      sum[block]=sum_tmp;
    }

    for (int block=0; block<Nblocks; block++) ret += sum[block];

    return ret;
  }

    // routines for Poisson equation
  void Ax_Poisson( const T* x, T* Ax) {
    const T *xc=x;
    const T *yc=x+GridSize;
    T *Ax_GridSize = Ax + GridSize;
    T G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my;
    int start, end;

#ifdef PARALLELIZE_CODE
      #pragma omp parallel for private(G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my, start, end)
#endif

    for (int block=0; block<Nblocks; block++) {
      start = block*GridSize/Nblocks;
      end= block==Nblocks-1 ? GridSize : (block+1)*GridSize/Nblocks;

      for (int p=start; p<end; p++) {
        if (region[p]) {
          if ( !(edgebits[p]&XPOS) && region[p+1] ) {
            G1px = xc[p+1] - xc[p];
            G2px = yc[p+1] - yc[p];
          }
          else {
            G1px = G2px = 0;
          }
          if ( !(edgebits[p]&XNEG) && region[p-1] ) {
            G1mx = xc[p-1] - xc[p]; 
            G2mx = yc[p-1] - yc[p];
          }
          else {
            G1mx = G2mx = 0;
          }

          if ( !(edgebits[p]&YPOS) && region[p+XSize] ) {
            G1py = xc[p+XSize] - xc[p];
            G2py = yc[p+XSize] - yc[p];
          }
          else {
            G1py = G2py = 0;
          }
          if ( !(edgebits[p]&YNEG) && region[p-XSize] ) {
            G1my = xc[p-XSize] - xc[p];
            G2my = yc[p-XSize] - yc[p];
          }
          else {
            G1my = G2my = 0;
          }

          Ax[p] = -(G1px + G1mx + G1py + G1my );
          Ax_GridSize[p] = -(G2px + G2mx + G2py + G2my );
        }
      }
    }

    return;
  }

  void Ax_Poisson_Dirichlet( const T* x, T* Ax ) {
    const T *xc=x;
    const T *yc=x+GridSize;
    T *Ax_GridSize = Ax + GridSize;
    T G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my;
    int start, end;
    
#ifdef PARALLELIZE_CODE
    #pragma omp parallel for private(G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my, start, end)
#endif

    for (int block=0; block<Nblocks; block++) {
      start = block*GridSize/Nblocks;
      end= block==Nblocks-1 ? GridSize : (block+1)*GridSize/Nblocks;

      for (int p=start; p<end; p++) {
        if (!region[p]) {
          if ( !(edgebits[p]&XPOS) && !region[p+1] ) {
            G1px = xc[p+1];
            G2px = yc[p+1];
          }
          else {
            G1px = G2px = 0;
          }
          if ( !(edgebits[p]&XNEG) && !region[p-1] ) {
            G1mx = xc[p-1];
            G2mx = yc[p-1];
          }
          else {
            G1mx = G2mx = 0;
          }
          if ( !(edgebits[p]&YPOS) && !region[p+XSize] ) {
            G1py = xc[p+XSize];
            G2py = yc[p+XSize];
          }
          else {
            G1py=G2py=0;
          }
          if ( !(edgebits[p]&YNEG) && !region[p-XSize] ) {
            G1my = xc[p-XSize];
            G2my = yc[p-XSize];
          }
          else {
            G1my = G2my = 0;
          }

          Ax[p]          = -(G1px + G1mx + G1py + G1my) + 4*xc[p];
          Ax_GridSize[p] = -(G2px + G2mx + G2py + G2my) + 4*yc[p];
        }
      }
    }

    return; 
  }

};

#include "conjugateGradient_no_virtual.C"

#endif
