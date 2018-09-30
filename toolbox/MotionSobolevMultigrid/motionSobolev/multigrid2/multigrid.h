#ifndef MULTIGRID_H
#define MULTIGRID_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "arrayOperations.h"
#include "MG_Jacobi.h"
#include "conjugateGradient.h"
#include "vector2D.h"
#include "../params.h"

//using namespace arrayOps;

template <class T> class multigrid {
  
 public:
  
  // Ax=b
  T **x; // grids to store x; x[0]-finest grid
  T **b; // grids to store b; b[0]-finest grid

  // number of levels in the pyramid
  int levels;
  // number of dimensions of the grid
  int dim;
  // sizes of each grid dimension in all levels; sizes[0]-finest grid
  int **sizes; //[MAX_LEVELS][MAX_DIM];
  // increments[l][d] is the value to add to pixel index to move by one-pixel 
  // in level in the positive dimension d direction
  int **increments; //[MAX_LEVELS][MAX_DIM];
  // number of pixels in each grid of the pyramid
  int *GridSize; //[MAX_LEVELS];
  
  // multigrid algorithm parameters
  int steps0;            // number of Vcycle iterations
  int steps1, steps2;    // steps1(2) - relaxation iterations in down-V (up-V)

  // for performing Jacobi relaxation steps (abstract class)
  MG_Jacobi<T> *mgJacobi;
  
  // conjugateGradient for solving in the lowest level (abstract class)
  conjugateGradient<T> *cgSolver;

  T *residual;                 // storing residual in the v-cycle

 protected:

  // temporary arrays for processing

  T *interp_data;              // storing the interpolated data



  // for binary representation
  int **binary_rep;

 public:
  multigrid() {
    x=b=0;
    residual=interp_data=0;
    dim=levels=0;
    mgJacobi=0;
    cgSolver=0;
    sizes=increments=0;
    GridSize=0;
    binary_rep=0;
  }

  int allocate(int dim, const int *sizes, int levels);
  void deallocate();

  // one full v-cycle starting from level start_level
  void vcycle(int start_level);
  // iterate several v-cycles at level start_level
  void solveUsingVCycle(int start_level, double error_tol, int bound=0);

  // full multigrid algorithm (FMG)
  void full_multigrid(double error_tol);
  void initializeFMG();

  // key components of a v-cycle
  template <class TT>
  void restrict_data(const TT *in, TT *out, int level_in);
  virtual void interpolate(const T *in, T *out, int level_in);
  virtual void relax(const T *rhs, T *soln, T* residual, int level, int steps);

  // compute residual at level specified and return its norm
  double computeError(int level);

  // debugging
  void plotImage(double *I, int level);
  void plotImage(vector2D<double> *I, int level, const char *str);

 protected:
  // routines relating to coordinates at various levels of the pyramid

  // get x,y, ... coordinates from the pixel index p
  void getcoords(int p, int level, int *coords) {
    for (int d=dim-1; d>=0; d--) {
      coords[d]=p/increments[level][d];
      p%=increments[level][d];
    }
    return;
  }

  void increment_coords(int *coords, int level) {
    for (int d=0; d<dim; d++) {
      if (coords[d]<sizes[level][d]-1) { coords[d]++; break; }
      else { coords[d]=0; }
    }
    return;
  }

  void init_coords(int *coords) {
    for (int d=0; d<dim; d++) coords[d]=0;
    return;
  }

  // get the coordinates & pixel index for coord in the next finer level
  void coordinate_up( const int *coord, int *coord_up,  int &p_up, int level_up ) {
    p_up=0;
    for (int d=0; d<dim; d++) {
      coord_up[d]=2*coord[d]+1;
      p_up+=coord_up[d]*increments[level_up][d];
    }
    return;
  }

  // get the coordinates & pixel index for coord in the next coarser level
  void coordinate_down( const int *coord, int *coord_down,  int &p_in, 
			int level_down) {
    p_in=0;
    for (int d=0; d<dim; d++) {
      //coord_down[d]=(int)floor( (coord[d]-1)/2.0 );
      coord_down[d]=coord[d]==0 ? -1 : (coord[d]-1)/2;
      p_in+=coord_down[d]*increments[level_down][d];
    }
    return;
  }
  
  // store binary representation for i (w/ bits # of bits, N=2^bits) in bin
  void getbin(int i, int bits, int N, int *bin) {
    N/=2;
    for (int d=bits-1; d>=0; d--) {
      bin[d]=i/N;
      i%=N;
      N/=2;
    }
    return;
  }
  
};

#include "multigrid.C"
#include "debug.C"

#endif
