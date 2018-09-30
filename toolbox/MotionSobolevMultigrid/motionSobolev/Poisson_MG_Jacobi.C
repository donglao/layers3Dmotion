#ifndef POISSON_MG_JACOBI_C
#define POISSON_MG_JACOBI_C

#include "Poisson_MG_Jacobi.h"

template <class T>
void Poisson_MG_Jacobi<T>::jacobiStep( const T *x, const T *b, int level )
{
  int GridSize_lev=this->GridSize[level];
  int Nneighbors;
  T sum_neighbors;
  int coords[this->dim];
  double *region_prob_lev=region_prob[level];
  double dt=0.2;
  double beta_lev = beta[level];
  int *increments_lev = this->increments[level];
  int *sizes_lev = this->sizes[level];
  T *x_update = this->x_update;
  int dim = this->dim;

  arrayOps::init_coords(dim, coords);

#ifdef PARALLELIZE_CODE
    #pragma parallel for private(Nneighbors, sum_neighbors, coords)
#endif

  //  for (int p=0; p<GridSize_lev; p++) {
  for (BLOCKS(block)) {

    BLOCK_ENDS( GridSize_lev );

    for (PIXELS(p)) {

#ifdef PARALLELIZE_CODE
      arrayOps::getcoords(p, dim, increments_lev, coords);
#endif

      sum_neighbors=Nneighbors=0;

      if ( region_prob_lev[p] >=0.5 ) {
        for (int d=0; d<dim; d++) {
          int incr = increments_lev[d];

          if (coords[d]<sizes_lev[d]-1) {
            if ( region_prob_lev[p+incr] >= 0.5 ) {
              Nneighbors++;
              sum_neighbors+=x[p+incr];
            }
          }
          if (coords[d]>0) {
            if ( region_prob_lev[p-incr] >= 0.5 ) {
              Nneighbors++;
              sum_neighbors+=x[p-incr];
            }
          }
        }

  //      diagA= Nneighbors > 0 ? Nneighbors : 1;
  //      this->x_update[p]= (b[p] + sum_neighbors)/diagA;
        x_update[p]= (1 - double(Nneighbors)*dt - beta_lev*dt)*x[p] + ( sum_neighbors + b[p] )*dt;
      }
      else {
        for (int d=0; d<dim; d++) {
          int incr = increments_lev[d];

          if (coords[d]<sizes_lev[d]-1) {
            Nneighbors++;
            sum_neighbors+=x[p+incr];
          }

          if (coords[d]>0) {
            Nneighbors++;
            sum_neighbors+=x[p-incr];
          }
        }

  //      diagA= Nneighbors > 0 ? Nneighbors : 1;
  //      this->x_update[p] = sum_neighbors/diagA;
        x_update[p]= (1 - double(Nneighbors)*dt)*x[p] + ( sum_neighbors )*dt;
  //      this->x_update[p]=0;
      }

#ifndef PARALLELIZE_CODE
      arrayOps::increment_coords(dim, sizes_lev, coords);
#endif

    }
  }
  
  return;
}


template <class T>
void Poisson_MG_Jacobi<T>::computeResidual( const T* x, const T* b, T* residual, int level )
{
  int GridSize_lev=this->GridSize[level];
  int Nneighbors;
  T sum_neighbors;
  int coords[this->dim];
  double *region_prob_lev=region_prob[level];
  double beta_lev = beta[level];
  int *increments_lev = this->increments[level];
  int *sizes_lev = this->sizes[level];
  int dim = this->dim;

  arrayOps::init_coords(dim, coords);

#ifdef PARALLELIZE_CODE
    #pragma parallel for private(Nneighbors, sum_neighbors, coords)
#endif


//  for (int p=0; p<GridSize_lev; p++) {
  for (BLOCKS(block)) {

    BLOCK_ENDS( GridSize_lev );

    for (PIXELS(p)) {

#ifdef PARALLELIZE_CODE    
      arrayOps::getcoords(p, dim, increments_lev, coords);
#endif

      sum_neighbors=Nneighbors=0;

      if ( region_prob_lev[p] >= 0.5 ) {
        for (int d=0; d<dim; d++) {
          int incr = increments_lev[d];

          if (coords[d]<sizes_lev[d]-1) {
            if ( region_prob_lev[p+incr] >= 0.5 ) {
              Nneighbors++;
              sum_neighbors+=x[p+incr];
            }
          }
          if (coords[d]>0) {
            if ( region_prob_lev[p-incr] >= 0.5 ) {
              Nneighbors++;
              sum_neighbors+=x[p-incr];
            }
          }
        }
        residual[p]=b[p] + sum_neighbors - x[p]*Nneighbors - x[p]*beta_lev;
      }
      else {
        for (int d=0; d<this->dim; d++) {
          int incr = increments_lev[d];

          if (coords[d]<sizes_lev[d]-1) {
            Nneighbors++;
            sum_neighbors+=x[p+incr];
          }
          if (coords[d]>0) {
            Nneighbors++;
            sum_neighbors+=x[p-incr];
          }
        }
        residual[p]=sum_neighbors - x[p]*Nneighbors;
      }

#ifndef PARALLELIZE_CODE
      arrayOps::increment_coords(dim, sizes_lev, coords);
#endif

    }
  }
  
  return;
}

#endif
