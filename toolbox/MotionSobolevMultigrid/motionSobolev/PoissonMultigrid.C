#ifndef POISSONMULTIGRID_C
#define POISSONMULTIGRID_C

#include "PoissonMultiGrid.h"


template <class T>
int PoissonMultiGrid<T>::allocate(int dim, const int *sizes, int levels, double _beta)
{
  if ( !multigrid<T>::allocate(dim, sizes, levels) ) { deallocate(); return 0; }

#define ALLOCATE_ARR( arr_type, arr, length )                           \
  if ( !(arr = new arr_type[length]) ) { deallocate(); return 0; }

  ALLOCATE_ARR( double*, region_prob, levels );
  ALLOCATE_ARR( double,  beta, levels );

  for (int i=0; i<levels; i++) {
    ALLOCATE_ARR( double, region_prob[i], this->GridSize[i] );
  }

  setBeta( _beta );

#undef ALLOCATE_ARR

  if ( !(this->mgJacobi=
	 new Poisson_MG_Jacobi<T>( dim, this->sizes, this->increments,
				    this->GridSize, region_prob, beta )) ||
       !(this->cgSolver= new conjugateGradient<T>() ) ||
       !(this->cgSolver->allocate(this->GridSize[levels-1])) ||
       !(this->cgSolver->Aptr=
	 new PoissonMG_CGLinearOperator<T>(dim, this->sizes[levels-1],
					  this->increments[levels-1], 
					  this->GridSize[levels-1], 
					  region_prob[levels-1], beta[levels-1]) ) ||
       !(this->Aptr_dir = new PoissonMGDirichlet_CGLinearOperator<T>(dim, this->sizes[levels-1],
            this->increments[levels-1], 
            this->GridSize[levels-1], 
            region_prob[levels-1], beta[levels-1]) ) )
    {
      deallocate(); return 0;
    }

//  if ( !ccd.allocate( XSize[0], YSize[0]) )
//    { deallocate(); return 0; }

  arrayOps::set( this->cgSolver->active, (char)1, this->GridSize[levels-1] );

  return 1;
}


template <class T>
void PoissonMultiGrid<T>::deallocate()
{
#define DEALLOCATE_ARR( arr )                   \
  delete[] arr; arr=0;

  for (int l=0; l<this->levels; l++) {
    DEALLOCATE_ARR( region_prob[l] );
  }
  DEALLOCATE_ARR( region_prob );
  DEALLOCATE_ARR( beta );

#undef DEALLOCATE_ARR

//  ccd.deallocate();
  delete this->Aptr_dir;

  //delete Aptr; Aptr=0;
  if (this->cgSolver!=0) {
    delete this->cgSolver->Aptr; 
    this->cgSolver->Aptr=0;
    this->cgSolver->deallocate();
  }
  delete this->cgSolver; this->cgSolver=0;
  if (this->mgJacobi!=0) this->mgJacobi->deallocate();
  multigrid<T>::deallocate();

  return;
}


template <class T>
void PoissonMultiGrid<T>::computeRegions(const double *region_prob_fine)
{
  arrayOps::copy(region_prob[0], region_prob_fine, this->GridSize[0]);

  for (int i=1; i<this->levels; i++) {
    this->restrict_data(region_prob[i-1], region_prob[i], i-1);
  }
  
  return;
}


template <class T> 
void PoissonMultiGrid<T>::computeB_CG_Dirichlet(const T *rhs, T *soln, int level)
{
  int GridSize_lev=this->GridSize[level];
  T sum_neighbors;
  int coords[this->dim];
  double *region_prob_lev=region_prob[level];
  T *x = this->cgSolver->x;
  T *b = this->cgSolver->b;
  int *increments_lev = this->increments[level];
  int *sizes_lev = this->sizes[level];
  int dim = this->dim;

  arrayOps::init_coords(dim, coords);


  for (int p=0; p<GridSize_lev; p++) {


    sum_neighbors=0;

    if (region_prob_lev[p]<0.5) {
      for (int d=0; d<dim; d++) {
        int incr = increments_lev[d];

        if (coords[d]<sizes_lev[d]-1) {
            if (region_prob_lev[p+incr]>=0.5) {
              sum_neighbors+=x[p+incr];
            }
        }
        if (coords[d]>0) {
          if (region_prob_lev[p-incr]>=0.5) {
            sum_neighbors+=x[p-incr];
          }
        }
      }
      b[p]=sum_neighbors;
    }

    arrayOps::increment_coords(dim, sizes_lev, coords);
  }

  return;
}


/*
 * Relaxation using Jacobi iterations.
 *
 * Modified to account for extension outside solved via CG
 *
 * soln - initialization that gets over written with the updated solution
 * rhs  - right hand side
 * residual - residual after completion
 * level - indicates the grid level to indicate size of grid to process
 * steps - number of iterations (steps==0 means use CG until convergence)
 */
template <class T>
void PoissonMultiGrid<T>::relax(const T *rhs, T *soln, T *residual, int level, int steps)
{
  int GridSize_lev=this->GridSize[level];
  double *region_prob_lev=region_prob[level];

  if (steps!=0) {
    this->mgJacobi->jacobiIterate( soln, rhs, steps, level );
  }
  else {
    arrayOps::copy( this->cgSolver->x, soln, GridSize_lev );
    arrayOps::copy( this->cgSolver->b, rhs, GridSize_lev );
    
    char *active = this->cgSolver->active;

    for (int p=0; p<GridSize_lev; p++)
      active[p] = region_prob_lev[p]>=0.5;
    // solve using conjugate gradient at the coarsest scale
    this->cgSolver->computeSolution( 1e-6 );
    
    for (int p=0; p<GridSize_lev; p++)
      active[p] = region_prob_lev[p]<0.5;

    computeB_CG_Dirichlet( rhs, soln, level );
    CGLinearOperator<T> *Aptr_tmp = this->cgSolver->Aptr;
    this->cgSolver->Aptr = this->Aptr_dir;
    this->cgSolver->computeSolution( 1e-6 );
    this->cgSolver->Aptr = Aptr_tmp;
  
    arrayOps::copy( soln, this->cgSolver->x, GridSize_lev );
  }

  // compute residual
  this->mgJacobi->computeResidual( soln, rhs, residual, level );

  return;
}


#endif