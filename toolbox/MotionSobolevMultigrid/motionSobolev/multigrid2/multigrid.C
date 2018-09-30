#ifndef MULTIGRID_C
#define MULTIGRID_C

#include "multigrid.h"


/*
 * dim - number of dimensions in the grid
 * sizes - array of dim elements with sizes of the dimensions
 * levels - number of grids
 */
template <class T>
int multigrid<T>::allocate(int dim, const int *sizes, int levels)
{
  if (dim<1 || levels<1) return 0;
  
  this->dim=dim;
  this->levels=levels;

#define ALLOCATE_ARR( arr_type, arr, length )				\
  if ( !(arr = new arr_type[length]) ) { deallocate(); return 0; }

  ALLOCATE_ARR( int*, this->sizes, levels );
  ALLOCATE_ARR( int*, increments,  levels );
  ALLOCATE_ARR( int,  GridSize,    levels );

  for (int i=0; i<levels; i++) {
    ALLOCATE_ARR( int, this->sizes[i], dim );
    ALLOCATE_ARR( int, increments[i],  dim );
  }

  for (int i=0; i<dim; i++) this->sizes[0][i]=sizes[i];
  
  for (int l=1; l<levels; l++)
    for (int i=0; i<dim; i++)
      this->sizes[l][i]=this->sizes[l-1][i]/2;
  
  for (int l=0; l<levels; l++) {
    increments[l][0]=1;
    for (int i=1; i<dim; i++)
      increments[l][i]=increments[l][i-1]*this->sizes[l][i-1];
    GridSize[l]=increments[l][dim-1]*this->sizes[l][dim-1];
  }

  ALLOCATE_ARR( T*, x,           levels );
  ALLOCATE_ARR( T*, b,           levels );
  ALLOCATE_ARR( T,  residual,    GridSize[0] );
  ALLOCATE_ARR( T,  interp_data, GridSize[0] );
  
  for (int l=0; l<levels; l++) {
    ALLOCATE_ARR( T,  x[l], GridSize[l] );
    ALLOCATE_ARR( T,  b[l], GridSize[l] );
  }

  ALLOCATE_ARR( int*, binary_rep, (int)pow( (double)2, (double)dim ) );
  for (int b=0; b<(int)pow((double)2,(double)dim); b++)  {
    ALLOCATE_ARR( int, binary_rep[b], dim );
    getbin(b, dim, (int)pow((double)2,(double)dim), binary_rep[b]);
  }

#undef ALLOCATE_ARR
  
  return 1;
}


template <class T>
void multigrid<T>::deallocate()
{

#define DEALLOCATE_ARR( arr ) \
  delete[] arr; arr=0;

  for (int l=0; l<levels; l++) {
    DEALLOCATE_ARR( x[l] );
    DEALLOCATE_ARR( b[l] );
    DEALLOCATE_ARR( sizes[l] ); 
    DEALLOCATE_ARR( increments[l] );
  }

  DEALLOCATE_ARR( x );
  DEALLOCATE_ARR( b );
  DEALLOCATE_ARR( sizes );
  DEALLOCATE_ARR( increments );

  DEALLOCATE_ARR( residual );
  DEALLOCATE_ARR( interp_data );

  levels=dim=0;
  
#undef DEALLOCATE_ARR

  return;
}


/*
 * Restrict data from the fine grid to the coarser grid using the full 
 * weighting restriction.  Applies in any dimensions.
 */
template <class T>
template <class TT>
void multigrid<T>::restrict_data(const TT *data, TT *data_down, int level)
{
  int level_down=level+1;
  int GridSize_down=GridSize[level_down];
  int coord[dim], coord_down[dim];
  int p, dp;
  
  init_coords( coord_down );

#ifdef PARALLELIZE_CODE
    #pragma parallel for private(coord, coord_down)
#endif

//  for (int p_down=0; p_down<GridSize_down; p_down++) {
  for (BLOCKS(block)) {

    BLOCK_ENDS( GridSize_down );

    for (PIXELS(p_down)) {
#ifdef PARALLELIZE_CODE
      getcoords( p_down, level_down, coord_down );
#endif

      coordinate_up( coord_down, coord, p, level );
    
      data_down[p_down]=2*dim*data[p];
      for (int d=0; d<dim; d++) {
        dp=increments[level][d];
        data_down[p_down]+=data[p-dp] +
        ( coord[d]==sizes[level][d]-1 ? data[p] : data[p+dp] );
      }
      data_down[p_down]/=4*dim;

#ifndef PARALLELIZE_CODE
      increment_coords(coord_down, level_down);
#endif
    }
  }
  
  return;
}


/*
 * Interpolate data from coarse grid to finer grid using linear interpolation.
 * Applies in any dimensions.
 */
template <class T>
void multigrid<T>::interpolate(const T *data, T *data_up, int level)
{
  int level_up=level-1;
  int GridSize_up=GridSize[level_up];
  int coord[dim], coord_up[dim];
//  int *bin;
  int N=(int)pow((double)2,(double)dim);

  init_coords( coord_up );

#ifdef PARALLELIZE_CODE
    #pragma parallel for private(coord, coord_up)
#endif

//  for (int p_up=0; p_up<GridSize_up; p_up++) {
  for (BLOCKS(block)) {

    BLOCK_ENDS( GridSize_up );

    for (PIXELS(p_up)) {

#ifdef PARALLELIZE_CODE
      getcoords( p_up, level_up, coord_up );
#endif
      
      int p;

      coordinate_down( coord_up, coord, p, level );
      data_up[p_up]=0;

      int Ninterp=0;

      // march through points on the coarse grid that are part of the box containing
      // p_up in the fine grid
    
      for (int i=0; i<N; i++) {
        //getbin( i, dim, N, bin );
        int *bin=binary_rep[i];

        int max_diff=-1, diff;

        for (int d=0; d<dim; d++) {
          diff=abs(2*bin[d]-(coord_up[d]-2*coord[d]-1));
          max_diff=max_diff > diff ? max_diff : diff;
        }

        if (max_diff<dim) {
          int p_box=0, crdpbin;
	
          for (int d=0; d<dim; d++) {
            crdpbin=coord[d]+bin[d];
            if (crdpbin<sizes[level][d] && crdpbin>=0) {
              p_box+=crdpbin*increments[level][d];
            }
            else {
              p_box=-1;
              break;
            }
          }
	
          if (p_box!=-1) { data_up[p_up]+=data[p_box]; Ninterp++; }
        }
      }
    
      data_up[p_up]=Ninterp==0 ? 0 : data_up[p_up]/Ninterp;

#ifndef PARALLELIZE_CODE
      increment_coords(coord_up, level_up);
#endif
      
    }
  }
  
  return;
}


/*
 * Relaxation using Jacobi iterations.
 *
 * soln - initialization that gets over written with the updated solution
 * rhs  - right hand side
 * residual - residual after completion
 * level - indicates the grid level to indicate size of grid to process
 * steps - number of iterations (steps==0 means use CG until convergence)
 */
template <class T>
void multigrid<T>::relax(const T *rhs, T *soln, T *residual, int level, int steps)
{
  int GridSize_lev=GridSize[level];
  
  if (steps!=0) {
    mgJacobi->jacobiIterate( soln, rhs, steps, level );
  }
  else {
    arrayOps::copy( cgSolver->x, soln, GridSize_lev );
    arrayOps::copy( cgSolver->b, rhs, GridSize_lev );
    
    // solve using conjugate gradient at the coarsest scale
    cgSolver->computeSolution( 1e-6 );
    arrayOps::copy(soln, cgSolver->x, GridSize_lev );
  }

  // compute residual
  mgJacobi->computeResidual( soln, rhs, residual, level );

  return;
}


/*
 * Compute residual at level level and return its l2 norm.
 */
template <class T>
double multigrid<T>::computeError(int level)
{
  mgJacobi->computeResidual( x[level], b[level], residual, level );

  return arrayOps::l2Norm( residual, GridSize[level] );
}


/*
 * Assumes initialization is stored in x[start_level]
 * and b[start_level] has appropriate values already stored
 */
template <class T>
void multigrid<T>::vcycle(int start_level)
{
  // downward part of the V-cycle
  for (int i=start_level; i<levels-1; i++) { 
 
    relax(b[i], x[i], residual, i, steps1);
    restrict_data(residual, b[i+1], i);
    arrayOps::set( x[i+1], 0.0, GridSize[i+1] );
  }

  // solve the coarsest grid until convergence
  relax(b[levels-1], x[levels-1], residual, levels-1, 0);
  
  // upward part of the V-cycle
  for (int i=levels-2; i>=start_level; i--) {
    interpolate(x[i+1], interp_data, i+1);
    arrayOps::add( x[i], x[i], 1, interp_data, 1, GridSize[i] );
    relax(b[i], x[i], residual, i, steps2);
  }

  return;
}


/*
 * Iterates v-cycles for steps0 iterations (if bound==1) or until error does
 * not change much specified by error_tol (if bound==0).
 * 
 * Assumes initialization is stored in x[start_level]
 * and b[start_level] has appropriate value.
 */
template <class T>
void multigrid<T>::solveUsingVCycle(int start_level, double error_tol, int bound)
{
  int iterations=0;
  double error, start_error;

  start_error=computeError(start_level);
  //  printf("start_error=%f\n", start_error);

  do {
    vcycle( start_level );
    error=arrayOps::l2Norm(residual, GridSize[start_level]);
    iterations++;

  //    printf("Iteration=%d, error=%f, exit error=%f\n", iterations, error, 
  //      error_tol*start_error);

    if (bound && iterations==steps0 ) break;

  } while (error>error_tol*start_error);

#ifdef DEBUG
  printf("Vcycle %d, Iterations=%d, start_error=%f, error=%f\n", start_level, iterations, start_error, error);
#endif

  return;
}

  
/*
 * Initializes b for full multigrid algorithm.
 *
 * Assumes that b[0] already has appropriate values stored.
 */
template <class T>
void multigrid<T>::initializeFMG()
{
  for (int i=1; i<levels; i++) {
    restrict_data(b[i-1], b[i], i-1);
  }
  return;
}


/*
 * Full multigrid algorithm.
 *
 * Assumes that b[0] already has appropriate values stored.
 * No initialization for x[0] is assumed.
 */
template <class T>
void multigrid<T>::full_multigrid(double error_tol)
{
  double error;

  //  initializeFMG();

  //  for (int i=0; i<levels; i++) plotImage( b[i], i ,"b" );

  arrayOps::set( x[levels-1], 0.0, GridSize[levels-1] );

  relax(b[levels-1], x[levels-1], residual, levels-1, 0);

#ifdef DEBUG
  printf("Coarse Solving Done; error=%f\n", arrayOps::l2Norm(residual, GridSize[levels-1]) );
#endif
  //  plotImage( x[levels-1], levels-1, "After solving coarsest grid" );

  for (int i=levels-2; i>=0; i--) {    
    interpolate(x[i+1], x[i], i+1);
    solveUsingVCycle(i, error_tol, 1);

  //    plotImage( x[i], i, "After v-cycle" );
  }
  
  return;
}

#endif
