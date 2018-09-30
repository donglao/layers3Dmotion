#ifndef POISSONMULTIGRID_H
#define POISSONMULTIGRID_H

#include "multigrid.h"
#include "Poisson_MG_Jacobi.h"
#include "PoissonMG_CGLinearOperator.h"
#include "PoissonMGDirchlet_CGLinearOperator.h"
#include "ConnectedComponentDetection.h"

template <class T>
class PoissonMultiGrid : public multigrid<T> {
 public:
  double **region_prob;            // probability that pixel belongs to region R
	//  ConnectedComponentDetection ccd; // for detecting connected components
  double *beta;                   // tweak parameter to avoid connected components

  CGLinearOperator<T> *Aptr_dir;

 public:
  PoissonMultiGrid() {
  	region_prob=0;
    beta=0;
  }

  int allocate(int dim, const int *sizes, int levels, double _beta);
  void deallocate();

  void computeRegions(const double *region_prob_fine);
  void computeB_CG_Dirichlet(const T *rhs, T *soln, int level);

//  void interpolate(const T *data, T *data_up, int level);
  void relax(const T *rhs, T *soln, T *residual, int level, int steps);

  void setBeta(double _beta) {
  	beta[0]=_beta;
  	for (int i=1; i<this->levels; i++)
  		beta[i] = beta[i-1]*4.0;
  	return;
  }

};

#include "PoissonMultiGrid.C"

#endif
