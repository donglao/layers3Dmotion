#ifndef POISSONMGDIRICHLET_CGLINEAROPERATOR
#define POISSONMGDIRICHLET_CGLINEAROPERATOR

#include "CGLinearOperator.h"

template <class T>
class PoissonMGDirichlet_CGLinearOperator : public CGLinearOperator<T> {
	public:
		int dim;
		int *sizes;
		int *increments;
		int GridSize;

		double *region_prob;
		double beta;

	public:
		PoissonMGDirichlet_CGLinearOperator( int _dim, int *_sizes, int *_increments, 
			int _GridSize, double *_region_prob, double _beta ) : 
		dim(_dim), sizes(_sizes), increments(_increments), GridSize(_GridSize), 
		region_prob(_region_prob), beta(_beta) { };


		void operator() ( const T *x, T *Ax )
		{
			T sum_neighbors;
			int coords[dim];

			arrayOps::init_coords(dim, coords);

			for (int p=0; p<GridSize; p++) {
                //getcoords(p, dim, increments, coords);

				sum_neighbors=0;

				if (region_prob[p]<0.5) {

					for (int d=0; d<dim; d++) {
						int incr=increments[d];

						if (coords[d]<sizes[d]-1) {
							if ( region_prob[p+incr] < 0.5 ) {
								sum_neighbors+= x[p+incr];
							}
						}
						if (coords[d]>0) {
							if ( region_prob[p-incr] < 0.5 ) {
								sum_neighbors+= x[p-incr];
							}
						}
					}

					Ax[p] = x[p]*4.0 - sum_neighbors;
				}

				arrayOps::increment_coords(dim, sizes, coords);
			}

			return;
		}

};

#endif