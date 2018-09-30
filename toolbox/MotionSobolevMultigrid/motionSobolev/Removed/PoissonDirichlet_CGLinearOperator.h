#ifndef POISSONDIRICHLET_CGLINEAROPERATOR
#define POISSONDIRICHLET_CGLINEAROPERATOR

#include "CGLinearOperator.h"

class PoissonDirichlet_CGLinearOperator : public CGLinearOperator<double> {
	public:
		int XSize, YSize, GridSize;

		char *region;

	public:
		PoissonDirichlet_CGLinearOperator( int _XSize, int _YSize, int _GridSize, char *_region ) : 
		XSize(_XSize), YSize(_YSize), GridSize(_GridSize), region(_region) { };

		void operator() ( const double *x, double *Ax ) {
			// assuming two components
			const double *xc=x;
			const double *yc=x+GridSize;
			double G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my;

			for (int x2=0, p=0; x2<YSize; x2++) {
				for (int x1=0; x1<XSize; x1++, p++) {
					if (!region[p]) {
						G1px = x1 < XSize - 1 && !region[p+1] ? xc[p+1] : 0;
						G2px = x1 < XSize - 1 && !region[p+1] ? yc[p+1] : 0;
						G1mx = x1 > 0         && !region[p-1] ? xc[p-1] : 0;
						G2mx = x1 > 0         && !region[p-1] ? yc[p-1] : 0;

						G1py = x2 < YSize - 1 && !region[p+XSize] ? xc[p+XSize] : 0;
						G2py = x2 < YSize - 1 && !region[p+XSize] ? yc[p+XSize] : 0;
						G1my = x2 > 0         && !region[p-XSize] ? xc[p-XSize] : 0;
						G2my = x2 > 0         && !region[p-XSize] ? yc[p-XSize] : 0;

						Ax[p]          = -(G1px + G1mx + G1py + G1my) + 4*xc[p];
						Ax[p+GridSize] = -(G2px + G2mx + G2py + G2my) + 4*yc[p];
					}
				}
			}

			return;
		}
};

#endif