#ifndef POISSONDIRICHLET_CGLINEAROPERATOR
#define POISSONDIRICHLET_CGLINEAROPERATOR

#include "CGLinearOperator.h"

#ifdef PARALLELIZE_CODE
#include <omp.h>
#endif

class PoissonDirichlet_CGLinearOperator : public CGLinearOperator<double> {
	public:
		int XSize, YSize, GridSize;

		char *region;
		unsigned char *edgebits;
		int Nblocks;

	public:
		PoissonDirichlet_CGLinearOperator( int _XSize, int _YSize, int _GridSize, char *_region, unsigned char *_edgebits, int _Nblocks ) : 
		XSize(_XSize), YSize(_YSize), GridSize(_GridSize), region(_region), edgebits(_edgebits), Nblocks(_Nblocks) { };

		void operator() ( const double *x, double *Ax ) {
			// assuming two components
			const double *xc=x;
			const double *yc=x+GridSize;
			double *Ax_GridSize = Ax + GridSize;
			double G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my;

#ifdef PARALLELIZE_CODE
			#pragma omp parallel for private(G1px, G2px, G1mx, G2mx, G1py, G2py, G1my, G2my)
#endif

			for (BLOCKS(block)) {
				BLOCK_ENDS(GridSize);

				for (PIXELS(p)) {
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

#endif