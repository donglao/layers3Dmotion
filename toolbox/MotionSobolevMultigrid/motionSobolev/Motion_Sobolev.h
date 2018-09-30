#ifndef MOTION_SOBOLEV
#define MOTION_SOBOLEV

#include <math.h>
#include <stdio.h>
#include "vector2D.h"
#include "PoissonMultigrid.h"
#include "ConnectedComponentDetection.h"
#include "params.h"
#include <float.h>
#include <string.h>
#include <omp.h>

typedef double ARRAYTYPE;

const static int FSize=3;


class Motion_Sobolev {

	public:
		int *XSize, *YSize, *GridSize;
		int levels;

		ARRAYTYPE ***I0, ***I1; // images, pyramid representation

		char **region;      // region indicator functions, pyramid
		char **region_warp; // array for forward warped region

		ARRAYTYPE ***forward_warp, ***backward_warp; // warps, pyramid

		double occlusion_threshold;

		ARRAYTYPE *L2grad, *SobolevGrad;

		// for conjugateGradient
		PoissonMultiGrid< vector2D<double> > mg;

		// for determining connected components
		ConnectedComponentDetection ccd;

		// for gradients in updating warps
		ARRAYTYPE *dfwx, *dfwy, *dbwx, *dbwy;

		ARRAYTYPE *residual, *residual_vec;

		unsigned char *edgebits;   //Bits describing neighborhoods of grid points
		enum {XPOS=1,XNEG=2,YPOS=4,YNEG=8,XEDGES=3,YEDGES=12,EDGES=15};

//		int Nblocks; /// for parallelization

	public:
		Motion_Sobolev() {
			XSize=YSize=GridSize=0;
			I0=I1=0;
			region=region_warp=0;
			forward_warp=backward_warp=0;
			L2grad=SobolevGrad=0;
			dfwx=dfwy=0; dbwx=dbwy=0;
			residual=0;
			edgebits=0;
//			Nblocks=1;
			residual_vec=0;
		};

		int allocate(int XSize0, int YSize0, int levels, double beta);
		void deallocate();

		void evolveForwardBackwardWarps(int l, int iters, double dt);
		void computeResidual(int l);
		void computeResidual2(int l);

		void initializeWarps(int l, ARRAYTYPE *fwarpi, ARRAYTYPE *bwarpi);
		void initializeWarps(int l);

	private:
		void markEdges();
		void createPyramid();
		void computeL2Gradient(int l);
		void computeGradTranslation(int l, ARRAYTYPE* trans);
		double computeSobolevDeformation(int l, ARRAYTYPE* trans, int usePrevInit);
//		double extendSobolevDeformation(int l, int usePrevInit);
		void MeanNormalizeL2Gradient(int l, ARRAYTYPE* trans);
		double computeEnergy(int l);
		void updateForwardMap( int l, double dt_scale );
		void updateBackwardMap( int l, int deforming, double dt, double& dt_scale, ARRAYTYPE* trans );

		// I1 - I0(w^{-1}(x)), x in w( region )
		void diff_I1_I0BackWarp( int l, int p, ARRAYTYPE *Iwarp ) {
			int XSize_l = XSize[l];
			int YSize_l = YSize[l];
			ARRAYTYPE **I0_l = I0[l];
			ARRAYTYPE **I1_l = I1[l];
			double xbw = backward_warp[l][0][p] > XSize_l-1 ? XSize_l-1 : ( backward_warp[l][0][p] < 0 ? 0 : backward_warp[l][0][p] );
			double ybw = backward_warp[l][1][p] > YSize_l-1 ? YSize_l-1 : ( backward_warp[l][1][p] < 0 ? 0 : backward_warp[l][1][p] );
			int xbwi = int(xbw);
			int ybwi = int(ybw);
			int pw = xbwi + ybwi * XSize_l;
			double dx = xbw-xbwi;
			double dy = ybw-ybwi;
			double I0pX, I0pY, I0pXpY;

			for (int f=0; f<FSize; f++) {
				I0pX   = edgebits[pw]&XPOS         ? I0_l[f][pw] : I0_l[f][pw+1];
				I0pY   = edgebits[pw]&YPOS         ? I0_l[f][pw] : I0_l[f][pw+XSize_l];
				I0pXpY = edgebits[pw]&(XPOS|YPOS)  ? I0_l[f][pw] : I0_l[f][pw+1+XSize_l];

				Iwarp[f] = (1-dx) * (1-dy) * I0_l[f][pw]  +
				              dx  * (1-dy) * I0pX         +
				           (1-dx) *    dy  * I0pY         +
				              dx  *    dy  * I0pXpY;
				Iwarp[f]=I1_l[f][p]-Iwarp[f];
			}

			return;
		};

		// I1(w(x)) - I0(x), x in region
		void diff_I1warp_I0( int l, int p, ARRAYTYPE *Iwarp ) {
			int XSize_l = XSize[l];
			int YSize_l = YSize[l];
			ARRAYTYPE **I0_l = I0[l];
			ARRAYTYPE **I1_l = I1[l];
			double xfw = forward_warp[l][0][p] > XSize_l-1 ? XSize_l-1 : (forward_warp[l][0][p] < 0 ? 0 : forward_warp[l][0][p]);
			double yfw = forward_warp[l][1][p] > YSize_l-1 ? YSize_l-1 : (forward_warp[l][1][p] < 0 ? 0 : forward_warp[l][1][p]);
			int xfwi = int(xfw);
			int yfwi = int(yfw);
			int pw = xfwi + yfwi * XSize_l;
			double dx = xfw-xfwi;
			double dy = yfw-yfwi;
			double I1pX, I1pY, I1pXpY;

			for (int f=0; f<FSize; f++) {
				I1pX   = edgebits[pw]&XPOS         ? I1_l[f][pw] : I1_l[f][pw+1];
				I1pY   = edgebits[pw]&YPOS         ? I1_l[f][pw] : I1_l[f][pw+XSize_l];
				I1pXpY = edgebits[pw]&(XPOS|YPOS)  ? I1_l[f][pw] : I1_l[f][pw+1+XSize_l];

				Iwarp[f] = (1-dx) * (1-dy) * I1_l[f][pw]  +
				              dx  * (1-dy) * I1pX         +
				           (1-dx) *    dy  * I1pY         +
				              dx  *    dy  * I1pXpY;
				Iwarp[f]=Iwarp[f] - I0_l[f][p];
			}

			return;
		};


		void computeImageGradient(int l, int f, ARRAYTYPE *gradI, int p ) {
			int XSize_l = XSize[l];
			int YSize_l = YSize[l];
			ARRAYTYPE *I1_l_f = I1[l][f];
			double I1pX, I1mX, I1pY, I1mY;

			I1pX = edgebits[p]&XPOS ? I1_l_f[p] : I1_l_f[p+1];
			I1mX = edgebits[p]&XNEG ? I1_l_f[p] : I1_l_f[p-1];
			I1pY = edgebits[p]&YPOS ? I1_l_f[p] : I1_l_f[p+XSize_l];
			I1mY = edgebits[p]&YNEG ? I1_l_f[p] : I1_l_f[p-XSize_l];

			gradI[0] = (I1pX-I1mX)/2;
			gradI[1] = (I1pY-I1mY)/2;

			return;
		}

		// compute w(R), i.e., forward warped region
		void computeWarpedRegion( int l );

		// compute det ( nabla w^{-1}(x) )
		double computeDetJacobianWarpInv(int l, int p) {
			int XSize_l = XSize[l];
			int YSize_l = YSize[l];
			char *region_warp_l = region_warp[l];
			ARRAYTYPE **bwarp = backward_warp[l];
			double J[2][2]; // Jacobian

			int pp_in, pm_in, incr;
			ARRAYTYPE wp, wm;

			for (int c=0; c<2; c++) { // warp component
				for ( int s=0; s<2; s++ ) { // x,y spatial locations
					pp_in = s==0 ? !(edgebits[p]&XPOS) : !(edgebits[p]&YPOS);
					pm_in = s==0 ? !(edgebits[p]&XNEG) : !(edgebits[p]&YNEG);
					incr  = s==0 ? 1 : XSize_l;

					wp = pp_in ? ( region_warp_l[p+incr] ? bwarp[c][p+incr] : bwarp[c][p] ) : bwarp[c][p];
					wm = pm_in ? ( region_warp_l[p-incr] ? bwarp[c][p-incr] : bwarp[c][p] ) : bwarp[c][p];
					J[c][s] = (wp-wm)/2;
				}
			}

			return J[0][0]*J[1][1] - J[0][1]*J[1][0];  // determinent
		}


		// small helper math functions

		double vnormsq( ARRAYTYPE *v ) {
			double norm = 0;

			for (int f=0; f<FSize; f++) norm += v[f]*v[f];

			return norm;
		}

		void add( ARRAYTYPE *sum, ARRAYTYPE *a, double scalar ) {
			for (int f=0; f<FSize; f++) sum[f]+=a[f]*scalar;
			return;
		}

		ARRAYTYPE maxabs( ARRAYTYPE *a, int N )
		{
			ARRAYTYPE max = a[0];

			for (int i=1; i<N; i++) {
				max = fabs(a[i]) > max ? fabs(a[i]) : max;
			}

			return max;
		}

};

#endif
