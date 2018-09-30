#include "mex.h"
#include "Motion_Sobolev.h"
#include <stdio.h>
#include <string.h>

/*
 * [fwarp, bwarp, region_warp, residual] 
 *   = mexMotionEstimation( I0, I1, region, dt, iters, levels, steps_v, steps_rel_u, steps_rel_d, occlusion_threshold, fwarpi, bwarpi )
 *
 *  I0, I1 - color images (3 components); double
 *  region - binary indicator for region (uint8)
 *  dt     - step size for warp evolution (not more than 0.5)
 *  iters  - number of iterations for warp evolution
 *  levels - number of levels in pyramid
 *  steps_v - # vcycles in multigrid
 *  steps_rel_u - relaxation steps in multigrid (upward direction)
 *  steps_rel_d - relaxation steps in multigrid (downward direction)
 *  occlusion_threshold - threshold for occlusion detection
 *  beta - smoothness (to avoid CC in MG implementation)
 *  fwarpi - warp initialization
 *  bwarpi - warp initialization
 * 
 *  fwarp  - forward warp displacement (defined on whole image domain): w(x) - x
 *  bwarp  - backward warp (defined on whole image domain) : w^{-1}(x)
 *  region_warp - warped region binary mask
 *  residual - | I1(w(x)) - I0(x) |^2 defined on original region
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int *size;
	ARRAYTYPE *I0, *I1;
	char *region;
	int XSize, YSize, GridSize, levels, iters;
	double occlusion_threshold, dt;
	Motion_Sobolev ms;
	ARRAYTYPE *fwarpi, *bwarpi;
	double beta = 1e-4;
	int steps_v, steps_rel_u, steps_rel_d;

	size     = mxGetDimensions( prhs[0] );
	XSize    = size[0];
	YSize    = size[1];
	GridSize = XSize*YSize;

	I0     = (ARRAYTYPE*)mxGetData(prhs[0]);
	I1     = (ARRAYTYPE*)mxGetData(prhs[1]);
	region =      (char*)mxGetData(prhs[2]);

	dt                  = (double)mxGetScalar(prhs[3]);
	iters               =    (int)mxGetScalar(prhs[4]);
	levels              =    (int)mxGetScalar(prhs[5]);
	steps_v             =    (int)mxGetScalar(prhs[6]);
	steps_rel_u         =    (int)mxGetScalar(prhs[7]);
	steps_rel_d         =    (int)mxGetScalar(prhs[8]);
	occlusion_threshold = (double)mxGetScalar(prhs[9]);
	beta                = (double)mxGetScalar(prhs[10]);

	double area=0;
	for (int p=0; p<XSize*YSize; p++)
		area += region[p]==1 ? 1 : 0;
	area /= 640*480; // scaling to make it scale invariant

	beta /= area*area;

	if ( !ms.allocate( XSize, YSize, levels, beta ) ) {
		mexPrintf("Allocation failure\n");
	}

#ifdef DEBUG
	printf( "XSize=%d, YSize=%d, dt=%f, iters=%d, levels=%d, occlusion_threshold=%f\n", 
		XSize, YSize, dt, iters, levels, occlusion_threshold );
	printf( "steps_v = %d, steps_rel_u = %d, steps_rel_d = %d\n", steps_v, steps_rel_u, steps_rel_d );

	for (int l=0; l<ms.mg.levels; l++) {
		printf("beta[%d] = %f\n", l, ms.mg.beta[l]);
	}
#endif

	ms.occlusion_threshold = occlusion_threshold;
	ms.mg.steps0 = steps_v;
	ms.mg.steps1 = steps_rel_u;
	ms.mg.steps2 = steps_rel_d;

	for (int f=0; f<3; f++) {
		memcpy( ms.I0[0][f], I0 + f*ms.GridSize[0], sizeof(*I0) * ms.GridSize[0] );
		memcpy( ms.I1[0][f], I1 + f*ms.GridSize[0], sizeof(*I1) * ms.GridSize[0] );
	}

	memcpy( ms.region[0], region, sizeof(*region)*ms.GridSize[0] );

	if (nrhs > 11) {
		fwarpi = (ARRAYTYPE*)mxGetData( prhs[11] );
		bwarpi = (ARRAYTYPE*)mxGetData( prhs[12] );

		ms.initializeWarps( 0, fwarpi, bwarpi );
	}
	else {
		ms.initializeWarps( 0 );
	}

	ms.evolveForwardBackwardWarps( 0, iters, dt );
	ms.computeResidual( 0 );

	mwSize dims[3];

	dims[0]=XSize; dims[1]=YSize; dims[2]=2;
	plhs[0] = mxCreateNumericArray( (mwSize)3, dims, mxDOUBLE_CLASS, mxREAL );
	plhs[1] = mxCreateNumericArray( (mwSize)3, dims, mxDOUBLE_CLASS, mxREAL );
	plhs[2] = mxCreateNumericArray( (mwSize)2, dims, mxINT8_CLASS, mxREAL );
	plhs[3] = mxCreateNumericArray( (mwSize)2, dims, mxDOUBLE_CLASS, mxREAL );

	double *fwarp       = (double*)mxGetData( plhs[0] );
	double *bwarp       = (double*)mxGetData( plhs[1] );
	char   *region_warp = (char*)mxGetData( plhs[2] );
	double *residual    = (double*)mxGetData( plhs[3] );

	int x, y;

	for (int p=0; p<ms.GridSize[0]; p++) {
		x = p % ms.XSize[0]; y = p / ms.XSize[0];
		//F[p]=ms.L2grad[p];
		//F[p+ms.GridSize[0]]=ms.L2grad[p+ms.GridSize[0]];
		//G[p]=ms.SobolevGrad[p];
		//G[p+ms.GridSize[0]]=ms.SobolevGrad[p+ms.GridSize[0]];
		region_warp[p]=ms.region_warp[0][p];

		// convert to displacements
		fwarp[p] = ms.forward_warp[0][0][p] - x;
		fwarp[p+ms.GridSize[0]] = ms.forward_warp[0][1][p] - y;
		bwarp[p] = ms.backward_warp[0][0][p] - x;
		bwarp[p+ms.GridSize[0]] = ms.backward_warp[0][1][p] - y;

		residual[p] = ms.residual[p];
	}

	ms.deallocate();

	return;
}