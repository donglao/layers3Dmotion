#include "mex.h"
#include "Smoothing.h"
#include <string.h>

void computePhiWarp( double *phi, double *bwarp, double *phi_warp, double *occl_map, double sigma_phi, int XSize, int YSize );

/*
 * label = mexComputePropagatedSegmentation( phi, bwarp, occlusion_map, sigma_phi )
 * 
 *   phi   - M x N x Nregions double matrix
 *   bwarp - M x N x 2 x Nregions double matirx (backward warp)
 *   sigma_phi - smoothing of sigma of phi's
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int *size;
	int Nregions;
	double *phi, *bwarp, *phi_warped, *occlusion_map;
	int XSize, YSize, GridSize;
	double sigma_phi;

	size=mxGetDimensions( prhs[1] );
	XSize=size[0]; YSize=size[1]; 
	GridSize = XSize*YSize;
	Nregions = size[3];

	phi = (double*)mxGetData( prhs[0] );
	bwarp = (double*)mxGetData( prhs[1] );
	occlusion_map = (double*)mxGetData( prhs[2] );
	sigma_phi = (double)mxGetScalar( prhs[3] );

	if ( ! ( phi_warped = new double[ GridSize * Nregions] ) ) {
		mexPrintf("Allocation error\n"); 
		return;
	}

	for (int i=0; i<Nregions; i++) {
		computePhiWarp( phi + i*GridSize, bwarp + i * 2 * GridSize, phi_warped + i*GridSize, occlusion_map + i*GridSize, sigma_phi, XSize, YSize );
	}

	mwSize dims[3];

	dims[0]=XSize; dims[1]=YSize; dims[2]=Nregions;
	plhs[0] = mxCreateNumericArray( (mwSize)2, dims, mxINT8_CLASS, mxREAL );
	plhs[1] = mxCreateNumericArray( (mwSize)3, dims, mxDOUBLE_CLASS, mxREAL );

	char *label = (char*) mxGetData( plhs[0] );
	double *phi_warp_out = (double*) mxGetData( plhs[1] );

	for (int p=0; p<GridSize; p++) {
		double max_phi_w = phi_warped[p];
		int max_label=0;
		for (int r=1; r<Nregions; r++) {
			if ( max_phi_w < phi_warped[p + r*GridSize] ) {
				max_phi_w = phi_warped[p + r*GridSize];
				max_label = r;
			}
		}
		label[p] = max_phi_w >= 0.5 ? max_label : -1;
	}

	for (int p=0; p<GridSize*Nregions; p++)
		phi_warp_out[p]=phi_warped[p];

	delete[] phi_warped;

	return;
}


void computePhiWarp( double *phi, double *bwarp, double *phi_warp, double *occl_map, double sigma_phi, int XSize, int YSize )
{
	int GridSize = XSize * YSize;
	double *bwarp_x = bwarp, *bwarp_y = bwarp + GridSize;

	for (int p=0; p<GridSize; p++) {
		int x = p % XSize, y= p /XSize;
		double xbw = bwarp_x[p]+x > XSize-1 ? XSize-1 : ( bwarp_x[p]+x < 0 ? 0 : bwarp_x[p]+x );
		double ybw = bwarp_y[p]+y > YSize-1 ? YSize-1 : ( bwarp_y[p]=y < 0 ? 0 : bwarp_y[p]+y );

		if ( !( xbw < 0 || ybw < 0 || xbw > XSize-1 || ybw > YSize - 1 ) ) {
			int xbwi = int(xbw);
			int ybwi = int(ybw);
			int pw = xbwi + ybwi * XSize;
			double dx = xbw-xbwi;
			double dy = ybw-ybwi;
			double I0pX, I0pY, I0pXpY, occpX, occpY, occpXpY, occ_interp;

			I0pX   = xbwi < XSize - 1                        ? phi[pw+1]       : phi[pw];
			I0pY   = ybwi < YSize - 1                        ? phi[pw+XSize]   : phi[pw];
			I0pXpY = xbwi < XSize - 1 && ybwi < YSize - 1    ? phi[pw+1+XSize] : phi[pw];

			occpX   = xbwi < XSize - 1                        ? occl_map[pw+1]       : occl_map[pw];
			occpY   = ybwi < YSize - 1                        ? occl_map[pw+XSize]   : occl_map[pw];
			occpXpY = xbwi < XSize - 1 && ybwi < YSize - 1    ? occl_map[pw+1+XSize] : occl_map[pw];

			occ_interp = (1-dx)  * (1-dy) * occl_map[pw]      +
                    		dx   * (1-dy) * occpX             +
                    	 (1-dx)  *    dy  * occpY             +
			                dx   *    dy  * occpXpY;

			phi_warp[p] = (1-dx) * (1-dy) * phi[pw]      +
			                 dx  * (1-dy) * I0pX         +
			              (1-dx) *    dy  * I0pY         +
			                 dx  *    dy  * I0pXpY;

       		//phi_warp[p] -= occ_interp;
			//phi_warp[p] = occ_interp > 0.5 ?  0 : phi_warp[p];
		}
		else {
			phi_warp[p] = 0;
		}
	}

	double *phi_warp_smooth = new double[GridSize];

//	CSmoothing::Blur( phi_warp_smooth, phi_warp, XSize, YSize, 1, sigma_phi );
//	memcpy( phi_warp, phi_warp_smooth, GridSize * sizeof(*phi_warp) );

	delete[] phi_warp_smooth;

	return;
}