#include "Motion_Sobolev.h"
#include <stdio.h>
#include <time.h>


int Motion_Sobolev::allocate(int XSize0, int YSize0, int levels, double beta)
{

#define ALLOCATE( array, arraytype, Nelements )       \
	if ( !( array = new arraytype[Nelements] ) ) {    \
		deallocate();                                 \
		return 0;                                     \
	}

	this->levels=levels;

	ALLOCATE( XSize,    int, levels );
	ALLOCATE( YSize,    int, levels );
	ALLOCATE( GridSize, int, levels );

	XSize[0]=XSize0; YSize[0]=YSize0; GridSize[0]=XSize0*YSize0;

	for (int l=1; l<levels; l++) {
		XSize[l]   = XSize[l-1]/2;
		YSize[l]   = YSize[l-1]/2;
		GridSize[l]= XSize[l]*YSize[l];
	}

	ALLOCATE( I0,            ARRAYTYPE**, levels );
	ALLOCATE( I1,            ARRAYTYPE**, levels );
	ALLOCATE( region,        char*,       levels );
	ALLOCATE( region_warp,   char*,       levels );
	ALLOCATE( forward_warp,  ARRAYTYPE**, levels );
	ALLOCATE( backward_warp, ARRAYTYPE**, levels );

	for	(int l=0; l<levels; l++) {
		ALLOCATE( I0[l],          ARRAYTYPE*, FSize );
		ALLOCATE( I1[l],          ARRAYTYPE*, FSize );
		ALLOCATE( region[l],      char,       GridSize[l] );
		ALLOCATE( region_warp[l], char,       GridSize[l] );

		for (int f=0; f<FSize; f++) {
			ALLOCATE( I0[l][f], ARRAYTYPE, GridSize[l] );
			ALLOCATE( I1[l][f], ARRAYTYPE, GridSize[l] );
		}

		ALLOCATE( forward_warp[l],  ARRAYTYPE*, 2 );
		ALLOCATE( backward_warp[l], ARRAYTYPE*, 2 );

		for (int c=0; c<2; c++) {
			ALLOCATE( forward_warp[l][c],  ARRAYTYPE, GridSize[l] );
			ALLOCATE( backward_warp[l][c], ARRAYTYPE, GridSize[l] );
		}
	}

	ALLOCATE( L2grad,      ARRAYTYPE, 2*GridSize[0] );
	ALLOCATE( SobolevGrad, ARRAYTYPE, 2*GridSize[0] );


#ifdef PARALLELIZE_CODE
//	Nblocks = omp_get_num_procs();
//	Nblocks=16;
	//	omp_set_num_threads(5);
#endif

	//printf("Nblocks=%d\n", Nblocks);
	int sizes[2];
	sizes[0]=XSize0; sizes[1]=YSize0;

	if (! mg.allocate( 2, sizes, levels, beta ) )
		{ deallocate(); return 0; }

	mg.steps0 = 2;
	mg.steps1 = mg.steps2 = 1;

	if ( !ccd.allocate( XSize[0], YSize[0]) )
		{ deallocate(); return 0; }

	ALLOCATE( dfwx, ARRAYTYPE, GridSize[0] );
	ALLOCATE( dfwy, ARRAYTYPE, GridSize[0] );

	ALLOCATE( dbwx, ARRAYTYPE, GridSize[0] );
	ALLOCATE( dbwy, ARRAYTYPE, GridSize[0] );

	ALLOCATE( residual, ARRAYTYPE, GridSize[0] );
	ALLOCATE( residual_vec, ARRAYTYPE, FSize*GridSize[0] );

	ALLOCATE( edgebits, unsigned char, GridSize[0] );

	markEdges();

#undef ALLOCATE

	return 1;
}


void Motion_Sobolev::deallocate()
{

#define DEALLOCATE( array )  \
	if ( array!=0 ) { delete[] array; array=0; }

	DEALLOCATE( XSize );
	DEALLOCATE( YSize );
	DEALLOCATE( GridSize );

	if ( I0 != 0 ) {
		for (int l=0; l<levels; l++) {
			if ( I0[l]!=0 ) {
				for (int f=0; f<FSize; f++) {
					DEALLOCATE( I0[l][f] );
				}
				DEALLOCATE( I0[l] );
			}
		}
		DEALLOCATE( I0 );
	}

	if ( I1 != 0 ) {
		for (int l=0; l<levels; l++) {
			if ( I1[l]!=0 ) {
				for (int f=0; f<FSize; f++) {
					DEALLOCATE( I1[l][f] );
				}
				DEALLOCATE( I1[l] );
			}
		}
		DEALLOCATE( I1 );
	}

	if ( region != 0 ) {
		for (int l=0; l<levels; l++) {
			DEALLOCATE( region[l] );
		}
		DEALLOCATE( region );
	}

	if ( region_warp != 0 ) {
		for (int l=0; l<levels; l++) {
			DEALLOCATE( region_warp[l] );
		}
		DEALLOCATE( region_warp );
	}

	if ( forward_warp !=0  ) {
		for (int l=0; l<levels; l++) {
			if ( forward_warp[l]!=0 ) {
				DEALLOCATE( forward_warp[l][0] );
				DEALLOCATE( forward_warp[l][1] );
			}
			DEALLOCATE( forward_warp[l] );
		}
		DEALLOCATE( forward_warp );
	}

	if ( backward_warp !=0  ) {
		for (int l=0; l<levels; l++) {
			if ( backward_warp[l]!=0 ) {
				DEALLOCATE( backward_warp[l][0] );
				DEALLOCATE( backward_warp[l][1] );
			}
			DEALLOCATE( backward_warp[l] );
		}
		DEALLOCATE( backward_warp );
	}

	DEALLOCATE( L2grad );
	DEALLOCATE( SobolevGrad );

	mg.deallocate();
	ccd.deallocate();

	DEALLOCATE( dfwx );
	DEALLOCATE( dfwy );

	DEALLOCATE( dbwx );
	DEALLOCATE( dbwy );

	DEALLOCATE( residual );
	DEALLOCATE( edgebits );
	DEALLOCATE( residual_vec );

#undef DEALLOCATE

	return;
}


void Motion_Sobolev::markEdges() {
  int p;
  int Y=XSize[0];

  memset(edgebits,0,GridSize[0]);

  for (p=0;             p<Y;           p++)   edgebits[p]|=YNEG;
  for (p=GridSize[0]-Y; p<GridSize[0]; p++)   edgebits[p]|=YPOS;
  for (p=0;             p<GridSize[0]; p+=Y)  edgebits[p]|=XNEG;
  for (p=Y-1;           p<GridSize[0]; p+=Y)  edgebits[p]|=XPOS;
}


void Motion_Sobolev::initializeWarps(int l) {
	ARRAYTYPE **bwarp=backward_warp[l];
	ARRAYTYPE **fwarp=forward_warp[l];
	int XSize_l = XSize[l];
	int YSize_l = YSize[l];

	for (int y=0, p=0; y<YSize_l; y++) {
		for (int x=0; x<XSize_l; x++, p++) {
			bwarp[0][p]=x; bwarp[1][p]=y;
			fwarp[0][p]=x; fwarp[1][p]=y;
		}
	}

	return;
}


void Motion_Sobolev::computeWarpedRegion( int l )
{
	ARRAYTYPE **bwarp   = backward_warp[l];
	char *region_l      = region[l];
	char *region_warp_l = region_warp[l];
	int XSize_l         = XSize[l];
	int YSize_l         = YSize[l];
	int GridSize_l      = GridSize[l];
	int xw, yw, xwi, ywi;
	double dxw, dyw;

#define INBOUNDS( x, y )      ( (x)>=0 && (x)<XSize_l && (y)>=0 && (y)<YSize_l )
#define INSIDE_REGION( x, y ) ( region_l[(x)+(y)*XSize_l] )

	for (int p=0; p<GridSize_l; p++) {
		xw = bwarp[0][p]; xwi = int(xw);
		yw = bwarp[1][p]; ywi = int(yw);
		
		dxw = xw - xwi;
		dyw = yw - ywi;

		// check whether any pixel in box of xw,yw is in region
		if ( ( INBOUNDS( xw,   yw   ) && INSIDE_REGION( xw,   yw ) )   ||
			 ( INBOUNDS( xw+1, yw   ) && INSIDE_REGION( xw+1, yw ) && fabs(dxw)>EPS ) ||
			 ( INBOUNDS( xw,   yw+1 ) && INSIDE_REGION( xw,   yw+1 ) && fabs(dyw)>EPS ) ||
			 ( INBOUNDS( xw+1, yw+1 ) && INSIDE_REGION( xw+1, yw+1 ) && fabs(dxw)>EPS && fabs(dyw)>EPS ) )
			region_warp_l[p] = 1;
		else
			region_warp_l[p] = 0;
	}

	return;
}


void Motion_Sobolev::computeL2Gradient(int l)
{
	ARRAYTYPE **bwarp   = backward_warp[l];
	char *region_warp_l = region_warp[l];
	int XSize_l    = XSize[l];
	int YSize_l    = YSize[l];
	int GridSize_l = GridSize[l];
	ARRAYTYPE diff[FSize];
	double diff_norm_sq;
	ARRAYTYPE gradIf[2];
	double detwinv;
	ARRAYTYPE *L2grad_x = L2grad, *L2grad_y = L2grad + GridSize_l;

	computeWarpedRegion( l );

#ifdef PARALLELIZE_CODE
	#pragma omp parallel for private( detwinv, diff_norm_sq, diff, gradIf )
#endif

	for (BLOCKS(block)) {
		BLOCK_ENDS(GridSize_l);

		for (PIXELS(p)) {
			if ( region_warp_l[p] ) {

		        // I1( p ) - I0( w^{-1}(p) )
				diff_I1_I0BackWarp( l, p, diff );
				diff_norm_sq = vnormsq( diff );

		        // compute detw^{-1}
				detwinv = computeDetJacobianWarpInv( l, p );

		        //printf("diff_norm_sq = %f, detwinv=%f\n", diff_norm_sq, detwinv);

				L2grad_x[p]=L2grad_y[p]=0;

				if ( diff_norm_sq < occlusion_threshold ) {
					for (int f=0; f<FSize; f++) {

						// grad I1[f][p]
						computeImageGradient( l, f, gradIf, p );
						// grad I1[f][p] * ( I1( p ) - I0( w^{-1}(p) ) ) * detwinv
						L2grad_x[p] += gradIf[0]*diff[f]*detwinv;
						L2grad_y[p] += gradIf[1]*diff[f]*detwinv;
					}
				}
			}
			else {
				L2grad_x[p]=L2grad_y[p]=0;
			}
		}
	}

	return;
}


void Motion_Sobolev::computeGradTranslation(int l, ARRAYTYPE* trans)
{
	int GridSize_l = GridSize[l];
	char *region_warp_l = region_warp[l];
	int count=0;
	ARRAYTYPE *L2grad_x = L2grad, *L2grad_y = L2grad + GridSize_l;

	trans[0]=trans[1]=0;

	for (int p=0; p<GridSize_l; p++) {
		if (region_warp_l[p]) {
			trans[0]+=L2grad_x[p];
			trans[1]+=L2grad_y[p];
			count++;
		}
	}

	if ( count!=0 ) { trans[0]/=count; trans[1]/=count; }

	return;
}



void Motion_Sobolev::MeanNormalizeL2Gradient(int l, ARRAYTYPE* trans)
{
	int GridSize_l      = mg.GridSize[l];
	int XSize_l         = mg.sizes[l][0];
	int YSize_l         = mg.sizes[l][1];
	char *region_warp_l = new char[GridSize_l];
	vector2D<double> *L2grad_l = mg.b[l];
	int area_cc, area;
	ARRAYTYPE trans_cc[2];

	// setup and find connected components
	ccd.XSize = XSize_l;
	ccd.YSize = YSize_l;
	ccd.GridSize = GridSize_l;
	for (int p=0; p<GridSize_l; p++)
		region_warp_l[p] = (char)(mg.region_prob[l][p]>=0.5);
	ccd.findConnectedComponents( region_warp_l );

	trans[0]=trans[1]=0;
	area=0;

	for (int c=1; c<=ccd.Ncomponents; c++) {
		trans_cc[0]=trans_cc[1]=0;
		area_cc=0;

		for (int p=0; p<GridSize_l; p++) {
			if (ccd.ccLabel[p]==c) {
				trans_cc[0]+=L2grad_l[p].x;
				trans_cc[1]+=L2grad_l[p].y;
				area_cc++;
			}
		}

		trans[0]+=trans_cc[0]; trans[1]+=trans_cc[1];
		area+=area_cc;

		trans_cc[0]/=area_cc; trans_cc[1]/=area_cc;

	//		printf( "trans_cc=(%f,%f), area_cc=%d\n", trans_cc[0], trans_cc[1], area_cc );

		for (int p=0; p<GridSize_l; p++) {
			if (ccd.ccLabel[p]==c) {
				L2grad_l[p].x-=trans_cc[0];
				L2grad_l[p].y-=trans_cc[1];
			}
		}
	}

	trans[0]/=area; trans[1]/=area;

	//	printf( "Ncc=%d, trans=(%f,%f)\n", ccd.Ncomponents, trans[0], trans[1] );

	delete[] region_warp_l;

	return;
}


double Motion_Sobolev::computeSobolevDeformation( int ll, ARRAYTYPE* trans, int usePrevInit )
{
	int GridSize_l = GridSize[ll];
	vector2D<double> *b0 = mg.b[0];
	ARRAYTYPE trans_tmp[2];
	ARRAYTYPE *L2grad_x = L2grad;
	ARRAYTYPE *L2grad_y = L2grad + GridSize_l;
	char *region_warp_l = region_warp[ll];

	for (int p=0; p<GridSize_l; p++) {
		b0[p].x=L2grad_x[p];
		b0[p].y=L2grad_y[p];
		mg.region_prob[0][p] = region_warp_l[p];
	}
	//	MeanNormalizeL2Gradient( 0, trans );
	//	printf( "level = 0, trans = (%f, %f)\n", trans[0], trans[1] );

	for (int l=1; l<mg.levels; l++) {
		mg.restrict_data( mg.b[l-1], mg.b[l], l-1 );
		mg.restrict_data( mg.region_prob[l-1], mg.region_prob[l], l-1 );
//		mg.plotImage( mg.region_prob[l], l );
//		MeanNormalizeL2Gradient( l, trans_tmp );
//		printf( "level = %d, trans = (%f, %f)\n", l, trans_tmp[0], trans_tmp[1] );
	}

	double residual_start, residual_end;

	residual_start = mg.computeError( 0 );

	if (usePrevInit) {
		for (int p=0; p<GridSize_l; p++) {
			mg.x[0][p].x = SobolevGrad[p];
			mg.x[0][p].y = SobolevGrad[p+GridSize_l];
		}
		mg.solveUsingVCycle( 0, 1e-3, 1 );
#ifdef DEBUG
		printf("V-cycle:\n");
#endif
	}
	else {
#ifdef DEBUG
		printf("Full Multigrid:\n");
#endif
		arrayOps::set( mg.x[0], 0.0, GridSize_l );
//		printf("residual = %f\n", mg.computeError( 0 ) );
//		mg.plotImage( mg.residual, 0 );
//		mg.plotImage( mg.b[0], 0 );

		mg.full_multigrid( 1e-3 );
//		printf("residual = %f\n", mg.computeError( 0 ) );
//		mg.plotImage( mg.x[0], 0, "Full Multigrid" );
	}

	residual_end = mg.computeError( 0 );

#ifdef DEBUG
	printf("residual_start = %f, residual_end = %f\n", residual_start, residual_end);
#endif

	for (int p=0; p<GridSize_l; p++) {
		SobolevGrad[p] = mg.x[0][p].x;
		SobolevGrad[p+GridSize_l] = mg.x[0][p].y;
	}

	return 0;
}

/*
double Motion_Sobolev::computeSobolevDeformation(int l, ARRAYTYPE* trans, int usePrevInit )
{
	int GridSize_l = GridSize[l];
	char *region_warp_l = region_warp[l];
	ARRAYTYPE *L2grad_x = L2grad, *L2grad_y = L2grad + GridSize_l;

	// compute F - avg(F)
	MeanNormalizeL2Gradient(l, trans);

	// setup CG
	cg.N = 2*GridSize_l;
	memcpy( cg.active, region_warp_l, sizeof(*cg.active)*GridSize_l );
	memcpy( cg.active + GridSize_l, region_warp_l, sizeof(*cg.active)*GridSize_l );

	memset( cg.active + 2*GridSize_l, 0, 2*(GridSize[0]-GridSize_l)*sizeof(*cg.active) );

	memcpy( cg.b, L2grad_x, GridSize_l*sizeof(*cg.b) );
	memcpy( cg.b+GridSize_l, L2grad_y, GridSize_l*sizeof(*cg.b) );
	
	if (usePrevInit) {
		memcpy( cg.x, SobolevGrad, 2*GridSize_l*sizeof(*cg.x) );
	}
	else {
		memset( cg.x, 0, 2*sizeof(*cg.x)*GridSize_l );
	}

	cg.Aptr = new Poisson_CGLinearOperator( XSize[l], YSize[l], GridSize[l], region_warp_l, edgebits, Nblocks );

#ifdef PARALLELIZE_CODE
	double begin, end, time_spent;

	begin = omp_get_wtime();
#else
	clock_t begin, end;
	double time_spent;

	begin = clock();
#endif

	// solve Poisson equation using CG
	cg.computeSolution( CG_ERR_TOL );

#ifdef PARALLELIZE_CODE
	end = omp_get_wtime();
	time_spent = end - begin;
#else
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
#endif

	// copy result to SobolevGrad
	ARRAYTYPE *SobGrad_x = SobolevGrad, *SobGrad_y = SobolevGrad + GridSize_l;

	for (int p=0; p<GridSize_l; p++) {
		if(region_warp_l[p]) {
			SobGrad_x[p]=cg.x[p];
			SobGrad_y[p]=cg.x[p+GridSize_l];
		}
	}

	delete (Poisson_CGLinearOperator*)cg.Aptr;

	return time_spent;
}
*/

/*
void Motion_Sobolev::MeanNormalizeL2Gradient(int l, ARRAYTYPE* trans)
{
	int GridSize_l = GridSize[l];
	int XSize_l = XSize[l];
	int YSize_l = YSize[l];
	char *region_warp_l = region_warp[l];
	ARRAYTYPE *L2grad_x = L2grad, *L2grad_y = L2grad + GridSize_l;
	int area_cc, area;
	ARRAYTYPE trans_cc[2];

	// setup and find connected components
	ccd.XSize = XSize_l;
	ccd.YSize = YSize_l;
	ccd.GridSize = GridSize_l;
	ccd.findConnectedComponents( region_warp_l );

	trans[0]=trans[1]=0;
	area=0;

	for (int c=1; c<=ccd.Ncomponents; c++) {
		trans_cc[0]=trans_cc[1]=0;
		area_cc=0;

		for (int p=0; p<GridSize_l; p++) {
			if (ccd.ccLabel[p]==c) {
				trans_cc[0]+=L2grad_x[p];
				trans_cc[1]+=L2grad_y[p];
				area_cc++;
			}
		}

		trans[0]+=trans_cc[0]; trans[1]+=trans_cc[1];
		area+=area_cc;

		trans_cc[0]/=area_cc; trans_cc[1]/=area_cc;

	//		printf( "trans_cc=(%f,%f), area_cc=%d\n", trans_cc[0], trans_cc[1], area_cc );

		for (int p=0; p<GridSize_l; p++) {
			if (ccd.ccLabel[p]==c) {
				L2grad_x[p]-=trans_cc[0];
				L2grad_y[p]-=trans_cc[1];
			}
		}
	}

	trans[0]/=area; trans[1]/=area;

	//	printf( "Ncc=%d, trans=(%f,%f)\n", ccd.Ncomponents, trans[0], trans[1] );

	return;
}
*/



void Motion_Sobolev::updateBackwardMap( int l, int deforming, double dt, double& dt_scale, ARRAYTYPE* trans )
{
	int GridSize_l = GridSize[l];
	int XSize_l = XSize[l];
	int YSize_l = YSize[l];
	ARRAYTYPE *bwarp_x = backward_warp[l][0];
	ARRAYTYPE *bwarp_y = backward_warp[l][1];
	ARRAYTYPE *Gx = SobolevGrad;
	ARRAYTYPE *Gy = SobolevGrad + GridSize_l;
	ARRAYTYPE maxG, maxGx, maxGy;
	double b1x, b2x, b1y, b2y;

#ifdef PARALLELIZE_CODE
		#pragma parallel for private(b1x, b2x, b1y, b2y)
#endif

// compute backward map differentials
//		for (int p=0; p<GridSize_l; p++) {

	for (BLOCKS(block)) {

		BLOCK_ENDS( GridSize_l );

		for (PIXELS(p)) {
			// upwinding scheme
			if ( -Gx[p]<0 ) {
				if (edgebits[p]&XPOS) {
					b1x = bwarp_x[p] - bwarp_x[p-1];
					b2x = bwarp_y[p] - bwarp_y[p-1];
				}
				else {
					b1x= bwarp_x[p+1] - bwarp_x[p];
					b2x= bwarp_y[p+1] - bwarp_y[p];
				}
			}
			else {
				if (edgebits[p]&XNEG) {
					b1x = bwarp_x[p+1] - bwarp_x[p];
					b2x = bwarp_y[p+1] - bwarp_y[p];
				}
				else {
					b1x= bwarp_x[p] - bwarp_x[p-1];
					b2x= bwarp_y[p] - bwarp_y[p-1];
				}
			}

			if ( -Gy[p]<0 ) {
				if (edgebits[p]&YPOS) {
					b1y = bwarp_x[p] - bwarp_x[p-XSize_l];
					b2y = bwarp_y[p] - bwarp_y[p-XSize_l];
				}
				else {
					b1y = bwarp_x[p+XSize_l] - bwarp_x[p];
					b2y = bwarp_y[p+XSize_l] - bwarp_y[p];
				}
			}
			else {
				if (edgebits[p]&YNEG) {
					b1y = bwarp_x[p+XSize_l] - bwarp_x[p];
					b2y = bwarp_y[p+XSize_l] - bwarp_y[p];
				}
				else {
					b1y= bwarp_x[p] - bwarp_x[p-XSize_l];
					b2y= bwarp_y[p] - bwarp_y[p-XSize_l];
				}
			}

			dbwx[p] = b1x * Gx[p] + b1y * Gy[p];
			dbwy[p] = b2x * Gx[p] + b2y * Gy[p];
		}
	}

	if (deforming) {
		maxGx = maxabs( Gx, GridSize_l );
		maxGy = maxabs( Gy, GridSize_l );
		maxG  = maxGx > maxGy ? maxGx : maxGy;

		dt_scale = fabs(maxG) > EPS ? dt/maxG : 0;
	}
	else {
		maxG = fabs(trans[0]) > fabs(trans[1]) ? fabs(trans[0]) : fabs(trans[1]);
		dt_scale = fabs(maxG) > EPS ? 0.2/maxG : 0;
	}

	//		printf("updating backward warps; dt_scale=%f\n", dt_scale);

	for (int p=0; p<GridSize_l; p++) {
		bwarp_x[p]+=dt_scale*dbwx[p];
		bwarp_y[p]+=dt_scale*dbwy[p];
	}

	return;
}


void Motion_Sobolev::updateForwardMap( int l, double dt_scale )
{
	int GridSize_l = GridSize[l];
	int XSize_l = XSize[l];
	int YSize_l = YSize[l];
	double fw_x, fw_y, delta_fwx, delta_fwy;
	int fw_x_i, fw_y_i;
	int pfw;
	double Ginterpx, Ginterpy;
	ARRAYTYPE *fwarp_x = forward_warp[l][0];
	ARRAYTYPE *fwarp_y = forward_warp[l][1];
	ARRAYTYPE *Gx = SobolevGrad;
	ARRAYTYPE *Gy = SobolevGrad + GridSize_l;

	//	printf("forward warps\n");

	// compute forward warp differential
#ifdef PARALLELIZE_CODE
		#pragma omp parallel for private(fw_x, fw_y, fw_x_i, fw_y_i, delta_fwx, delta_fwy, pfw, Ginterpx, Ginterpy)
#endif

	for (BLOCKS(block)) {

		BLOCK_ENDS( GridSize_l );

		for (PIXELS(p)) {
			fw_x = fwarp_x[p];
			fw_y = fwarp_y[p];

			fw_x = fw_x < 0         ?         0 : fw_x;
			fw_x = fw_x > XSize_l-1 ? XSize_l-1 : fw_x;

			fw_y = fw_y < 0         ?         0 : fw_y;
			fw_y = fw_y > YSize_l-1 ? YSize_l-1 : fw_y;

			fw_x_i = int(fw_x);
			fw_y_i = int(fw_y);

			delta_fwx = fw_x - fw_x_i;
			delta_fwy = fw_y - fw_y_i;

			pfw = fw_x_i + XSize_l * fw_y_i;

			Ginterpx = (1-delta_fwx) * (1-delta_fwy) * Gx[pfw]                                  +
               			   delta_fwx * (1-delta_fwy) * ( edgebits[pfw]&XPOS ? Gx[pfw] : Gx[pfw+1] ) +
			           (1-delta_fwx) *    delta_fwy  * ( edgebits[pfw]&YPOS ? Gx[pfw] : Gx[pfw+XSize_l]) +
			               delta_fwx *    delta_fwy  * ( edgebits[pfw]&(XPOS|YPOS) ? Gx[pfw] : Gx[pfw+1+XSize_l]);

			Ginterpy = (1-delta_fwx) * (1-delta_fwy) * Gy[pfw]                                  +
                		   delta_fwx * (1-delta_fwy) * ( edgebits[pfw]&XPOS ? Gy[pfw] : Gy[pfw+1] ) +
			           (1-delta_fwx) *    delta_fwy  * ( edgebits[pfw]&YPOS ? Gy[pfw] : Gy[pfw+XSize_l]) +
			               delta_fwx *    delta_fwy  * ( edgebits[pfw]&(XPOS|YPOS) ? Gy[pfw] : Gy[pfw+1+XSize_l] );

			dfwx[p]=Ginterpx;
			dfwy[p]=Ginterpy;
		}
	}

	//	printf("updating forward warps\n");

	// update forward warp
	for (int p=0; p<GridSize_l; p++) {
		fwarp_x[p] -= dt_scale * dfwx[p];
		fwarp_y[p] -= dt_scale * dfwy[p];
	}

	return;
}


void Motion_Sobolev::evolveForwardBackwardWarps(int l, int iters, double dt)
{
	int GridSize_l = GridSize[l];
	int XSize_l = XSize[l];
	int YSize_l = YSize[l];
	ARRAYTYPE trans[2];

	ARRAYTYPE *bwarp_x = backward_warp[l][0];
	ARRAYTYPE *bwarp_y = backward_warp[l][1];
	ARRAYTYPE *fwarp_x = forward_warp[l][0];
	ARRAYTYPE *fwarp_y = forward_warp[l][1];
	ARRAYTYPE *Gx = SobolevGrad;
	ARRAYTYPE *Gy = SobolevGrad + GridSize_l;


	double dt_scale, energy, energy_prev, diff_energy;

	int deform_init_valid = 0;
	int count_iters       = 0;
	int stop_translating  = 0;
	int deforming         = 0;

	double cg_time=0, cg_time_2 = 0;

	do {

	//	printf("computing gradients\n");
		computeL2Gradient( l );
		computeGradTranslation ( l, trans );
		// && (fabs(trans[0]) > TRANSLATION_EPS || fabs(trans[1]) > TRANSLATION_EPS)
       // stop_translating=1;

		if ( !stop_translating ) {
#ifdef DEBUG
			printf("translating; t = (%f, %f) ...  \n", trans[0], trans[1]);
#endif
			for (int p=0; p<GridSize_l; p++) {
				Gx[p]=trans[0];
				Gy[p]=trans[1];
			}
			deforming = 0;
		}
		else {
#ifdef DEBUG
			printf("deforming ... \n");
#endif
			cg_time += computeSobolevDeformation( l, trans, deform_init_valid );
			deform_init_valid = 1;
			deforming         = 1;
			stop_translating  = 0;
		}

		// update warps in Sobolev Gradient direction

	//		printf("backward warp computation\n");

		updateBackwardMap( l, deforming, dt, dt_scale, trans );

		updateForwardMap( l, dt_scale );

		// compute quantitities related to convergence

		if ( count_iters == 0 ) {
			energy_prev = DBL_MAX;
			diff_energy = DBL_MAX;
		}
		else {
			energy_prev = energy;
		}

		energy = computeEnergy( l );

		if (count_iters !=0 ) diff_energy = (energy_prev - energy) / energy_prev;

#ifdef WRITE_PRINT_STATEMENTS

		if ( count_iters ==0 )
			printf("iters = %d,\tenergy = %f,\tenergy_prev = inf,\tdiff_energy = inf\n", count_iters, energy );
		else 
			printf("iters = %d,\tenergy = %f,\tenergy_prev = %f,\tdiff_energy = %f\n", count_iters, energy, energy_prev, diff_energy );

		if (deforming) printf("Deforming\n"); else printf("Translating\n");
#endif

		count_iters++;

		if ( (deforming && diff_energy < MAX_ENERGY_DIFF_MOTION) || count_iters == iters ) {
			if ( diff_energy < 0 ) {
				// backtrack to previous state since we went to higher energy state

				for (int p=0; p<GridSize_l; p++) {
					fwarp_x[p] += dt_scale * dfwx[p];
					fwarp_y[p] += dt_scale * dfwy[p];
					bwarp_x[p] -= dt_scale*dbwx[p];
					bwarp_y[p] -= dt_scale*dbwy[p];
				}
				
				energy = energy_prev;
			}
			
			break;
		}

		if ( !deforming && diff_energy < TRANSLATION_EPS ) {
			if ( diff_energy < 0 ) {
				// backtrack to previous state since we went to higher energy state

				for (int p=0; p<GridSize_l; p++) {
					fwarp_x[p] += dt_scale * dfwx[p];
					fwarp_y[p] += dt_scale * dfwy[p];
					bwarp_x[p] -= dt_scale*dbwx[p];
					bwarp_y[p] -= dt_scale*dbwy[p];
				}
				
				energy = energy_prev;
			}
			stop_translating = 1;
		}

	} while ( 1 );
#ifdef WRITE_PRINT_STATEMENTS
	printf( "total cg_time = %f + %f = %f\n", cg_time, cg_time_2, cg_time + cg_time_2 );
#endif
	return;
}


void Motion_Sobolev::initializeWarps(int l, ARRAYTYPE *fwarpi, ARRAYTYPE *bwarpi)
{
	int XSize_l = XSize[l];
	int YSize_l = YSize[l];
	int GridSize_l = GridSize[l];
	ARRAYTYPE **fwarp=forward_warp[l];
	ARRAYTYPE **bwarp=backward_warp[l];

	for (int y=0, p=0; y<YSize_l; y++) {
		for (int x=0; x<XSize_l; x++, p++) {
			fwarp[0][p]=fwarpi[p] + x;
			fwarp[1][p]=fwarpi[p+GridSize_l] + y;
			bwarp[0][p]=bwarpi[p] + x;
			bwarp[1][p]=bwarpi[p+GridSize_l] + y;
		}
	}

	return;
}


void Motion_Sobolev::computeResidual(int l)
{
	int XSize_l    = XSize[l];
	int YSize_l    = YSize[l];
	int GridSize_l = GridSize[l];
	char *region_l = region[0];
	ARRAYTYPE diff[FSize];
	double diff_norm_sq;


	for (int p=0; p<GridSize_l; p++) {
		if (region_l[p]) {
			double xfw = forward_warp[l][0][p];
			double yfw = forward_warp[l][1][p];

			if ( xfw > XSize_l - 1 || yfw > YSize_l - 1 || xfw < 0 || yfw < 0 ) {
				residual[p] = occlusion_threshold;
			}
			else {
				diff_I1warp_I0( l, p, diff );
				diff_norm_sq = vnormsq( diff );
				residual[p] = diff_norm_sq; // > occlusion_threshold ? occlusion_threshold : diff_norm_sq;
			}
		}
		else {
			residual[p]=0;
		}
	}

	for (int p=GridSize_l; p<GridSize[0]; p++) {
		residual[p]=0;
	}

	return;
}


/*
 * computes ( I1(w(x)) - I0(x) ) * 1_{ R \ occlusion } (x)
 */
void Motion_Sobolev::computeResidual2(int l)
{
	int XSize_l    = XSize[l];
	int YSize_l    = YSize[l];
	int GridSize_l = GridSize[l];
	char *region_l = region[0];
	ARRAYTYPE diff[FSize];
	double diff_norm_sq;


	for (int p=0; p<GridSize_l; p++) {
		if (region_l[p]) {
			double xfw = forward_warp[l][0][p];
			double yfw = forward_warp[l][1][p];

			if ( xfw > XSize_l - 1 || yfw > YSize_l - 1 || xfw < 0 || yfw < 0 ) {
				for (int f=0; f<FSize; f++)
					residual_vec[p+f*GridSize[0]] = 100; //occlusion_threshold; //0;
			}
			else {
				diff_I1warp_I0( l, p, diff );
				diff_norm_sq = vnormsq( diff );
				for (int f=0; f<FSize; f++)
					residual_vec[p+f*GridSize[0]] = diff[f]; //diff_norm_sq >= occlusion_threshold ? 0 : diff[f];
			}
		}
		else {
			for (int f=0; f<FSize; f++)
				residual_vec[p+f*GridSize[0]] = 0;
		}
	}

	for (int p=GridSize_l; p<GridSize[0]; p++) {
		for (int f=0; f<FSize; f++)
			residual_vec[p+f*GridSize[0]] = 0;
	}

	return;
}


double Motion_Sobolev::computeEnergy(int l)
{
	int XSize_l    = XSize[l];
	int YSize_l    = YSize[l];
	int GridSize_l = GridSize[l];
	char *region_l = region[0];
	ARRAYTYPE diff[FSize];
	double diff_norm_sq;
	double energy;
	double sum[Nblocks];

	energy=0;

	for (BLOCKS(block)) sum[block]=0;

#ifdef PARALLELIZE_CODE
    #pragma omp parallel for private (diff, diff_norm_sq)
#endif

	for (BLOCKS(block)) {
		BLOCK_ENDS( GridSize_l );

		double sum_tmp=0;
		double xfw, yfw;

		for (PIXELS(p)) {
			if (region_l[p]) {
				xfw = forward_warp[l][0][p];
				yfw = forward_warp[l][1][p];

				if ( xfw > XSize_l - 1 || yfw > YSize_l - 1 || xfw < 0 || yfw < 0 ) {
					sum_tmp += occlusion_threshold;
				}
				else {
					diff_I1warp_I0( l, p, diff );
					diff_norm_sq = vnormsq( diff );
					if (diff_norm_sq < occlusion_threshold)
						sum_tmp += diff_norm_sq;
					else 
						sum_tmp +=occlusion_threshold;
				}

			}
		}

		sum[block]=sum_tmp;
	}

	for (BLOCKS(block)) energy+=sum[block];

	return energy;
}