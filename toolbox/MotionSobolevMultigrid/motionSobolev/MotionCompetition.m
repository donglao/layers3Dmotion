function MotionCompetition( I0, I1, label, nIters )
% label - 0 to nLabels - 1

	addpath('./flow-code-matlab/');

	band_length = 1;
	iters_flow  = 10000;
	dt_motion = 0.49;
	occlusion_threshold = 0.1;
	levels = 1;

	reg_phi = 0.01;
	dt_phi = 0.4;

	nLabels = double( max( label(:) ) + 1 );

	[m,n,o] = size(I0);

	phi            = zeros( m, n, nLabels );
	bwarp_all      = zeros( m, n, 2, nLabels );
	fwarp_all      = zeros( m, n, 2, nLabels );
	residual_all   = zeros( m, n, nLabels );
	region_ext_all = zeros( m, n, nLabels );
	region_warp_all = zeros( m, n, nLabels );

	for l = 1 : nLabels,
		phi(:,:,l) = double( label==l-1 );
	end

	for i = 1 : nIters,

		bPlot = (mod(i,10) == 0);

		[val, L] = max( phi, [], 3);

		for l = 1 : nLabels,
			region = double(L==l);
			[region_ext, band, ind]=extendRegion( region, 1, band_length );

			region_ext_all(:,:,l) = region_ext;

			[fwarp, bwarp, region_warp, residual] = mexMotionEstimation( I0, I1, int8(region_ext), dt_motion, iters_flow, levels, occlusion_threshold, fwarp_all(:,:,:,l), bwarp_all(:,:,:,l) );

			region_warp_all(:,:,l) = region_warp;

%	        [fwarp, bwarp, region_warp, residual] = mexMotionEstimation( I0, I1, int8(region_ext), dt_motion, iters_flow, levels, occlusion_threshold );

			fwarp_all(:,:,:,l)  = fwarp;
			bwarp_all(:,:,:,l)  = bwarp;
			residual_all(:,:,l) = residual; %min( residual, occlusion_threshold );
		end

		if bPlot
		subplot(321); show_mask_on_image2( I0, double(label==1), 0.3, 'red' ); axis image; title('initial');
		subplot(322); show_mask_on_image2( I0, double(L==2), 0.3, 'red' ); axis image;
		maxflow_0 = max( max( max( abs(fwarp_all(:,:,:,1)) ) ) );
		maxflow_1 = max( max( max( abs(fwarp_all(:,:,:,2)) ) ) );
		minflow_0 = min( min( min( abs(fwarp_all(:,:,:,1)) ) ) );
		minflow_1 = min( min( min( abs(fwarp_all(:,:,:,2)) ) ) );
		str0 = sprintf('maxflow=%f, minflow=%f', maxflow_0, minflow_0);
		str1 = sprintf('maxflow=%f, minflow=%f', maxflow_1, minflow_1);
		subplot(323); imagesc(residual_all(:,:,1)); colorbar; axis image;
		subplot(324); imagesc(residual_all(:,:,2)); colorbar; axis image;

		fwarp_flow = fwarp_all(:,:,:,1).*repmat( phi(:,:,1), [1, 1, 2] ) + fwarp_all(:,:,:,2).*repmat( phi(:,:,2), [1, 1, 2] );
		bwarp_flow = bwarp_all(:,:,:,1).*repmat( region_warp_all(:,:,1), [1, 1, 2] ) + bwarp_all(:,:,:,2).*repmat( region_warp_all(:,:,2), [1, 1, 2] );

		subplot(325); imagesc( flowToColor( fwarp_flow ) ); axis image; colorbar;
		title(str0);
		subplot(326); imagesc(flowToColor( bwarp_flow ) ); axis image; colorbar;
		title(str1);

		flow1 = sqrt( fwarp_all(:,:,1,1).^2 + fwarp_all(:,:,2,1).^2 ); flow1 = flow1(:)';
		flow2 = sqrt( fwarp_all(:,:,1,2).^2 + fwarp_all(:,:,2,2).^2 ); flow2 = flow2(:)';

%		subplot(427); imagesc(flowToColor(  ) ); axis image; colorbar;
%		title(str0);
%    	subplot(428); imagesc(flowToColor(  ); axis image; colorbar;
%		title(str1);

%		subplot(427); hist( flow1, 50 );
%		subplot(428); hist( flow2, 50 );

%		figure(2000);
%		subplot(121); show_mask_on_image2( I1, double( region_warp_all(:,:,1) ), 0.3, 'red' ); axis image; title('warped region');
%		subplot(122); show_mask_on_image2( I1, double( region_warp_all(:,:,2) ), 0.3, 'red' ); axis image; title('warped region');

		%subplot(427); imagesc(phi(:,:,1)); axis image; colorbar;
		%subplot(428); imagesc(phi(:,:,2)); axis image; colorbar;
		drawnow;
		end


		energy =0;
		for l = 1 : nLabels,
			rs = residual_all(:,:,l).*phi(:,:,l) ;
			energy = energy + sum(rs(:));
		end

		fprintf( 1, 'iter = %d,\tenergy = %f\n', i, energy );

		% compute max force

		maxF = 0;

		for l1 = 1 : nLabels,
			for l2 = 1 : nLabels,
				if l1 ~= l2,
					band = region_ext_all(:,:,l1) & region_ext_all(:,:,l2);
					ind = find(band==1);

					%keyboard;

					F = residual_all( ind + (l1-1)*m*n ) - residual_all( ind + (l2-1)*m*n );
					maxF = max( max(abs(F)), maxF );
				end
			end
		end

		dt_phi_scale = min( dt_phi / maxF, 0.25/reg_phi );

		for l1 = 1 : nLabels,
			d_phi = zeros( m, n );
			for l2 = 1 : nLabels,
				if l1 ~= l2,
					band = region_ext_all(:,:,l1) & region_ext_all(:,:,l2);
					ind = find(band==1);

					d_phi(ind) = residual_all( ind + (l1-1)*m*n ) - residual_all( ind + (l2-1)*m*n );
				end
			end
			phi(:,:,l1) = phi(:,:,l1) + dt_phi_scale * ( -d_phi + reg_phi * del2( phi(:,:,l1) ) );
			phi(:,:,l1) = max( phi(:,:,l1), 0 );
			phi(:,:,l1) = min( phi(:,:,l1), 1 );
		end

	end
