function [fwarp,bwarp,region_warp,residual] = pyramidMotion( J0, J1, region, settings, plotlevs )
% J0, J1 - rgb normalized between 0 & 1
% region - indicator function of region
% 
	region = double(region);
	[m,n,o] = size(J0);
settings.levels = max( floor(log2(min(m,n)))-3, 1 );
	for l = settings.levels : -1 : 1,
		scale = 2^(-l+1);
		size_l = floor([m,n]*scale);
		J0_l = imresize( J0, size_l ); 
		J0_l( J0_l < 0 ) = 0;  J0_l( J0_l > 1 ) = 1;
		J1_l = imresize( J1, size_l );
		J1_l( J1_l < 0 ) = 0;  J1_l( J1_l > 1 ) = 1;

		region_l = int8(imresize( region, size_l ) >=0.5);

		levels_mg = settings.levels - l + 1;
  		beta = settings.beta * 4^(-l+1);
%         beta = settings.beta;
		if l == settings.levels,
			[fwarp,bwarp,region_warp,residual] = mexMotionEstimation( J0_l, J1_l, region_l, settings.dt, settings.iters, levels_mg, settings.steps_v, settings.steps_rel_u, settings.steps_rel_d, settings.occlusion_threshold, beta );
		else
			fwarp_up =zeros([size_l 2]);
			bwarp_up =zeros([size_l 2]);
			fwarp_up(:,:,1) = 2*imresize( fwarp(:,:,1), size_l );
			fwarp_up(:,:,2) = 2*imresize( fwarp(:,:,2), size_l );
			bwarp_up(:,:,1) = 2*imresize( bwarp(:,:,1), size_l );
			bwarp_up(:,:,2) = 2*imresize( bwarp(:,:,2), size_l );

			[fwarp,bwarp,region_warp,residual] = mexMotionEstimation( J0_l, J1_l, region_l, settings.dt, settings.iters, levels_mg, settings.steps_v, settings.steps_rel_u, settings.steps_rel_d, settings.occlusion_threshold, beta, fwarp_up, bwarp_up );
		end

% 		if plotlevs,
% 			plotWarps( J0_l, J1_l, fwarp, bwarp, region_l, region_warp );
% 		end
		
	end