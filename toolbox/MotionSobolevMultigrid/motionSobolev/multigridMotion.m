function [fwarp,bwarp,region_warp,residual] = multigridMotion( J0, J1, region, settings, plotwarps, fwarp_i, bwarp_i )
% J0, J1 - rgb normalized between 0 & 1
% region - indicator function of region
% if last two arguments are omitted, the those are initilalized with zero

	if nargin < 7,
		[nRows, nCols, nDep] = size(J0);
		fwarp_i = zeros( [nRows, nCols 2] );
		bwarp_i = zeros( [nRows, nCols 2] );
    end
    
%     [fwarp,bwarp,region_warp,residual] = mexMotionEstimation( J0, J1, region, settings.dt, settings.iters, settings.levels, settings.steps_v, settings.steps_rel_u, settings.steps_rel_d, settings.occlusion_threshold, settings.beta, fwarp_i, bwarp_i );
[fwarp,bwarp,region_warp,residual] = mexMotionEstimation( J0, J1, region, settings.dt, settings.iters, settings.levels, settings.steps_v, settings.steps_rel_u, settings.steps_rel_d, settings.occlusion_threshold, settings.beta, fwarp_i, bwarp_i );

% 	if plotwarps,
% 		plotWarps( J0, J1, fwarp, bwarp, region, region_warp );
% 	end
end