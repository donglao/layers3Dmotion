function [flow_f, flow_b, region_warp, residual] = warp_matlab(J0, J1, region, flow_f, flow_b)
settings.steps_v     =1;    % multigrid, # vcycles
settings.steps_rel_u =1;    % multigrid, # upward relaxation steps
settings.steps_rel_d =1;    % multigrid, # downward relaxation steps
settings.occlusion_threshold = 0.1;
settings.beta        =1e-10;   % smoothness in flow; lower is smoother
settings.dt          =0.48;   % step for warp updates
settings.iters       =1000;   % max iterations
plotlevs             =1;      % whether you want to plot result

if nargin > 3
    flow_f(:,:,[1,2])=flow_f(:,:,[2,1]);
    flow_b(:,:,[1,2])=flow_b(:,:,[2,1]);
    flow_f = flow_repaint(flow_f);
    flow_b = flow_repaint(flow_b);
end
region = int8(region);
J0 = im2double(J0);
J1 = im2double(J1);
[h, w] = size (J0(:,:,1));
settings.levels = max( floor(log2(min(h,w)))-3, 1 );
    if nargin > 3
        [flow_f, flow_b, region_warp,residual] = multigridMotion( J0, J1, region, settings, plotlevs, flow_f, flow_b );
    else
        [flow_f, flow_b, region_warp,residual] = pyramidMotion( J0, J1, region, settings, plotlevs );
    end
flow_f(:,:,[1,2])=flow_f(:,:,[2,1]);
flow_b(:,:,[1,2])=flow_b(:,:,[2,1]);
end

function flow = flow_repaint(flow, region)
if nargin == 2
    region = double(region);
    region(region == 0) = nan;
    for i = 1:2
        flow(:,:,i) = inpaint_nans(flow(:,:,i).*region, 5);
    end
else
    for i = 1:2
        flow(:,:,i) = inpaint_nans(flow(:,:,i), 5);
    end
end
end