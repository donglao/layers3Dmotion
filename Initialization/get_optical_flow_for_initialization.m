function [f_flow, b_flow] = get_optical_flow_for_initialization( org )
addpath( genpath( './Initialization/3rdparty/flow_code_v2' ));
parfor i = 1:length(org)-1
    f_flow{i} = estimate_flow_interface(org{i}, org{i + 1}, 'classic+nl-fast');
    b_flow{i} = estimate_flow_interface(org{i + 1}, org{i}, 'classic+nl-fast');
end
end