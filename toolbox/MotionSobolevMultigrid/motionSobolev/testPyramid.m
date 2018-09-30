dData = '~/work/Datasets_Vision/moseg_dataset/people2/';
dData = '~/work/Datasets_Vision/davis/JPEGimages/1080p/cows/';
dDataGT = '~/work/Datasets_Vision/davis/Annotations/1080p/cows/';

addpath('./flow-code-matlab/');
addpath('./multigrid2/');

settings.steps_v     =1;    % multigrid, # vcycles
settings.steps_rel_u =1;    % multigrid, # upward relaxation steps
settings.steps_rel_d =1;    % multigrid, # downward relaxation steps
settings.occlusion_threshold = 0.1;
settings.beta        =1e-5;   % smoothness in flow; lower is smoother
settings.dt          =0.48;   % step for warp updates
settings.iters       =1000;   % max iterations
plotlevs             =1;      % whether you want to plot result

I0=imread( strcat( dData, '00000.jpg') );
I1=imread( strcat( dData, '00001.jpg') );
region = imread( strcat( dDataGT,'00000.png' ) );

[nRows,nCols,nFeatures] = size(I0);

settings.levels = max( floor(log2(min(nRows,nCols)))-3, 1 );

J0 = double( I0 )/255;
J1 = double( I1 )/255;

%region = 255 - region;
%region = 255*ones( size(region));
region = int8(region == 255);
label  = double(region == 1);

se =  strel( 'square', 15*2 );
%region = int8( imdilate( label, se ) );

% initialization for flow
fwarp_i = zeros( [nRows, nCols 2] );
bwarp_i = zeros( [nRows, nCols 2] );

% call this one if you have a good initialization for the flow
[fwarp,bwarp,region_warp,residual] = multigridMotion( J0, J1, region, settings, plotlevs, fwarp_i, bwarp_i );

% call this one if you don't have an initialization
[fwarp,bwarp,region_warp,residual] = pyramidMotion( J0, J1, region, settings, plotlevs );