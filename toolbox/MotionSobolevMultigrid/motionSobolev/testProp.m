dData = '~/work/Datasets_Vision/moseg_dataset/marple7/';

addpath('./flow-code-matlab/');
addpath('../propagation/');

levels=2;
occlusion_threshold = 0.08;
dt=0.2;
iters = 250;
sigma_phi = 0.85;

if 0,
I0=imread( strcat( dData, 'marple7_001.jpg') );
I1=imread( strcat( dData, 'marple7_003.jpg') );

region = imread( strcat( dData,'marple7_001.pgm' ) );

I0 = imresize( I0, 0.5 );
I1 = imresize( I1, 0.5 );
region = imresize( region, 0.5 );
%region = 255 - region;
region0 = uint8(region == 255);
region1 = uint8(region ~= 255);

J0 = double( I0 )/255;
J1 = double( I1 )/255;
end

dWork = '/Users/sundarg/work/motion_hong/';
dCode = [dWork 'FastLabelScaleSpace_h/'];
dMex  = [dWork 'FastLabelScaleSpace_h/mex/'];
dMotion = [dWork 'FastLabelScaleSpace_h/motionSobolev/'];
dData = ['~/work/Datasets_Vision/moseg_dataset/marple7/'];
dUtil  = [dWork 'FastLabelScaleSpace_h/util/'];

addpath(dCode);
addpath(dMex);
addpath(dMotion);
addpath(dUtil);

%run([dMex 'mexCompile.m'])

fGT1    = '~/work/Datasets_Vision/moseg_dataset/marple7/marple7_001.pgm';
fIM1    = '~/work/Datasets_Vision/moseg_dataset/marple7/marple7_001.jpg';
fIM2    = '~/work/Datasets_Vision/moseg_dataset/marple7/marple7_002.jpg';

G1 = imread([fGT1]);
I1 = imread([fIM1]);
I2 = imread([fIM2]);

G1 = imresize(G1, 0.5, 'nearest');
I1 = imresize(I1, 0.5, 'bicubic');
I2 = imresize(I2, 0.5, 'bicubic');

colorTransform = makecform('srgb2lab');
I1c = lab2double(applycform(I1, colorTransform));
I2c = lab2double(applycform(I2, colorTransform));
I1c = double(I1c);
I2c = double(I2c);

[nRow, nCol, nVec] = size(I1c);

I1_LAB = zeros(nRow, nCol, nVec);
I2_LAB = zeros(nRow, nCol, nVec);

for v = 1 : nVec
    minVal  = min( min( [I1c(:,:,v), I2c(:,:,v)] ) );
    maxVal  = max( max( [I1c(:,:,v), I2c(:,:,v)] ) );

    I1_LAB(:,:,v)  = (I1c(:,:,v) - minVal) / (maxVal - minVal);
    I2_LAB(:,:,v)  = (I2c(:,:,v) - minVal) / (maxVal - minVal);
end

I1 = double(I1);
I2 = double(I2);

minVal  = min( [I1(:) ; I2(:)] );
maxVal  = max( [I1(:) ; I2(:)] );
I1_RGB = (I1 - minVal ) / (maxVal - minVal);
I2_RGB = (I2 - minVal ) / (maxVal - minVal);

l       = unique(G1);
L0      = zeros(nRow, nCol);
nLabel  = length(l);

for i = 1 : nLabel
    idx = find(G1 == l(i));
    L0(idx) = i;
end

region0 = int8(L0==1);
region1 = int8(L0==2);

[fwarp0,bwarp0,region_warp0,residual0] = mexMotionEstimation( I1_RGB, I2_RGB, region0, dt, iters, levels, occlusion_threshold );
[fwarp1,bwarp1,region_warp1,residual1] = mexMotionEstimation( I1_RGB, I2_RGB, region1, dt, iters, levels, occlusion_threshold );

[m,n,o]=size(I1);

phi = zeros( m, n, 2 );
phi(:,:,1) = double( region0 );
phi(:,:,2) = double( region1 );

bwarp_full = zeros( m, n, 2, 2 );
bwarp_full(:,:,:,1) = bwarp0;
bwarp_full(:,:,:,2) = bwarp1;

%h = fspecial('gaussian', 2*sigma+1, sigma);

occlusion_map = zeros( m, n, 2 );
%occlusion_map(:,:,1) = imfilter( residual0, h ) > occlusion_threshold;
%occlusion_map(:,:,2) = imfilter( residual1, h ) > occlusion_threshold;
occlusion_map(:,:,1) = residual0 >= occlusion_threshold;
occlusion_map(:,:,2) = residual1 >= occlusion_threshold;

%keyboard;

[label,phi_warped] = mexComputePropagatedSegmentation( phi, bwarp_full, occlusion_map, sigma_phi );

figure;

subplot(131); show_mask_on_image2( uint8(I1), double(region0), 0.3, 'red' ); axis image; title('I0 with initial region');
subplot(132); imagesc( uint8(I2) ); axis image; title('I1');
subplot(133); imagesc( label ); axis image; title('label; propagated segmentation');

figure;

subplot(131); show_mask_on_image2( uint8(I2), double(label==0), 0.3, 'red' ); axis image;
title('label 0 propagated');
subplot(132); show_mask_on_image2( uint8(I2), double(label==1), 0.3, 'red' ); axis image;
title('label 1 propagated');
subplot(133); show_mask_on_image2( uint8(I2), double(label==-1), 0.3, 'red' ); axis image;
title('disocclusion');

% fill disocclusion

option.plot                 = 1;    % flag for plotting intermediate results
option.plotStep             = 3;    % flag for plotting intermediate results
option.iter                 = 2000;  % number of iterations at optimisation
option.Histogram.bin        = 10;   % number of bins for histogram
option.ScaleSpace.scale     = 3;	% number of scales at scale space
option.Subpixel.tau         = 0.5;	% step size of the subpixel on indicator
option.Regularity.lambda    = 0.1;	% weight on the laplacian of indicator
option.Energy.weight        = 1;	% weight on the scale at energy computation
option.Histogram.local      = 20;   % size of region for local histogram

option.Motion.scale               = 1;
option.Motion.level               = 1;
option.Motion.Occlusion.threshold = 0.08;
option.Motion.dt 				  = 0.49;
option.Motion.iters 			  = 10000;

option.Motion.Ambiguity.k = 0.1;  % threshold for the standard deviation
option.Motion.Ambiguity.r = 10;     % window size for computing s.t.d.
option.Motion.Ambiguity.b = 0.05;  % threshold for the minimum residual

option.Motion.Ambiguity.k = 0.1;  % threshold for the standard deviation
option.Motion.Ambiguity.r = 3;     % window size for computing s.t.d.
option.Motion.Ambiguity.b = 0.08;  % threshold for the minimum residual

Lnew = getRegionForce(I1_LAB, I2_LAB, I1_RGB, I2_RGB, label, option);

figure;
subplot(121); show_mask_on_image2( uint8(I2), double(Lnew==0), 0.3, 'red' ); axis image;
subplot(122); show_mask_on_image2( uint8(I2), double(Lnew==1), 0.3, 'red' ); axis image;
%imagesc(Lnew); axis image; colormap gray;