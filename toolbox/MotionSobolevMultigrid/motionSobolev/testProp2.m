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

fGT1    = strcat( dData, 'marple7_001.pgm' );
fIM1    = strcat( dData, 'marple7_001.jpg' );
fIM2    = strcat( dData, 'marple7_003.jpg' );

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

Lnew = propagateSegmentation(I1_LAB, I2_LAB, I1_RGB, I2_RGB, L0, option);

figure;
subplot(121); show_mask_on_image2( uint8(I2), double(Lnew==0), 0.3, 'red' ); axis image; title('propagated label 0 in I2');
subplot(122); show_mask_on_image2( uint8(I2), double(Lnew==1), 0.3, 'red' ); axis image; title('propagated label 1 in I2');
