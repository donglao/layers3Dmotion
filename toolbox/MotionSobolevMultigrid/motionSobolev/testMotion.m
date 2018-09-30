% dData = '~/work/Datasets_Vision/moseg_dataset/marple7/';

addpath('./flow-code-matlab/');
addpath('./multigrid2/');

levels=6;
occlusion_threshold = 0.1;
beta=1e-7;
dt=0.48;
iters = 1000;

% I0=imread( strcat( dData, 'marple7_001.jpg') );
% I1=imread( strcat( dData, 'marple7_003.jpg') );

% region = imread( strcat( dData,'marple7_001.pgm' ) );

I0 = org{10};
I1 = org{11};
region = labels{10};


%I0 = imresize( I0, 0.5 );
%I1 = imresize( I1, 0.5 );
%region = imresize( region, 0.5 );

%region = 255 - region;
%region = 255*ones( size(region));
region = int8(region >0);
%region = int8( ones( size(region) ) );

label = double(region == 1);

se =  strel( 'square', 15*2 );
%region = int8( imdilate( label, se ) );


J0 = im2double( I0 );
J1 = im2double( I1 );

[fwarp,bwarp,region_warp,residual] = mexMotionEstimation( J0, J1, region, dt, iters, levels, occlusion_threshold, beta );

%Iwarp = imageWarp( J1, fwarp );
%fwarp_img=flowToColor( fwarp );
%figure; imagesc( fwarp_img ); axis image;

%keyboard;

%[fwarp,bwarp,region_warp,residual] = mexMotionEstimation( J0, J1, region, dt, iters, levels, occlusion_threshold, fwarp, bwarp );
fwarp_masked = fwarp;
fwarp_masked(:,:,1) = fwarp_masked(:,:,1) .* double( region );
fwarp_masked(:,:,2) = fwarp_masked(:,:,2) .* double( region );

bwarp_masked = bwarp;
bwarp_masked(:,:,1) = bwarp_masked(:,:,1) .* double( region_warp );
bwarp_masked(:,:,2) = bwarp_masked(:,:,2) .* double( region_warp );

fwarp_img=flowToColor( fwarp_masked );
bwarp_img=flowToColor( bwarp_masked );

figure;

subplot(221); imagesc(I0); axis image; title('I0');
subplot(222); imagesc(I1); axis image; title('I1');
subplot(223); imagesc(fwarp_img); axis image; title('Forward warp displacement');
subplot(224); imagesc(bwarp_img); axis image; title('Backward warp displacement');

figure;

subplot(121); show_mask_on_image2( I0, double(region), 0.2, 'red' ); axis image;
title('Image 1; initial region');
subplot(122); show_mask_on_image2( I1, double(region_warp), 0.2, 'red' ); axis image;
title('Image 2; Sobolev warped region');

figure;

subplot(121); imagesc(residual); axis image; colormap gray; colorbar;
title('Residual defined on region');
subplot(122); imagesc(double(residual>occlusion_threshold)); axis image; colormap gray; colorbar;
title('Occlusion map');