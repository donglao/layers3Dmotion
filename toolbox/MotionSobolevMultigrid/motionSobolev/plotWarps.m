function plotWarps( I0, I1, fwarp, bwarp, region, region_warp )

	figure;

	fwarp_img = convertFlowToColorMasked( fwarp, region );
	bwarp_img = convertFlowToColorMasked( bwarp, region_warp );

	subplot(321); imagesc(uint8(I0*255)); axis image; title('I0');
	subplot(322); imagesc(uint8(I1*255)); axis image; title('I1');
	subplot(323); imagesc(fwarp_img); axis image; title('Forward warp displacement');
	subplot(324); imagesc(bwarp_img); axis image; title('Backward warp displacement');

	subplot(325); show_mask_on_image2( I0, double(region), 0.2, 'red' ); axis image;
	title('Image 1; initial region');
	subplot(326); show_mask_on_image2( I1, double(region_warp), 0.2, 'red' ); axis image;
	title('Image 2; warped region');