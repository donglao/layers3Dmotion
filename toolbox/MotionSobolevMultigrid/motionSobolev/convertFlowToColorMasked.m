function warp_img = convertFlowToColorMasked( warp, region )

	warp_masked = warp;
	warp_masked(:,:,1) = warp(:,:,1) .* double( region );
	warp_masked(:,:,2) = warp(:,:,2) .* double( region );
	warp_img=flowToColor( warp_masked );