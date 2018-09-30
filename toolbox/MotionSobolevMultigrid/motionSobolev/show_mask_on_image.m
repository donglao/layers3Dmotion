function show_mask_on_image(I, mask, alpha, color, mask2, alpha2, color2)
% 	SHOW_MASK_ON_IMAGE   Short description
% 		SHOW_MASK_ON_IMAGE(I, MASK, ALPHA, COLOR)
% 

  if nargin < 3; alpha = 1; end
  if nargin < 4; color = 'red'; end
		
  [M N D] = size(I);
	
  if D > 1
    imshow(rgb2gray(I));
  else
    imagesc(I); colormap gray;
  end

  hold on;
  h=imagesc(solid_color_image([M N], color));
  set(h, 'AlphaData', alpha*mask);


  if nargin > 4
    h = imagesc(solid_color_image([M N], color2));
    set(h, 'AlphaData', alpha2*mask2);
  end

end %  function	 ctionshowmask()