function show_mask_on_image2(I, mask, alpha, color, mask2, alpha2, color2)
% 	SHOW_MASK_ON_IMAGE   Short description
% 		SHOW_MASK_ON_IMAGE(I, MASK, ALPHA, COLOR)
% 

  if nargin < 3; alpha = 1; end
  if nargin < 4; color = 'red'; end
		
  [M N D] = size(I);
	
%  if D > 1
%    imshow(rgb2gray(I));
%  else
%    imagesc(I); colormap gray;
%  end


  h1 = imagesc(solid_color_image([M N], color));
  % set(h1, 'AlphaData', 1-alpha*mask);
  hold on;
  
  if nargin > 4
     h2 = imagesc(solid_color_image([M N], color2));
     set(h2, 'AlphaData', 0.5*(mask2&mask)+(1-mask));
  end
  
  
  % if D > 1
  %   h3 = imshow(rgb2gray(I));
  % else
  %   h3 = imagesc(I); colormap gray;
  % end

  h3 = imagesc(I);colormap gray;
  %set(h3, 'AlphaData', alpha*mask);

  if nargin > 4,
    set(h3, 'AlphaData', 1-alpha2*(mask2|mask));
  else
    set(h3, 'AlphaData', 1-alpha*mask);
  end

  hold off;
end %  function	 ctionshowmask()