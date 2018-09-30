function im = solid_color_image(imsize, color)
% 	SOLID_COLOR_IMAGE   Short description
% 		IM = SOLID_COLOR_IMAGE(SIZE, COLOR)
% 

switch lower(color)
	case 'red'
		im = cat(3, ones(imsize), zeros(imsize), zeros(imsize));
	case 'green'
		im = cat(3, ones(imsize), ones(imsize), zeros(imsize));
	case 'blue'
		im = cat(3, zeros(imsize), zeros(imsize), ones(imsize));
	case 'yellow'
		im = cat(3, ones(imsize), ones(imsize), zeros(imsize));
	otherwise
		im = cat(3, zeros(imsize), zeros(imsize), zeros(imsize));
end

end %  function