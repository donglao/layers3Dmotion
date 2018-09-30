function uvo = estimate_flow_interface(im1, im2, method, params) 

%ESTIMATE_FLOW_INTERFACE  Optical flow estimation with various methods
%
% Demo program
%     [im1, im2, tu, tv] = read_flow_file('middle-other', 1);
%     uv = estimate_flow_interface(im1, im2, 'classic+nl-fast');
%     [aae stdae aepe] = flowAngErr(tu, tv, uv2(:,:,1), uv2(:,:,2), 0)
%
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
% Contact: dqsun@cs.brown.edu
% $Date: $
% $Revision: $
%
% Copyright 2007-2010, Brown University, Providence, RI. USA
% 
%                          All Rights Reserved
% 
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.     
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.        
%
% For commercial uses contact the Technology Venture Office of Brown University
% 
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.        


% Read in arguments
if nargin < 3
    method = 'classic+nl-fast';
end;

if (~isdeployed)
    addpath(genpath('utils'));
end

% Load default parameters
ope = load_of_method(method);

if nargin > 3
    ope = parse_input_parameter(ope, params);    
end;

% Uncomment this line if Error using ==> \  Out of memory. Type HELP MEMORY for your option.
%ope.solver    = 'pcg';  

if size(im1, 3) > 1
    tmp1 = double(rgb2gray(uint8(im1)));
    tmp2 = double(rgb2gray(uint8(im2)));
    ope.images  = cat(length(size(tmp1))+1, tmp1, tmp2);
else
    
    if isinteger(im1);
        im1 = double(im1);
        im2 = double(im2);
    end;
    ope.images  = cat(length(size(im1))+1, im1, im2);
end;

% Use color for weighted non-local term
if ~isempty(ope.color_images)    
    if size(im1, 3) > 1        
        % Convert to Lab space       
        im1 = RGB2Lab(im1);          
        for j = 1:size(im1, 3);
            im1(:,:,j) = scale_image(im1(:,:,j), 0, 255);
        end;        
    end;    
    ope.color_images   = im1;
end;

% Compute flow field
uv  = compute_flow(ope, zeros([size(im1,1) size(im1,2) 2]));

if nargout == 1
    uvo = uv;
end;

% function uv = compute_flow(this, init, gt)
% %
% %COMPUTE_FLOW   Compute flow field
% %   UV = COMPUTE_FLOW(THIS[, INIT]) computes the flow field UV with
% %   algorithm THIS and the optional initialization INIT.
% %  
% %   This is a member function of the class 'hs_optical_flow'. 
% %
% %   Author: Deqing Sun, Department of Computer Science, Brown University
% %   Contact: dqsun@cs.brown.edu
% %   $Date: 2007-10-30 $
% %   $Revision: $
% %
% % Copyright 2007-2010, Brown University, Providence, RI. USA
% % 
% %                          All Rights Reserved
% % 
% % All commercial use of this software, whether direct or indirect, is
% % strictly prohibited including, without limitation, incorporation into in
% % a commercial product, use in a commercial service, or production of other
% % artifacts for commercial purposes.     
% %
% % Permission to use, copy, modify, and distribute this software and its
% % documentation for research purposes is hereby granted without fee,
% % provided that the above copyright notice appears in all copies and that
% % both that copyright notice and this permission notice appear in
% % supporting documentation, and that the name of the author and Brown
% % University not be used in advertising or publicity pertaining to
% % distribution of the software without specific, written prior permission.        
% %
% % For commercial uses contact the Technology Venture Office of Brown University
% % 
% % THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% % THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% % FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% % BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% % DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% % PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% % ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% % THIS SOFTWARE.  
% 
%   % Frame size
%   sz = [size(this.images, 1), size(this.images, 2)];
% 
%   % If we have no initialization argument, initialize with all zeros
%   if (nargin < 2)
%     uv = zeros([sz, 2]);
%   else
%     uv = init;
%   end
% 
%   % Perform structure texture decomposition to get the texture component
%   if this.texture == true;
%       % Texture constancy
%       
%       if size(this.images, 4) ==1
%           images  = structure_texture_decomposition_rof( this.images); %, 1/8, 100, this.alp);
%       else
%           
%           for i = 1:size(this.images,4)
%               images(:,:,i,:)  = structure_texture_decomposition_rof( squeeze(this.images(:,:,i,:)));
%           end;
%           
%       end;
% 
%   else
%       images  = scale_image(this.images, 0, 255);      
%   end;
% 
%   % Automatic determine pyramid level
%   this.pyramid_levels  =  1 + floor( log(min(size(images, 1), size(images,2))/16) / log(this.pyramid_spacing) );  
%   
%   % Construct image pyramid, using setting in Bruhn et al in  "Lucas/Kanade.." (IJCV2005') page 218
% 
%   factor            = sqrt(2);  % sqrt(3)
%   smooth_sigma      = sqrt(this.pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger   
%   f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);    
% 
%   pyramid_images    = compute_image_pyramid(images, f, this.pyramid_levels, 1/this.pyramid_spacing);
%    
%   % Iterate through all pyramid levels starting at the top
%   for l = this.pyramid_levels:-1:1
%          
%     if this.display
%         disp(['Pyramid level: ', num2str(l)])    
%     end
%     
%     % Scale flow to the current pyramid level
%     uv    =  resample_flow(uv, [size(pyramid_images{l}, 1) size(pyramid_images{l}, 2)]);
%    
%     % Generate copy of algorithm with single pyramid level and the appropriate subsampling
%     small = this;
%     small.pyramid_levels = 1;
%     small.images         = pyramid_images{l};
%     
%     % Run flow method on subsampled images
%     uv = compute_flow_base(small, uv);    
%  
%   end
%         
%   if this.display
%       fprintf('energy of solution \t%3.3e\n', -evaluate_log_posterior(small, uv));
%   end;
%       
%   if ~isempty(this.median_filter_size)
%       uv(:,:,1) = medfilt2(uv(:,:,1), this.median_filter_size, 'symmetric');
%       uv(:,:,2) = medfilt2(uv(:,:,2), this.median_filter_size, 'symmetric');
%   end;
%   
%   function uv = compute_flow_base(this, uv)
% 
% %COMPUTE_FLOW_BASE   Base function for computing flow field
% %   UV = COMPUTE_FLOW_BASE(THIS, INIT) computes the flow field UV with
% %   algorithm THIS and the initialization INIT.
% %  
% %   This is a member function of the class 'hs_optical_flow'. 
% %
% %   Author: Deqing Sun, Department of Computer Science, Brown University
% %   Contact: dqsun@cs.brown.edu
% %   $Date: 2007-11-30 $
% %   $Revision: $
% %
% % Copyright 2007-2010, Brown University, Providence, RI. USA
% % 
% %                          All Rights Reserved
% % 
% % All commercial use of this software, whether direct or indirect, is
% % strictly prohibited including, without limitation, incorporation into in
% % a commercial product, use in a commercial service, or production of other
% % artifacts for commercial purposes.     
% %
% % Permission to use, copy, modify, and distribute this software and its
% % documentation for research purposes is hereby granted without fee,
% % provided that the above copyright notice appears in all copies and that
% % both that copyright notice and this permission notice appear in
% % supporting documentation, and that the name of the author and Brown
% % University not be used in advertising or publicity pertaining to
% % distribution of the software without specific, written prior permission.        
% %
% % For commercial uses contact the Technology Venture Office of Brown University
% % 
% % THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% % THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% % FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% % BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% % DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% % PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% % ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% % THIS SOFTWARE.  
% 
%   % Iterate flow computation
%   for i = 1:this.max_warping_iters             
%    
%     % Compute linear flow operator
%     [A, b, parm, iterative] = ...
%         flow_operator(this, uv);
%     
%     % Invoke the selected linear equation solver
%     switch (lower(this.solver))
%       case 'backslash'
%         x = reshape(A \ b, size(uv));
%       case 'sor'
%         [x, flag, res, n] = sor(A', b, 1.9, this.sor_max_iters, 1E-2, uv(:));
%         x = reshape(x, size(uv));
%         fprintf('%d %d %d  ', flag, res, n);
%       case 'bicgstab'
%         x = reshape(bicgstab(A, b, 1E-3, 200, [], [], uv(:), parm), size(uv));
%       case 'pcg'
%           [x flag] = pcg(A,b, [], 100);  %100           
%       otherwise
%         error('Invalid solver!')
%     end
%     
%     % Print status information
%     if this.display
%         disp(['--Iteration: ', num2str(i), '    (', ...
%              num2str(norm(x(:))), ')'])
%     end;
%     
%     % Terminate iteration early if flow doesn't change substantially
%     if (length(this.lambda) == 1 && norm(x(:)) < 1E-3)
%       break
%     end
%    
%     % If limiting the incremental flow to [-1, 1] is requested, do so
%     if (this.limit_update)
%       x(x > 1)  = 1;
%       x(x < -1) = -1;
%     end
%     
%     uv = uv + x;
%         
%     % Perform median filtering to remove outliers
%     if ~isempty(this.median_filter_size)
%         for m = 1:this.mf_iter; % extensive MF filtering 
%             uv(:,:,1) = medfilt2(uv(:,:,1), this.median_filter_size, 'symmetric');
%             uv(:,:,2) = medfilt2(uv(:,:,2), this.median_filter_size, 'symmetric');
%         end;
%     end;
%     
%     
%     % Terminate early if the flow_operator doesn't require multiple
%     % interations 
% %     if (~iterative)
% %       break;
% %     end
%     
%   end
