function [result, occlusion] = computeresidual(im1,im2,flow,occ_t)

[h, w]=size(flow(:,:,1));
[x, y] = meshgrid(1:w, 1:h);
w_x=x+flow(:,:,1);
w_y=y+flow(:,:,2);
%%gray scale
% im1=im2double(rgb2gray(im1));
% im2=im2double(rgb2gray(im2));
% im2=interp2(im2, x+flow(:,:,1), y+flow(:,:,2));
% result=abs(im1-im2);
%%%%rgb scale
result=zeros(h,w);
for i=1:3
    temp=interp2(double(im2(:,:,i)), w_x, w_y);
    result=result+(temp-double(im1(:,:,i))).^2;
end
% temp=result;
% temp(isnan(temp))=0;
% background=mean(mean(temp));
% result(isnan(result))=background;
% result=sqrt(result); %l2 norm

occlusion = zeros(h,w);

if nargin == 4

% occlusion=zeros(h,w);
occlusion(result>occ_t)=1;
occlusion=logical(occlusion);
end
end