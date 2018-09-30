function [layer_raw, time, nan_map, jacobian] = layer_interp(img,flow,mask)
[h, w] = size(img(:,:,1));
[x, y] = meshgrid(1:w, 1:h);
w_x=x+flow(:,:,1);
w_y=y+flow(:,:,2);
layer_raw=zeros(h,w,3);
img(isnan(img))=0;
% img(img<0)=0;
% mask(mask==0)= nan;
for i=3:-1:1
    layer_raw(:,:,i)=interp2(img(:,:,i).*mask, w_x, w_y);
end
time=interp2(double(mask), w_x, w_y);
time(isnan(time)) = 0;
temp = isnan(layer_raw);
layer_raw(temp)=0;
nan_map = temp(:,:,1);
[A, B] = gradient(w_x);
[C, D] = gradient(w_y);
jacobian = abs(A.*D - B.*C);
jacobian(jacobian>= 10) = 10;
end