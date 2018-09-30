function [representation, time_final] = layerbuild(org, labels, flow_forward, flow_backward, expansion)
if expansion<1
    error('Expansion rate must be larger or equal than 1...');
end
n = length(org);
mid = round(n/2);
[h,w] = size(org{mid}(:,:,1));
h_mod = round(h*expansion);
w_mod = round(w*expansion);

flow=flow_forward;
for i = 1:mid-1
    flow{i}=flow_backward{i};
end
parfor i = 1:n
    [img, mask] = img_resize(org{i},labels{i},h_mod,w_mod);
    motion = flow_resize(flow{i},h_mod,w_mod);
    [layer_raw{i}, time{i}, nan_map{i}, Jacobian{i}] = layer_interp(img,motion,mask);
end
representation = zeros(h_mod,w_mod,3);
time_final = zeros(h_mod,w_mod);
null_value_detector = zeros(h_mod, w_mod);
for i = 1:n
    temp = layer_raw{i};
    for j = 1:3
    temp(:,:,j) = temp(:,:,j).*Jacobian{i};
    end
    representation = representation + temp;
    time_final = time_final + time{i}.*Jacobian{i};
    null_value_detector = null_value_detector.*nan_map{i};
end
for i = 1:3
    temp = representation(:,:,i)./time_final;
    temp(logical(null_value_detector)) = -1;
    representation(:,:,i) = temp;
end
end

function [output, mask] = img_resize(org,label,h_mod,w_mod)
[h,w]=size(org(:,:,1));
label = double(label);
org = im2double(org);
a = round((h_mod - h)/2)+1;
b = round((w_mod - w)/2)+1;
mask = zeros(h_mod,w_mod);
mask(a:a+h-1,b:b+w-1) = label;
output = ones(h_mod,w_mod,3)*nan;
label(label==0) = nan;
[temp(:,:,1),temp(:,:,2),temp(:,:,3)] = deal(label);
output(a:a+h-1,b:b+w-1,:) = org.*temp;
end

function output = flow_resize(flow,h_mod,w_mod)
[h,w]=size(flow(:,:,1));
output = ones(h_mod,w_mod,2)*nan;
a = round((h_mod - h)/2)+1;
b = round((w_mod - w)/2)+1;
output(a:a+h-1,b:b+w-1,:) = flow;
output(:,:,2) = inpaint_nans(output(:,:,2),5);
output(:,:,1) = inpaint_nans(output(:,:,1),5);
end