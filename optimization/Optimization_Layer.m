function [l, f_flow, b_flow, layers, time] = Optimization_Layer(o,l,f_flow,b_flow)
scale = 1.1;
[h,w,~] = size(o{1});
h_mod = round(h*scale);
w_mod = round(w*scale);
n = length(o);
[~, k] = size(f_flow);
for i = n:-1:1
    o{i} = im2double(o{i});
    outlier{i} = zeros(h_mod,w_mod);
    oo{i} = resize_simple(o{i},h_mod,w_mod,-1);
    ll{i} = resize_simple(l{i},h_mod,w_mod,-1);
end
for ite = 1:2
    if ite == 1
        [ll, f_flow, b_flow, layers, time, outlier] = update_scale(o,l,f_flow,b_flow,scale,outlier);
    else
        [ll, f_flow, b_flow, layers, time, outlier] = update_scale(oo,ll,f_flow,b_flow,1,outlier);
    end
end
a = round((h_mod - h)/2)+1;
b = round((w_mod - w)/2)+1;
for i = 1:n
    l{i} = ll{i}(a:a+h-1,b:b+w-1);
    outlier{i} = outlier{i}(a:a+h-1,b:b+w-1);
    for j = 1:k
        b_flow{i,j} = b_flow{i,j}(a:a+h-1,b:b+w-1,:);
        f_flow{i,j} = f_flow{i,j}(a:a+h-1,b:b+w-1,:);
    end
end
end

function [l, f_flow, b_flow, layers, time, outlier] = update_scale(o,l,f_flow,b_flow,scale, outlier)
[h,w,~] = size(o{1});
n = length(o);
[~, k] = size(f_flow);
if scale > 1
    h_mod = round(h*scale);
    w_mod = round(w*scale);
    parfor i = 1:n
        o{i} = resize_simple(o{i},h_mod,w_mod,-1);
        l{i} = resize_simple(l{i},h_mod,w_mod,-1);
        for j = 1:k
            b_flow{i,j} = resize_fill(b_flow{i,j},h_mod,w_mod);
            f_flow{i,j} = resize_fill(f_flow{i,j},h_mod,w_mod);
        end
    end
end
[layers, time, f_flow, b_flow] = Optimize_Representations(o, l, f_flow, b_flow, outlier);
[l, f_flow, b_flow, outlier, layers] = Update_Flow_and_Segmentation(o, l, f_flow, b_flow, layers, time, 10, 5, 0.1);
end