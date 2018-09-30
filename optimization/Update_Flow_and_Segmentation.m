function [labels, fwarp, bwarp, outlier, layers] = Update_Flow_and_Segmentation(org, labels, fwarp, bwarp, layers, time, beta, alpha, gamma)
k = length(layers);
occlusion_threshold = 0.1;
[h, w] = size( org{1}( :,:,1 ) );
for i = 1:k
    temp = layers{i};
    temp(isnan(temp)) = -1;
    layers{i} = temp;
end
stepsize = 3; %gradient descent step size updating segmentation
ite = 5; %number of iterations
tic;
layer_update_index = k;
for iteration = 1:ite
    disp('Refining optical flow...');
    % only update optical flow once
    [fwarp, bwarp, motion_term, reliability, outlier] = refine_flow(org, layers, labels, fwarp, bwarp, time, occlusion_threshold, ite == 1);
    disp('Computing intensity and smoothness term...');
    [intensity_term, smoothness_term] = intensity_and_smoothness(org, labels, k, stepsize, [20,30] , 4);
    disp('Updating segmentation...');
    parfor i = 1:length(org)
        outside = (labels{i} == -1);
        intensity_updated = zeros(h,w,k);
        for j = 1:k
            intensity_updated(:,:,j) = intensity_term{i}(:,:,j)./reliability{i};
        end
        seg_update = motion_term{i} + iteration*beta*intensity_updated + iteration*alpha*smoothness_term{i};
        % the more iterations done, the more the model relies on intensity
        % and smoothness. Therefore in the first couple of iterations
        % motion contributes more.
        [~, temp] = min(permute(seg_update, [3,1,2]));
        label_updated = permute(temp, [2,3,1]);
        temp_label = labels{i};
        temp_label(temp_label<1)=1;
        candidate= get_edge(edge(temp_label), stepsize);
        labels{i}(candidate) = 0;
        labels{i} = uint8((1 - candidate) .* double(labels{i})) + uint8(candidate .* label_updated);
        labels{i}(outside) = - 1;
        outlier{i}(outside) = 0;
    end
    if layer_update_index >= 1
        disp('Updating 2-D representation...');
        layers = update_template(layers, labels, motion_term, fwarp, bwarp, gamma, layer_update_index);
        layer_update_index = layer_update_index - 1;
    end
    if ite == iteration        
        parfor i = 1:length(org)
            outlier{i} = (org{i}(:,:,1)>=0).*outlier{i};
            new_label{i} = uint8(org{i}(:,:,1)>=0);
            for j = 2:k
                
                temp = remove_holes_and_small_regions(labels{i} == j,round(h*w/500),'weak');
                new_label{i}(logical(temp)) = j;
            end
        end
        labels = new_label;
    end
end
toc;
end

function [fwarp, bwarp, motion_term, reliability, outlier] = refine_flow(org, layers, labels, fwarp, bwarp, time, occlusion_threshold, refine)
mid = round(length(org)/2);
[h, w] = size( org{1}( :,:,1 ) );
[motion_term, reliability, outlier] = deal(cell(length(org),1));
[x, y] = meshgrid(1:w, 1:h);
parfor i = 1:length(org)
    if i >= mid
        b = fwarp(i,:);
        f = bwarp(i,:);
    else
        f = fwarp(i,:);
        b = bwarp(i,:);
    end
    [residual, time_rev]= deal(zeros(h,w,length(layers)));
    for j = 1:length(layers)
        if refine ==1
            [f{j}, b{j}, ~, ~] = warp_matlab(org{i}, layers{j}, (labels{i}==j), f{j}, b{j});
        end
        res = computeresidual(org{i}, layers{j}, f{j});
        res(isnan(res))=100;
        residual(:,:,j) = res;
        w_x=x+f{j}(:,:,1);
        w_y=y+f{j}(:,:,2);
        time_rev(:,:,j) = interp2(time{j}, w_x, w_y);
    end
    if i >= mid
        fwarp(i,:) = b;
        bwarp(i,:) = f;
    else
        fwarp(i,:) = f;
        bwarp(i,:) = b;
    end
    motion_term{i} = residual;
    outlier{i} = permute((min(permute(residual,[3,1,2]))>occlusion_threshold),[2,3,1]);
    time_rev(time_rev<0.1) = 100; % By this setting, the reliability only depends on neighboring layers
    reliability{i} = permute(min(permute(time_rev,[3,1,2])),[2,3,1]);
    reliability{i}(reliability{i}<0.5) = 0.5; % set the lower bound of the reliability
end
end


function [intensity_term, smoothness_term] = intensity_and_smoothness(org, labels, k, edge_width, size_neighborhood, delta)
[h, w] = size( org{1}( :,:,1 ) );
[intensity_term, smoothness_term] = deal(cell(length(org),1));
h_n = size_neighborhood(1);
w_n = size_neighborhood(2);
parfor i = 1:length(org)
    [temp1, temp2]= deal( zeros(h,w,k));
    for j = 1:k
        % Higher intensity_likelihood, more likely to be outside the region
        region = (labels{i}==j);
        boundary = get_edge(edge(region), edge_width, 'index');
        temp1(:,:,j) =  mexIntensity(org{i}, region, boundary, 0.002, h_n, w_n) ...
            + mexIntensity(org{i}, region, boundary, 0.002, h_n*2, w_n*2) ...
            + mexIntensity(org{i}, region, boundary, 0.004, h_n, w_n) ...
            + mexIntensity(org{i}, region, boundary, 0.006, h_n, w_n); % computing local histogram
        temp2(:,:,j) = -imgaussfilt(double(labels{i}==j), delta);
    end
    intensity_term{i} = temp1;
    smoothness_term{i} = temp2;
end
end

function template = update_template(template, label, residual, f_flow, b_flow, area_weight, layer_index)
[h, w] = size(label{1});
[x, y] = meshgrid(1:w, 1:h);
mid = round(length(label)/2);
saving = cell(length(label),1);
parfor i = 1:length(label)
    res = residual{i};
    if i >= mid
        b = f_flow(i,:);
    else
        b = b_flow(i,:);
    end
    base = zeros(h,w);
    for j = 1:length(template)
        temp = res(:,:,j);
        base = base + temp.*(label{i}==j);
        temp(label{i}==j) = inf;
        res(:,:,j) = temp;
    end
    temp = min(permute(res, [3,1,2]));
    gap = permute(temp, [2,3,1]) - base;
    temp_save = zeros(h,w,length(template));
    for j = 1:length(template)
        w_x=x+b{j}(:,:,1);
        w_y=y+b{j}(:,:,2);
        temp = gap.*(label{i}==j);
        temp_save(:,:,j) = interp2(temp, w_x, w_y);
        temp_save(isnan(temp_save)) = 1;
    end
    saving{i} = temp_save;
end
gap = zeros(h,w,length(template));
for i = 1:length(label)
    gap = gap + saving{i};
end
if nargin == 7
    j = layer_index;
    for i = 1:3
        temp2 = template{j}(:,:,i);
        temp2(gap(:,:,j)<area_weight) = -1;
        template{j}(:,:,i) = temp2;
    end
else
    for j = 1:length(template)
        for i = 1:3
            temp2 = template{j}(:,:,i);
            temp2(gap(:,:,j)<area_weight) = -1;
            template{j}(:,:,i) = temp2;
        end
    end
end
end



function boundary = get_edge(input, n, options)
% nargin == 3: output the index of edge
if nargin < 2
    n=1;
end
[h, w]=size(input);
aaa=toeplitz([ones(1,n+1), zeros(1,h-n-1)]);
bbb=toeplitz([ones(1,n+1), zeros(1,w-n-1)]);
result=logical(aaa*double(input)*bbb);

if nargin == 3
    [boundary(:,1), temp] = find(result);
    boundary(:,2) = temp;
else
    boundary = result;
end
end