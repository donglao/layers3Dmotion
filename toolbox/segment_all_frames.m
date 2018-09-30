function final_label = segment_all_frames(org, labels, forwardflow, backwardflow, interval)
%This part is modified to pick up the right number of interval
max_interval = get_interval(forwardflow);
if nargin < 5
interval = 11;
elseif interval > max_interval
interval = max(max_interval,11);
end
disp(interval);
mid = round(interval/2);
[h,w] = size(org{1}(:,:,1));
labels{length(org)} = zeros(h,w);
for i = 1:length(org)
    if sum(sum(labels{i})) == 0
        labels{i} = zeros(h,w);
    end
end
step = round(interval/3);
group = ceil((length(org)-interval)/step) + 1;
final_label = labels;
for index = 0:group-1
    if index < group-1
        k = step*index+1;
    else k = length(org)-interval + 1;
    end
    [l, o, f, b] = deal(cell(1,interval));
    for i = k:k+interval-1
        if i < length(org)
            o{i-k+1} = org{i};
            f{i-k+1} = forwardflow{i};
            b{i-k+1} = backwardflow{i};
            l{i-k+1} = final_label{i};
        else o{i-k+1} = org{i};
            l{i-k+1} = final_label{i};
        end
    end
    i = round(interval/2) + k - 1;
    l{round(interval/2)} = final_label{i};
    [l, layers] =  Layer_Segmentation(o, l);
%     for i = 1:length(layers)
%         imwrite(layers{i},['layer_',num2str(i),'.png']);
%     end
    parfor i = 1:interval
        temp_l{i} = ones(h,w);
        for j = 2:max(max(l{i}))
            temp = (l{i}==j);
            temp_l{i}(logical(temp)) = j;
        end
    end 
    for i = mid:mid+step
        final_label{i+k-1} = temp_l{i} - 1;
%         imwrite((temp_l{i}-1)/(max(max(l{i}))-1), [num2str(i+k-1),'.png']);
    end
    if index == 0
        for i = 1:mid-1
            final_label{i+k-1} = temp_l{i} - 1;
%             imwrite((temp_l{i}-1)/(max(max(l{i}))-1), [num2str(i+k-1),'.png']);
        end
    end
    if index == group-1
        for i = mid+1:interval
            final_label{i+k-1} = temp_l{i} - 1;
%             imwrite((temp_l{i}-1)/(max(max(l{i}))-1), [num2str(i+k-1),'.png']);
        end
    end
end
end

function interval = get_interval(forwardflow)
[h,w,~] = size(forwardflow{1});
h_movement = 0;
w_movement = 0;
for i = 1:length(forwardflow)
    h_movement = h_movement + abs(mean(mean(forwardflow{i}(:,:,2))));
    w_movement = w_movement + abs(mean(mean(forwardflow{i}(:,:,1))));
end
split = max(4*h_movement/h , 4*w_movement/w);
if split<1
    split = 1;
end
interval = round((length(forwardflow)+1)/split);
end