function [final_label, layers, f_flow, b_flow, time] = Layer_Segmentation(org, labels)
mid = round(length(org)/2);
k = max(max(max(labels{mid}))+1,2);
[h, w] = size( org{1}( :,:,1 ) );
[f_flow, b_flow] = deal(cell(1,k));
for i = 1:k
    [f_flow{1,i}, b_flow{1,i}] = deal(zeros(h,w,2));
end
l{1} = labels{mid}+1;
% initialize layers
temp = im2double(org{mid});
for j = 1:k
    temp2 = double(l{1}==j);
    temp2(temp2==0) = nan;
    for i = 1:3
        layers{j}(:,:,i) = temp(:,:,i).*temp2;
    end
end

for expansion = 1:mid-1
    first = max(1,mid-expansion);
    last = min(length(org),mid+expansion);
    [temp_f, temp_b] = deal(cell(last-first+1,1));
    label_expand{last-first+1} = zeros(h,w);
    label_expand{1} = zeros(h,w);
    for i = 2:(last-first)
        o{i} = org{i+first-1};
        for j = 1:k
            temp_f(i,j) = f_flow(i-1,j);
            temp_b(i,j) = b_flow(i-1,j);
        end
        label_expand{i} = uint8(l{i-1});
    end
    disp('Propagating Initialized region...');
    tic
    temp_org_a = org{first+1};
    temp_org_b = org{first};
    temp_org_c = org{last-1};
    temp_org_d = org{last};
    temp_label_a = label_expand{2};
    temp_label_b = label_expand{last-first};
    temp_bflow_a = b_flow(1,:);
    temp_fflow_a = f_flow(1,:);
    temp_bflow_b = b_flow(last-first-1,:);
    temp_fflow_b = f_flow(last-first-1,:);
    parfor flow_index = 1:k*2
        i = mod(flow_index, k)+1; % i: layer number
        if flow_index <= k
            [b1, f1, temp] = warp_matlab(temp_org_a,temp_org_b, (temp_label_a==i));  %testing without initialization
            label_temp{flow_index} = (temp==1);
            flow_temp1{flow_index} = compose_flow_single( temp_bflow_a{i}, b1);
            flow_temp2{flow_index} = compose_flow_single(f1, temp_fflow_a{i});
        else
            [fe, be, temp] = warp_matlab(temp_org_c,temp_org_d, (temp_label_b==i));
            label_temp{flow_index} = (temp==1);
            flow_temp1{flow_index} = compose_flow_single(be, temp_bflow_b{i});
            flow_temp2{flow_index} = compose_flow_single(temp_fflow_b{i}, fe);
        end
    end
    for flow_index = 1:k*2
        i = mod(flow_index, k)+1;
        if flow_index <= k
            temp = label_temp{flow_index};
            label_expand{1}(temp==1) = i;
            temp_b{1,i} = flow_temp1{flow_index};
            temp_f{1,i} = flow_temp2{flow_index};
        else
            temp = label_temp{flow_index};
            label_expand{last-first+1}(temp==1) = i;
            temp_b{last-first+1,i} = flow_temp1{flow_index};
            temp_f{last-first+1,i} = flow_temp2{flow_index};
        end
    end
    toc
    % trying better initialization
    for i = 2:k
        if fmeasure((labels{first}==i-1),(label_expand{1}==i))>0.5
            label_expand{1}(labels{first}==i-1) = i;
        end
        if fmeasure((labels{last}==i-1),(label_expand{expansion*2+1}==i))>0.5
            label_expand{expansion*2+1}(labels{last}==i-1) = i;
        end
    end
    %
    o{1}=org{first};
    label_expand{1} = uint8(label_expand{1});
    label_expand{expansion*2+1} = uint8(label_expand{expansion*2+1});
    o{last - first + 1} = org{last};
    [l, f_flow, b_flow, layers, time]= Optimization_Layer(o, label_expand, temp_f, temp_b);
    
    %     figure; imagesc(l{round(length(l)/2)});
    k = length(layers);
end
for i = 1:length(l)
    final_label{i} = l{i}-1;
end
end

function f = fmeasure(r1, r2)
r1 = logical(r1);
r2 = logical(r2);
t = sum(sum(r1.*r2));
a = t/sum(sum(r1));
b = t/sum(sum(r2));
if a+b >0
    f = 2*a*b/(a+b);
else f = 0;
end
end