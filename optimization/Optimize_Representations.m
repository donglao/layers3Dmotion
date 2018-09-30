function  [layers, time, fwarp, bwarp] = Optimize_Representations(org, labels, fwarp, bwarp, outlier, scale)
if nargin < 6
    scale = 1;
end
[h,w] = size(labels{1});
mid = round(length(org)/2);
[~,num_of_layers] = size(fwarp);
[l, flow_forward, flow_backward] = deal(cell(length(org),1));
[layers, time] = deal(cell(num_of_layers, 1));
disp('Optimizing 2D representation appearance...');
tic
for i = 1:num_of_layers
    for j = 1:length(org)
    if nargin < 5
        bad = (labels{j}==0);
    else
        bad = outlier{j};
    end
        l{j} = logical((labels{j} == i) + bad.*get_edge(edge(labels{j}==i),3));
        flow_forward{j} = fwarp{j,i};
        flow_backward{j} = bwarp{j,i};
    end
    [layers{i}, time{i}] = layerbuild(org, l, flow_forward, flow_backward, scale);
    [a,b] = find(layers{i}(:,:,1)>0); % a row, b column
    row_shift = round((h+1-max(a)-min(a))/2);
    column_shift = round((w+1-max(b)-min(b))/2);
    if isempty(row_shift)
        row_shift = 0;
    end
    if isempty(column_shift)
        column_shift = 0;
    end
for k = 1:3
    layers{i}(:,:,k) = circshift(layers{i}(:,:,k),[row_shift, column_shift]); % put the representation in the middle
end
for j = 1:length(org)
    if j>=mid
        fwarp{j,i}(:,:,1) = fwarp{j,i}(:,:,1) - column_shift;
        fwarp{j,i}(:,:,2) = fwarp{j,i}(:,:,2) - row_shift;
        bwarp{j,i}(:,:,1) = bwarp{j,i}(:,:,1) + column_shift;
        bwarp{j,i}(:,:,2) = bwarp{j,i}(:,:,2) + row_shift;
    else
        fwarp{j,i}(:,:,1) = fwarp{j,i}(:,:,1) + column_shift;
        fwarp{j,i}(:,:,2) = fwarp{j,i}(:,:,2) + row_shift;
        bwarp{j,i}(:,:,1) = bwarp{j,i}(:,:,1) - column_shift;
        bwarp{j,i}(:,:,2) = bwarp{j,i}(:,:,2) - row_shift;
    end
end
end
toc
end

function boundary = get_edge(input, n)
if nargin < 2
    n=1;
end
[h, w]=size(input);
a=toeplitz([ones(1,n+1), zeros(1,h-n-1)]);
b=toeplitz([ones(1,n+1), zeros(1,w-n-1)]);
boundary=logical(a*double(input)*b);
end