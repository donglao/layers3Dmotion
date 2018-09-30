function labels = get_initialized_label( fwd, bwd )
addpath( genpath( './Initialization/3rdparty/edges-master' ));
model=load('./Initialization/3rdparty/edges-master/models/forest/modelBsds'); 
model=model.model;
model.opts.nms=-1; model.opts.nThreads=4;
model.opts.multiscale=0; model.opts.sharpen=2;
[h,w,~] = size(fwd{1});
%% set up opts for spDetect (see spDetect.m)
T=0.4;
opts = spDetect;
opts.nThreads = 4;  % number of computation threads
opts.k = 512;       % controls scale of superpixels (big k -> big sp)
opts.alpha = .5;    % relative importance of regularity versus data terms
opts.beta = .9;     % relative importance of edge versus color terms
opts.merge = 0;     % set to small value to merge nearby superpixels at end
parfor i = 2 : length(fwd) 
     I = flowToColor( fwd{i}-bwd{i-1} );
            [E,~,~,segs]=edgesDetect(I,model);
            tic, [S,V] = spDetect(I,E,opts); toc
            tic, [~,~,U]=spAffinities(S,E,segs,opts.nThreads); toc
            U(U<T)=0;
            labels{i} = bwlabel(U <= T); 
            labels{i}( isnan( fwd{i}(:,:,1) ) ) = 0;
            list=unique(labels{i});
            frequency = histc(labels{i}(:),list);
            frequency (list == 0) = 0;
            [ a , b ] = max(frequency);
            labels{i}( labels{i} == list(b) ) = 0;
            for j = 1 : length(list)
                if frequency(j) < 100
                     labels{i}( labels{i} == j) = 0;
                end
            end
            labels{i}( labels{i} > 0.5 ) = 1;
            labels{i} = imgaussfilt(labels{i});
            labels{i} = double( logical( labels{i} ) );
            labels{i} = bwareaopen(labels{i} , 100);
end
[labels{1}, labels{length(fwd)+1}] = deal(false(h,w));
end