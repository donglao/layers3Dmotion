% Extending Layered Models to 3D Motion demo program
%
% References:
% -----------
% Dong Lao & Ganesh Sundaramoorthi "Extending Layered Models to 3D Motion" 
% 15th European Conference on Computer Vision, 2018 
% 
% Authors: Dong Lao and Ganesh Sundaramoorthi, King Abdullah University of
% Science and Technology
% Contact: dong.lao@kaust.edu.sa
% $Date: 8/25/2018$
% 
% Copyright 2015-2018, King Abdullah University of Science and Technology
% (KAUST), Saudi Arabia
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
% supporting documentation, and that the name of the author and KAUST
% not be used in advertising or publicity pertaining to distribution 
% of the software without specific, written prior permission.        
% 
% THE AUTHOR AND DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, 
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
% ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR 
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
% RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
% CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
% CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.    



clear all;
clc;
% run mex_compile.m
load bear_demo.mat;
addpath('./Initialization')
load bear_demo.mat;
[forwardflow, backwardflow] = get_optical_flow_for_initialization( org );
save bear_demo2.mat forwardflow backwardflow;
% load bear_demo2.mat forwardflow backwardflow;
labels = get_initialized_label( forwardflow, backwardflow );

addpath(genpath('./toolbox'))
addpath(genpath('./optimization'))
[labels_updated, layers] =  Layer_Segmentation(org, labels);
% figure;
% for i = 1:length(labels_updated)
%     subplot(3,length(labels_updated),i);
%     imagesc(org{i});
%     subplot(3,length(labels_updated),length(labels_updated)+i);
%     imagesc(labels{i});
%     subplot(3,length(labels_updated),length(labels_updated)*2+i);
%     imagesc(labels_updated{i});
% end

for i = 1:length(labels_updated)
    figure;
    imagesc(org{i});
    hold on;
    contour(labels_updated{i});
end
figure;
subplot(2,1,1);
imagesc(layers{1});
subplot(2,1,2);
imagesc(layers{2});

    
    