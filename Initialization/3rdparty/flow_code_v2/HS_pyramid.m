[h,w] = size(img1(:,:,1));
levels = fix(log2(min(h,w)));
flow = zeros(h,w);
for k = levels:-1:0
    h_m = round(h/(2^k));
    w_m = round(w/(2^k));
    img1_m = imresize(img1, [h_m, w_m]);
    img2_m = imresize(img2, [h_m, w_m]);
    flow = imresize(flow, [h_m, w_m]);
    if k == levels
        flow = HS(img1_m,img2_m,100,100);
    else 
        flow = HS(img1_m,img2_m,100,100,flow(:,:,1),flow(:,:,2));
    end
end
