function [seg,L,seg_vals,seg_lab_vals,seg_edges,seg_img] = make_superpixels(img_loc,para_ms,para_gbis)

% Perform multiple over-segmentations by Mean Shift and FH with
% varying parameters

img = imread(img_loc);
[X,Y,Z] = size(img);
lab_img = colorspace('Lab<-', img);
lab_vals = reshape(lab_img, X*Y, Z);%��ÿһ���ռ���������

%%% do segmentation
for k = 1:para_ms.K
    [~, L{k}, seg{k}, seg_vals{k}, seg_lab_vals{k}, seg_edges{k}] = ...
        msseg(double(img),lab_vals,para_ms.hs{k},para_ms.hr{k},para_ms.M{k});
end

for i = 1:para_gbis.K
    k = i + para_ms.K;
    [seg{k}, L{k}, seg_vals{k}, seg_lab_vals{k}, seg_edges{k}] = ...
        gbis(img,lab_vals,para_gbis.sigma{i},para_gbis.k{i},para_gbis.minsize{i});
end

%%% make mean color image for display
for k = 1:(para_ms.K + para_gbis.K)
    Mimg = zeros(X*Y,Z); 
    for i = 1:size(seg{k},2)
        for j = 1:Z
            Mimg(seg{k}{i},j) = seg_vals{k}(i,j)/255;
        end
    end
    seg_img{k} = reshape(Mimg,[X,Y,Z]);
end

