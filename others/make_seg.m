function [seg seg_vals seg_lab_vals seg_edges] = make_seg(img,L,lab_vals,full)

%%% input
% Img - image matrix
% L - label map
% lab_vals - [X*Y, Z] color feature vector in Lab color space
%%% output
% seg - pixel indexes in each superpixel
% seg_vals - mean RGB color in each superpixel
% seg_lab_vals - mean Lab color in each superpixel


img = double(img);
[X,Y,Z] = size(img); 
nseg = max(L(:)); 
vals = reshape(img,X*Y,Z);

if full == 1
    [x y] = meshgrid(1:nseg,1:nseg);
    seg_edges = [x(:) y(:)];
else
    [~,edges] = lattice(X,Y,0);
    d_edges = edges(L(edges(:,1))~=L(edges(:,2)),:);
    all_seg_edges = [L(d_edges(:,1)) L(d_edges(:,2))]; all_seg_edges = sort(all_seg_edges,2);
    
    tmp = zeros(nseg,nseg);
    tmp(nseg*(all_seg_edges(:,1)-1)+all_seg_edges(:,2)) = 1;
    [edges_x edges_y] = find(tmp==1); seg_edges = [edges_x edges_y];
end

seg_vals = zeros(nseg,Z);
seg_lab_vals = zeros(nseg,size(lab_vals,2));
for i=1:nseg
    seg{i} = find(L(:)==i);
    seg_vals(i,:) = mean(vals(seg{i},:));
    seg_lab_vals(i,:) = mean(lab_vals(seg{i},:));
end

