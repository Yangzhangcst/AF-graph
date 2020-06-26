function labels = APclustering(feat)
% % 
histogram = feat.mlab;
[width,height,bands] = size(histogram);
H_total = zeros(1,width)';
histogram1 = zeros(width,1,bands);
for i = 1:bands         
    temp = histogram(:,:,i);
    histogram1(:,:,i) = mean(temp,2);         
    %%%% Use 3D Histogram if more than one slice is given%%%%%
    H_total = H_total + histogram1(:,:,i);  
end
% % 
% histogram1 = histogram_smoothing(H_total(H_total>0), 20, 20);
idx = 1:length(H_total);
histogram1 = [idx' H_total];
% % 
binary = affinity_calculation(histogram1, 3, 1);
labels = length(binary(2:end));
if labels==0, labels=1; end
if labels>10, labels=10; end