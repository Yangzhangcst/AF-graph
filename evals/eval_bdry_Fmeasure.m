function [FM] = eval_bdry_Fmeasure(labels,gt_imgs,maxDist)
%
% INPUT
%	inFile  : Can be one of the following:
%             - a soft or hard boundary map in image format.
%             - a collection of segmentations in a cell 'segs' stored in a mat file
%             - an ultrametric contour map in 'doubleSize' format, 'ucm2'
%               stored in a mat file with values in [0 1].
%
%	gtFile	: File containing a cell of ground truth boundaries
%   prFile  : Temporary output for this image.
%	nthresh	: Number of points in PR curve.
%   MaxDist : For computing Precision / Recall.
%   thinpb  : option to apply morphological thinning on segmentation
%             boundaries.
%
% OUTPUT
%	thresh		Vector of threshold values.
%	cntR,sumR	Ratio gives recall.
%	cntP,sumP	Ratio gives precision.
%
%   FM = 2*P*R/(P+R) ;
%

maxDist = 0.0075;

bmap = logical(seg2bdry(labels,'imageSize'));

for i = 1: size(gt_imgs,2)
    gt_bdry{i} = logical(seg2bdry(gt_imgs{i},'imageSize'));
end

% thin the thresholded pb to make sure boundaries are standard thickness

bmap = double(bwmorph(bmap, 'thin', inf));    % OJO

% accumulate machine matches, since the machine pixels are
% allowed to match with any segmentation
accP = zeros(size(bmap));
cntR = 0;
sumR = 0;

% compare to each seg in turn
for i = 1:size(gt_bdry,2),
    % compute the correspondence
    [match1,match2] = correspondPixels(bmap, double(gt_bdry{i}), maxDist);
    % accumulate machine matches
    accP = accP | match1;
    % compute recall
    sumR = sumR + sum(gt_bdry{i}(:));
    cntR = cntR + sum(match2(:)>0);
end

% compute precision
sumP = sum(bmap(:));
cntP = sum(accP(:));

R = sumR/cntR;
P = sumP/cntP;

FM = 2*P*R/(P+R);



