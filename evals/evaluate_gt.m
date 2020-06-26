function [evals num_err_pix num_band_pix num_all_pix] = evaluate_gt(mask, ref_name, seg_name)

trimap = imread(ref_name);
seg = imread(seg_name);

[h,w] = size(trimap);
if h ~= size(mask,1)
    mask = imresize(mask,[h,w]);
    mask(find(mask>=1.5)) = 2;
    mask(find(mask< 1.5)) = 1;
end;

band_idx = find(trimap(:) == 128);

num_band_pix = size(band_idx,1);

num_err_pix = 0;
num_all_pix = 0;
for i=1:num_band_pix
    if mask(band_idx(i))==1 && seg(band_idx(i))==255
        num_err_pix = num_err_pix + 1;
    elseif mask(band_idx(i))==2 && seg(band_idx(i))==0
        num_err_pix = num_err_pix + 1;
    end;
    if seg(band_idx(i))==255 | seg(band_idx(i))==0
        num_all_pix = num_all_pix + 1;
    end;
end;

evals = num_err_pix/num_all_pix;