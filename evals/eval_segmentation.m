function out_vals = eval_segmentation(label_img,gt_imgs)

% Four metrics are used: PRI, VoI, GCE, and BDE

out_vals.PRI = 0; out_vals.VoI = 0; out_vals.GCE = 0; out_vals.BDE = 0;
for i=1:size(gt_imgs,2)
    out_vals.BDE = out_vals.BDE + compare_image_boundary_error(label_img,gt_imgs{i});
    [curRI,curGCE,curVOI] = compare_segmentations(label_img,gt_imgs{i});
    out_vals.PRI = out_vals.PRI + curRI;
    out_vals.VoI = out_vals.VoI + curVOI;
    out_vals.GCE = out_vals.GCE + curGCE;
end
out_vals.PRI = out_vals.PRI/size(gt_imgs,2);
out_vals.VoI = out_vals.VoI/size(gt_imgs,2);
out_vals.GCE = out_vals.GCE/size(gt_imgs,2);
out_vals.BDE = out_vals.BDE/size(gt_imgs,2);