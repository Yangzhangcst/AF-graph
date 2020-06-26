function [gt_imgs gt_cnt] = view_gt_segmentation(bsdsRoot,img,out_path,only_name,save)

gt_imgs = readSegs(bsdsRoot,'color',str2num(only_name));
gt_path = fullfile(out_path, only_name, 'gt'); 
if ~exist(gt_path), mkdir(gt_path);end

gt_cnt = [];
for i=1:size(gt_imgs,2)
    if save == 1
        [imgMasks,segOutline,imgMarkup]=segoutput(img,gt_imgs{i});
        imwrite(imgMarkup,fullfile(gt_path, [only_name, '_', int2str(i), '.bmp'])); 
    end
    gt_cnt(i) = max(gt_imgs{i}(:));
end
