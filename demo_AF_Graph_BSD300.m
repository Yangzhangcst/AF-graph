close all;clear;

addpath 'others'
addpath 'evals'
addpath 'KSC'
addpath 'SSC'
addpath 'functions'

%% set parameters for bipartite graph
para.alpha = 0.001; % affinity between pixels and superpixels
para.beta  =  20;   % scale factor in superpixel affinity
para.nb = 1; % number of neighbors for superpixels
para.rho = 1;
para.Nimgs = 300; % number of images in BSDS300
para.alphak = 1e-5;
para.betak = 1e-6;
para.save = 0; % save results

%% read image
bsdsRoot='BSD';
load_file='bsd_300_feat';
outputpath = 'results';
fid = fopen(fullfile('Nsegs_bsd300_l.txt'),'r');
[BSDS_INFO] = fscanf(fid,'%d %d \n',[2,para.Nimgs]);
fclose(fid);

Neg_all = zeros(para.Nimgs,1);
PRI_all = zeros(para.Nimgs,1);
VoI_all = zeros(para.Nimgs,1);
GCE_all = zeros(para.Nimgs,1);
BDE_all = zeros(para.Nimgs,1);

%%
for idxI = 1:para.Nimgs
    % read number of segments
    tic; Nseg = BSDS_INFO(2,idxI);   
    out_path= fullfile(outputpath,'BSDS300');
    if ~exist(out_path,'dir'), mkdir(out_path);  end
    
    %% locate image
    img_name = int2str(BSDS_INFO(1,idxI));
    img_loc = fullfile(bsdsRoot,'images','test',[img_name,'.jpg']);
    if ~exist(img_loc,'file')
        img_loc = fullfile(bsdsRoot,'images','train',[img_name,'.jpg']);
    end
    img = im2double(imread(img_loc)); [X,Y,~] = size(img);
    load_name = fullfile(load_file,[img_name '.mat']); load(load_name)
    
    %% Construct graph
    Np = X*Y;   Nsp = 0;
    for k = 1:length(seg)
        Nsp = Nsp + size(seg{k},2);
    end

    W_Y = sparse(Nsp,Nsp); edgesXY = [];  j = 1;
    for k = 1:length(seg) % for each over-segmentation
        feature = feat{k}.mlab;
        % subspace-preserving representation
        tmp1 = dimReduction_PCA(feature,0); % dimension reduction (if necessary) 
        tmp1 = cnormalize_inplace(tmp1);
        R = SP_mat_func(tmp1', 3, 1e-6); % second parameter is sparsity
        R(1:length(feature)+1:end) = 0;
        A = abs(R) + abs(R)';
        % affinity node selection (global nodes)
        nGCluster(idxI,k) = APclustering(feat{k});
        index_tmp = SpectralClustering(A, nGCluster(idxI,k), 'Eig_Solver', 'eigs');
        % superpixel division
        local_nodes  = find(index_tmp == mode(index_tmp));
        global_nodes = find(index_tmp ~= mode(index_tmp));
        
        feature(:,all(feature == 0, 1))=[];
        [fm,fn] = size(feature);
        feature=(feature-repmat(mean(feature),fm,1))./repmat(std(feature),fm,1); 

        % adjacent graph
        w = makeweights(seg_edges{k},feature,para.beta);
        W_local = adjacency(seg_edges{k},w);
        W = W_local;
       
        % KSC_graph
        W_KSC = KSC_graph(feature,para);
        
        % update graph
        W = assignGraphValue(W,W_KSC,global_nodes);
        W = sparse(W);
        
        Nk = size(seg{k},2); % number of superpixels in over-segmentation k
        W_Y(j:j+Nk-1,j:j+Nk-1) = prune_knn(W,para.nb);
        % affinities between pixels and superpixels
        for i = 1:Nk
            idxp = seg{k}{i}; % pixel indices in superpixel i
            Nki = length(idxp);
            idxsp = j + zeros(Nki,1);
            edgesXY = [edgesXY; [idxp, idxsp]];
            j = j + 1;
        end
    end
    W_XY = sparse(edgesXY(:,1),edgesXY(:,2),para.alpha,Np,Nsp);
    % affinity between a superpixel and itself is set to be the maximum 1.
    W_Y(1:Nsp+1:end) = 1;  B = [W_XY;W_Y];
 
    %% Graph patition
    label_img = Tcut(B,Nseg,[X,Y]); clear B; ti = toc;

    % save segmentation
    view_segmentation(img,label_img(:),out_path,img_name,para.save);

    % evaluate segmentation
    [gt_imgs, gt_cnt] = view_gt_segmentation(bsdsRoot,img,out_path,img_name,para.save);
    out_vals = eval_segmentation(label_img,gt_imgs); 
    fprintf('%6s: %2d %9.6f, %9.6f, %9.6f, %9.6f %.4fs\n', img_name, Nseg,...
        out_vals.PRI, out_vals.VoI, out_vals.GCE, out_vals.BDE, ti);
    PRI_all(idxI) = out_vals.PRI;
    VoI_all(idxI) = out_vals.VoI;
    GCE_all(idxI) = out_vals.GCE;
    BDE_all(idxI) = out_vals.BDE;
end
%%
fprintf('Mean: %14.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
fid_out = fopen(fullfile(outputpath,'BSDS300','evaluation.txt'),'w');
for idxI=1:para.Nimgs
    fprintf(fid_out,'%6d %9.6f, %9.6f, %9.6f, %9.6f \n', BSDS_INFO(1,idxI),...
        PRI_all(idxI), VoI_all(idxI), GCE_all(idxI), BDE_all(idxI));
end
fprintf(fid_out,'Mean: %10.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
fclose(fid_out);
