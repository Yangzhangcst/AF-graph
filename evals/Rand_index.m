function ri=Rand_index(sampleLabels1,sampleLabels2)

% compare_segmentations
%
%   Computes several simple segmentation benchmarks. Written for use with
%   images, but works for generic segmentation as well (i.e. if the
%   sampleLabels inputs are just lists of labels, rather than rectangular
%   arrays).
%
%   The measures:
%       Rand Index
%
%   The Rand Index can be easily extended to the Probabilistic Rand Index
%   by averaging the result across all human segmentations of a given
%   image: 
%       PRI = 1/K sum_1^K RI( seg, humanSeg_K ).
%   With a little more work, this can also be extended to the Normalized
%   PRI. 
%
%   Inputs:
%       sampleLabels1 - n x m array whose entries are integers between 1
%                       and K1
%       sampleLabels2 - n x m (sample size as sampleLabels1) array whose
%                       entries are integers between 1 and K2 (not
%                       necessarily the same as K1).
%   Outputs:
%       ri  - Rand Index
%       gce - Global Consistency Error
%       vi  - Variation of Information
%
%   NOTE:
%       There are a few formulas here that look less-straightforward (e.g.
%       the log2_quotient function). These exist to handle corner cases
%       where some of the groups are empty, and quantities like 0 *
%       log(0/0) arise...
%
%   Oct. 2006 
%       Questions? John Wright - jnwright@uiuc.edu

[imWidth,imHeight]=size(sampleLabels1);
[imWidth2,imHeight2]=size(sampleLabels2);
N=imWidth*imHeight;
if (imWidth~=imWidth2)||(imHeight~=imHeight2)
    disp( 'Input sizes: ' );
    disp( size(sampleLabels1) );
    disp( size(sampleLabels2) );
    error('Input sizes do not match in compare_segmentations.m');
end;

% make the group indices start at 1
if min(min(sampleLabels1)) < 1
    sampleLabels1 = sampleLabels1 - min(min(sampleLabels1)) + 1;
end
if min(min(sampleLabels2)) < 1
    sampleLabels2 = sampleLabels2 - min(min(sampleLabels2)) + 1;
end

segmentcount1=max(max(sampleLabels1));
segmentcount2=max(max(sampleLabels2));

% compute the count matrix
%  from this we can quickly compute rand index, GCE, VOI, ect...
n=zeros(segmentcount1,segmentcount2);

for i=1:imWidth
    for j=1:imHeight
        u=sampleLabels1(i,j);
        v=sampleLabels2(i,j);
        n(u,v)=n(u,v)+1;
    end;
end;

ri = rand_index(n);

return;


% the performance measures

% the rand index, in [0,1] ... higher => better
%  fast computation is based on equation (2.2) of Rand's paper.
function ri = rand_index(n)
N = sum(sum(n));
n_u=sum(n,2);
n_v=sum(n,1);
N_choose_2=N*(N-1)/2;
ri = 1 - ( sum(n_u .* n_u)/2 + sum(n_v .* n_v)/2 - sum(sum(n.*n)) )/N_choose_2;

% log2_quotient 
%   helper function for computing the mutual information
%   returns a matrix whose ij entry is
%     log2( a_ij / b_ij )   if a_ij, b_ij > 0
%     0                     if a_ij is 0
%     log2( a_ij + 1 )      if a_ij > 0 but b_ij = 0 (this behavior should
%                                         not be encountered in practice!
function lq = log2_quotient( A, B )
lq = log2( (A + ((A==0).*B) + (B==0)) ./ (B + (B==0)) );