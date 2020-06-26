function A = prune_knn(W,nb)

% W - affinity matrix
% nb - the number of nearest neighbors 
% A - symmetric (weighted) knn graph

n = size(W,1);

[~,idx] = sort(W,'descend');

idxR = idx(1:nb,:);
idxC = ones(nb,1)*[1:n];
A = sparse(idxR(:),idxC(:),1,n,n);
A = double(A|A');
A = A.*W;

