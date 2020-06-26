function  M = clusteringcenter(X, K, seed)

if ~exist('seed', 'var'), seed = 1; end

[D,N] = size(X);

M = zeros(D, K);
Dist = zeros(N, K);
M(:, 1) = X(:,seed);
Dist(:, 1) = sum((X - repmat(M(:, 1), 1, N)).^2, 1)';
for ii = 2:K
  % maximum, minimum dist
  mindist = min(Dist(:,1:ii-1), [], 2);
  [~, jj] = max(mindist);
  M(:, ii) = X(:, jj);
  Dist(:, ii) = sum((X - repmat(M(:, ii), 1, N)).^2, 1)';
end