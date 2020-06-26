function [R, M] = k_meanspp(X, K, seed)
%KMEANS:  K-means clustering
%  idx = KMEANS(X, K) returns M with K columns, one for each mean.  Each
%      column of X is a datapoint.  K is the number of clusters
%  [idx, mu] = KMEANS(X, K) also returns mu, a row vector, R(i) is the
%      index of the cluster datapoint X(:, i) is assigned to.
%  idx = KMEANS(X,K) returns idx where idx(i) is the index of the cluster
%      that datapoint X(:,i) is assigned to.
%  [idx,mu] = KMEANS(X,K) also returns mu, the K cluster centers.
%
%  KMEANS(X, K, SEED) uses SEED (default 1) to randomise initial assignments.

if ~exist('seed', 'var'), seed = 1; end

[D,N] = size(X);

M = zeros(D, K);
Dist = zeros(N, K);
M(:, 1) = X(:,seed);
Dist(:, 1) = sum((X - repmat(M(:, 1), 1, N)).^2, 1)';
for ii = 2:K
  % maximum, minimum dist
  mindist = min(Dist(:,1:ii-1), [], 2);
  [junk, jj] = max(mindist);
  M(:, ii) = X(:, jj);
  Dist(:, ii) = sum((X - repmat(M(:, ii), 1, N)).^2, 1)';
  index = Roulettemethod(Dist(:,ii));
  M(:, ii) = X(:,index);
end

X2 = sum(X.^2,1)';
converged = 0;
iter = 0;
R = zeros(N, 1);
while ~converged && iter<100
  distance = repmat(X2,1,K) - 2 * X' * M + repmat(sum(M.^2, 1), N, 1);
  [junk, newR] = min(distance, [], 2);
  if norm(R-newR) == 0
    converged = 1;
  else
    R = newR;
  end
  total = 0;
  for ii = 1:K
    ix = find(R == ii);
    M(:, ii) = mean(X(:, ix), 2);
    total = total + sum(distance(ix, ii));
  end
  iter = iter+1;
end
end


function [index] = Roulettemethod(distance_matrix)

% Find shortest distance between one sample and its closest cluster centroid
[min_distance,~] = min(distance_matrix,[],2);

% Normalize for further operations
min_distance = min_distance ./ sum(min_distance);

% Construct roulette according to min_distance
temp_roulette = cumsum(sqrt(dot(min_distance',min_distance',1)));
temp_roulette = temp_roulette';
% Generate a random number for selection
% temp_rand = rand();
num = size(temp_roulette,1);
temp_rand = mean(temp_roulette);
% Find the corresponding index
for i = 1:num
    if((i == 1) && temp_roulette(i,1) > temp_rand)
        index = 1;
    elseif((temp_roulette(i,1) >= temp_rand) && (temp_roulette(i-1,1) < temp_rand))
        index = i;
    end
end
end
