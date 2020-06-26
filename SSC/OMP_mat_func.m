function C = OMP_mat_func(X, K, thr)
%OMP_MAT_FUNC Perform OMP for self representation.
%   This code implements the subspace clustering algorithm described in
% 
%   Chong You, Daniel Robinson, Rene Vidal,
%   "Scalable Sparse Subspace Clustering by Orthogonal Matching Pursuit", 
%   CVPR 2016.
% 
% 	It perform OMP for each column of data X = [x_1, \dots, x_N]
%   using all other columns as a dicitonary. 
%   i.e., for each j = 1, \dots, N, compute the following by OMP:
%   \min_{c_j} \| x_j - X c_j \|_F^2 s.t. \|c_j\|_0 \le K, c_jj = 0.
%   The output C is given by [c_1, \dots, c_N].

% Input Arguments
% X                 -- data matrix D by N where each column is a data point.
% K                 -- termination by checking the number of nonzero 
%                      entries in c_j, i.e. OMP terminates if \|c_j\| >= K
% thr               -- termination by checking the reconstruction error, 
%                      i.e. OMP terminates if \| x_j - X c_j \|_2^2 < thr

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com


MEMORY_TOTAL = 0.1 * 10^9; % memory available for double precision.
[~, N] = size(X);

Xn = X; % assume that data is column normalized. If not, uncomment the following.
Xn = cnormalize(X);

S = ones(N, K); % Support set
Ind = repmat((1:N)', 1, K);
Val = zeros(N, K); % Nonzero Value 
t_vec = ones(N, 1) * K;

res = Xn; % residual
for t = 1:K
    blockSize = round(MEMORY_TOTAL / N);
    counter = 0;
    while(1)
        mask = [counter+1 : min(counter + blockSize, N)];
        I = abs(X' * res(:, mask));
        I(counter+1:N+1:end) = 0; % set diagonal = 0
        [~, J] = max(I, [], 1);
        S(mask, t) = J;
        counter = counter + blockSize;
        if counter >= N
            break;
        end
    end
    
    if t ~= K % not the last step. compute residual
        for iN = 1:N	 
            if t_vec(iN) == K % termination has not been reached
                B = Xn(:, S(iN, 1:t));
                res(:, iN) = Xn(:, iN) - B* (B \ Xn(:, iN));
                if sum( res(:, iN).^2 ) < thr
                    t_vec(iN) = t;
                end 
            end
        end
    end
    if sum(t_vec == K) == 0
        break;
    end
%     fprintf('Step %d in %d\n', t, K);
end

% compute coefficients
for iN = 1:N	 
    Val(iN, 1:t_vec(iN)) = (X(:, S(iN, 1:t_vec(iN))) \ X(:, iN))'; % use X rather than Xn
end

C = sparse( S(:), Ind(:), Val(:), N, N );

