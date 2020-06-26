function Xp = dimReduction_PCA(X, dim)
%DIMREDUCTION_PCA Dimension reduction by PCA.
% 	XP = dimReduction_PCA(X,DIM) computes DIM-dimensional embedding by PCA.
%   Let X = U D V' be its SVD where the singular values in D are ordered. 
%   XP is computed as the first DIM rows of D * V'.

% Input Arguments
% X                 -- data matrix of size D by N.
% dim               -- dimension
% Output Arguments
% Xp                -- data matrix of size dim by N.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

LARGE_SCALE = 5000; % if min(D, N) > LARGE_SCALE then the data is too large.

[m, N] = size(X);

if dim == 0 || dim >= size(X',2) || dim >= N
    Xp = X;
    return;
end

if m <= N
    if m > LARGE_SCALE, error(['Error in ''' mfilename ''': matrix too big\n']); end;
    ddata = X * X';
    [U, eigval] = eig(ddata);
    [~, order] = sort(-diag(eigval));
    Xp = U(:, order(1:dim))' * X;
else
    if N > LARGE_SCALE, error(['Error in ''' mfilename ''': matrix too big\n']); end;
    ddata = X' * X;
    [V, eigval] = eig(ddata);
    [~, order] = sort(-diag(eigval));
    Xp = (eigval(order(1:dim), order(1:dim)) .^.5) * V(:, order(1:dim))';
end

end
            
