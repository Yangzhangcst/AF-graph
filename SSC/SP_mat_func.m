function C = SP_mat_func(X, K, thr)

warning off

MEMORY_TOTAL = 0.1 * 10^9; % memory available for double precision.
[~, N] = size(X);

% assume that data is column normalized. If not, uncomment the following.
X = smoothKDE(X',10,10);
Xn = cnormalize(X);

S = ones(N, K); % Support set
B = ones(N, K);
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
        
        [~, J] = sort(I,'descend');
        
        counter = counter + blockSize;
        if counter >= N
            break;
        end
    end
    
    if t ~= K+1 % not the last step. compute residual
        for iN = 1:N	 
            if t_vec(iN) == K % termination has not been reached
                SJ = J(1:t,iN);
                A = Xn(:,SJ)\Xn(:, iN);
                [~, JA] = sort(abs(A),'descend');
                SJ = SJ(JA(1:t,:));
                A = A(JA(1:t,:));
                res(:, iN) = Xn(:, iN) - Xn(:, SJ)*A;
                if sum( res(:, iN).^2 ) < thr
                    t_vec(iN) = t;
                end
                S(iN,1:t) = SJ;
                B(iN,1:t) = A;
                clear SJ JA
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
%     Val(iN, 1:t_vec(iN)) = (X(:, S(iN, 1:t_vec(iN))) \ X(:, iN))'; % use X rather than Xn
    C = X(:, S(iN, 1:t_vec(iN))) \ X(:, iN);
    [~, JA] = sort(abs(C),'descend');
    C = C(JA(1:t_vec(iN),:));
    Val(iN, 1:t_vec(iN)) = C';
end
C = sparse( S(:), Ind(:), Val(:), N, N );
