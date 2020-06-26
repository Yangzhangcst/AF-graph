function [X, out]= solveF(X,fun,opts,varargin)
if isempty(X)
    error('input X is an empty matrix');
else
    [n, k] = size(X);
end

if isfield(opts, 'xtol')
    if opts.xtol < 0 || opts.xtol > 1
        opts.xtol = 1e-6;
    end
else
    opts.xtol = 1e-6;
end

if isfield(opts, 'gtol')
    if opts.gtol < 0 || opts.gtol > 1
        opts.gtol = 1e-6;
    end
else
    opts.gtol = 1e-6;
end

if isfield(opts, 'ftol')
    if opts.ftol < 0 || opts.ftol > 1
        opts.ftol = 1e-12;
    end
else
    opts.ftol = 1e-12;
end

% parameters for control the linear approximation in line search
if isfield(opts, 'rho')
   if opts.rho < 0 || opts.rho > 1
        opts.rho = 1e-4;
   end
else
    opts.rho = 1e-4;
end

% factor for decreasing the step size in the backtracking line search
if isfield(opts, 'eta')
   if opts.eta < 0 || opts.eta > 1
        opts.eta = 0.1;
   end
else
    opts.eta = 0.2;
end

% parameters for updating C by HongChao, Zhang
if isfield(opts, 'gamma')
   if opts.gamma < 0 || opts.gamma > 1
        opts.gamma = 0.85;
   end
else
    opts.gamma = 0.85;
end

if isfield(opts, 'tau')
   if opts.tau < 0 || opts.tau > 1e3
        opts.tau = 1e-3;
   end
else
    opts.tau = 1e-3;
end

% parameters for the  nonmontone line search by Raydan
if ~isfield(opts, 'STPEPS')
    opts.STPEPS = 1e-10;
end

if isfield(opts, 'nt')
    if opts.nt < 0 || opts.nt > 100
        opts.nt = 5;
    end
else
    opts.nt = 5;
end

if isfield(opts, 'projG')
    switch opts.projG
        case {1,2}; otherwise; opts.projG = 1;
    end
else
    opts.projG = 1;
end

if isfield(opts, 'iscomplex')
    switch opts.iscomplex
        case {0, 1}; otherwise; opts.iscomplex = 0;
    end
else
    opts.iscomplex = 0;
end

if isfield(opts, 'mxitr')
    if opts.mxitr < 0 || opts.mxitr > 2^20
        opts.mxitr = 1000;
    end
else
    opts.mxitr = 1000;
end

if ~isfield(opts, 'record')
    opts.record = 0;
end


%-------------------------------------------------------------------------------
% copy parameters
xtol = opts.xtol;
gtol = opts.gtol;
ftol = opts.ftol;
rho  = opts.rho;
STPEPS = opts.STPEPS;
eta   = opts.eta;
gamma = opts.gamma;
iscomplex = opts.iscomplex;
record = opts.record;

nt = opts.nt;   crit = ones(nt, 3);

invH = true; if k < n/2; invH = false;  eye2k = eye(2*k); end

%% Initial function value and gradient
% prepare for iterations
 [F,  G] = feval(fun,X, varargin{:}); 
%[F,G]=fun(varargin{:});
out.nfe = 1;  
GX = G'*X;

if invH
    GXT = G*X';  H = 0.5*(GXT - GXT');  RX = H*X;
else
    if opts.projG == 1
        U =  [G, X];    V = [X, -G];       VU = V'*U;
    elseif opts.projG == 2
        GB = G - 0.5*X*(X'*G);
        U =  [GB, X];    V = [X, -GB];       VU = V'*U;
    end
    %U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];   
    %VX = VU(:,k+1:end); %VX = V'*X;
    VX = V'*X;
end
dtX = G - X*GX;     nrmG  = norm(dtX, 'fro');
    
Q = 1; Cval = F;  tau = opts.tau;

%% Print iteration header if debug == 1
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '----------- Gradient Method with Line search ----------- \n');
    fprintf(fid, '%4s %8s %8s %10s %10s\n', 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff');
    %fprintf(fid, '%4d \t %3.2e \t %3.2e \t %5d \t %5d	\t %6d	\n', 0, 0, F, 0, 0, 0);
end

%% main iteration
for itr = 1 : opts.mxitr
    XP = X;     FP = F;   GP = G;   dtXP = dtX;
     % scale step size

    nls = 1; deriv = rho*nrmG^2; %deriv
    while 1
        % calculate G, F,
        if invH
            [X, infX] = linsolve(eye(n) + tau*H, XP - tau*RX);
        else
            [aa, infR] = linsolve(eye2k + (0.5*tau)*VU, VX);
            X = XP - U*(tau*aa);
        end
        %if norm(X'*X - eye(k),'fro') > 1e-6; error('X^T*X~=I'); end
        if ~isreal(X) && ~iscomplex ; error('X is complex'); end
   [F,  G] = feval(fun,X, varargin{:}); 
%[F,G]=fun(varargin{:});
        
        
        out.nfe = out.nfe + 1;
        
        if F <= Cval - tau*deriv || nls >= 5
            break;
        end
        tau = eta*tau;          nls = nls+1;
    end  
    
    GX = G'*X;
    if invH
        GXT = G*X';  H = 0.5*(GXT - GXT');  RX = H*X;
    else
        if opts.projG == 1
            U =  [G, X];    V = [X, -G];       VU = V'*U;
        elseif opts.projG == 2
            GB = G - 0.5*X*(X'*G);
            U =  [GB, X];    V = [X, -GB];     VU = V'*U; 
        end
        %U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];
        %VX = VU(:,k+1:end); % VX = V'*X;
        VX = V'*X;
    end
    dtX = G - X*GX;    nrmG  = norm(dtX, 'fro');
    
    S = X - XP;         XDiff = norm(S,'fro')/sqrt(n);
    tau = opts.tau; FDiff = abs(FP-F)/(abs(FP)+1);
    
    if iscomplex
        %Y = dtX - dtXP;     SY = (sum(sum(real(conj(S).*Y))));
        Y = dtX - dtXP;     SY = abs(sum(sum(conj(S).*Y)));
        if mod(itr,2)==0; tau = sum(sum(conj(S).*S))/SY; 
        else tau = SY/sum(sum(conj(Y).*Y)); end    
    else
        %Y = G - GP;     SY = abs(sum(sum(S.*Y)));
        Y = dtX - dtXP;     SY = abs(sum(sum(S.*Y)));
        %alpha = sum(sum(S.*S))/SY;
        %alpha = SY/sum(sum(Y.*Y));
        %alpha = max([sum(sum(S.*S))/SY, SY/sum(sum(Y.*Y))]);
        if mod(itr,2)==0; tau = sum(sum(S.*S))/SY;
        else tau  = SY/sum(sum(Y.*Y)); end
        
    end
    tau = max(min(tau, 1e20), 1e-20);
    
    if (record >= 1)
        fprintf('%4d  %3.2e  %4.3e  %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr, tau, F, nrmG, XDiff, FDiff, nls);
        %fprintf('%4d  %3.2e  %4.3e  %3.2e  %3.2e (%3.2e, %3.2e)\n', ...
        %    itr, tau, F, nrmG, XDiff, alpha1, alpha2);
    end
    
    crit(itr,:) = [nrmG, XDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
  
    if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol || all(mcrit(2:3) < 10*[xtol, ftol])  
        if itr <= 2
            ftol = 0.1*ftol;
            xtol = 0.1*xtol;
            gtol = 0.1*gtol;
        else
            out.msg = 'converge';
            break;
        end
    end
    
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + F)/Q;
 end

if itr >= opts.mxitr
    out.msg = 'exceed max iteration';
end

out.feasi = norm(X'*X-eye(k),'fro');
if  out.feasi > 1e-13
    X = MGramSchmidt(X);
    [F,  G] = feval(fun,X, varargin{:}); 
%[F,G]=fun(varargin{:});
    out.nfe = out.nfe + 1;
    out.feasi = norm(X'*X-eye(k),'fro');
end

out.nrmG = nrmG;
out.fval = F;
out.itr = itr;




end
% function [G,F]=fun1(P,alpha,Y,Q,L)
% G=2*L*P-2*alpha*Y*Q';
% F=trace(P'*L*P)+alpha*norm(Y-P*Q);
% end
% 
% function [G,F]=fun2(P,alpha,Y,Q,L)
% G=-2*P'*Y;
% F=norm(Y-P*Q);
% end