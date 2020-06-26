% [idx,netsim,dpsim,expref]=apcluster(s,p)

%   Parameter    Value
%   'maxits'     Any positive integer. This specifies the
%                maximum number of iterations performed by
%                affinity propagation. Default: 500.
%   'convits'    Any positive integer. APCLUSTER decides that
%                the algorithm has converged if the estimated
%                cluster centers stay fixed for convits
%                iterations. Increase this value to apply a
%                more stringent convergence test. Default: 50.
%   'lam'        A real number that is less than 1 and
%                greater than or equal to 0.5. This sets the
%                damping level of the message-passing method,
%                where values close to 1 correspond to heavy
%                damping which may be needed if oscillations
%                occur.
%   'nonoise'    No value needed. Degenerate input similarities
%                (eg, where the similarity of i to k equals the
%                similarity of k to i) can prevent convergence.
%                To avoid this, APCLUSTER adds a small amount
%                of noise to the input similarities. This flag
%                turns off the addition of noise.
%
function [idx,unconverged]=apcluster(s,p)

% Handle arguments to function
if nargin<2 error('Too few input arguments');
else
    maxits=80; convits=50; lam=0.8; nonoise=0;
end

% Construct similarity matrix
N=length(p);
S=-Inf*ones(N,N); 
for j=1:size(s,1) 
    S(s(j,1),s(j,2))=s(j,3); 
end

% In case user did not remove degeneracies from the input similarities,
% avoid degenerate solutions by adding a small amount of noise to the
% input similarities
if ~nonoise
    rns=randn('state'); randn('state',0);
    S=S+(eps*S+realmin*100).*rand(N,N);
    randn('state',rns);
end

% Place preferences on the diagonal of S
for i=1:N S(i,i)=p(i); end

% Allocate space for messages, etc
A=zeros(N,N); R=zeros(N,N);

% Execute parallel affinity propagation updates
e=zeros(N,convits); dn=0; i=0;
while ~dn
    i=i+1; 
    % Compute responsibilities
    Rold=R;
    AS=A+S; [Y,I]=max(AS,[],2); 
    for k=1:N AS(k,I(k))=-realmax; end
    [Y2,I2]=max(AS,[],2);
    R=S-repmat(Y,[1,N]);
    for k=1:N R(k,I(k))=S(k,I(k))-Y2(k); end
    R=(1-lam)*R+lam*Rold; % Damping

    % Compute availabilities
    Aold=A;
    Rp=max(R,0);
    for k=1:N Rp(k,k)=R(k,k); end
    A=repmat(sum(Rp,1),[N,1])-Rp;
    dA=diag(A); A=min(A,0); for k=1:N A(k,k)=dA(k); end
    A=(1-lam)*A+lam*Aold; % Damping

    % Check for convergence
    E=((diag(A)+diag(R))>0); e(:,mod(i-1,convits)+1)=E; K=sum(E);
    if i>=convits || i>=maxits
        se=sum(e,2);
        unconverged=(sum((se==convits)+(se==0))~=N);
        if (~unconverged&&(K>0))||(i==maxits) dn=1; end
    end
end
idx=find(diag(A+R)>0); 

