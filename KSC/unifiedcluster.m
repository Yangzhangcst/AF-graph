function [Z] = unifiedcluster(K,alpha,beta)

K1 = K*K';
K2 = K1 ./ max(K1(:));% linear kernel

addpath('KSC/qpc')
[~,n]=size(K2);
Z=eye(n);
for i=1:200
    Zold=Z;
    Z = (Z+Z')/2;
    D = diag(sum(Z));
    L = D-Z; 
    F=eig1(L, n, 0);

    parfor ij=1:n
        all = veccomp2(ij,n,F);
        all=all.*all;
        H=2*alpha*eye(n)+2*K2;
        H=(H+H')/2;
        ff=beta/2*all'-2*K2(:,ij);
        Z(:,ij) = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
    end
    err = norm(Z-Zold)/norm(Zold);
    if err<1e-3,  break;  end

end
end

function [all]=veccomp2(ij,n,F)
    for ji=1:n
        all(ji)=norm(F(ij,:)-F(ji,:));
    end
end
