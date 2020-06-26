function [GCE_score] = GCE( assig1, assig2 )

k1 = max( assig1 );         % number clusters
k2 = max( assig2 );
n = length( assig1 );          % nuber data points

confusion = zeros( k1, k2 );     % build confusion matrix

for i1 = 1:k1;
for i2 = 1:k2;
  confusion( i1, i2 ) = length( find( (assig1 == i1) & (assig2 == i2 )));
end;
end;

GCE_score1 = 0; GCE_score2 = 0;
for i1 = 1:k1;
for i2 = 1:k2;
    idx1 = find(assig1== i1);
    idx2 = find(assig2== i2);
    
    diff_idx1 = setdiff(idx1,idx2);
    diff_idx2 = setdiff(idx2,idx1);
    
    GCE_score1 = GCE_score1 + size(diff_idx1,1)*confusion( i1, i2 )/size(idx1,1);
    GCE_score2 = GCE_score2 + size(diff_idx2,1)*confusion( i1, i2 )/size(idx2,1);
    clear idx1 idx2
end;
end;
GCE_score = min(GCE_score1,GCE_score2)/n;
