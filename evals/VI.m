function [ distance, H1given2, H2given1, I12 ] = VI( assig1, assig2 )

%function distance = VI( assig1, assig2 )
%
% Computes the Variation of Information distance between two clusterings
% given by assig1,2. The values in assig1 must range between 1 and k1, 
% the number of clusters in assig1. Similarly for assig2.
%
%
%assig1, assig2(1,n) = vectors of integers from 1:k that represent the
%                      assignment to clusters. must have the same length
%

k1 = max( assig1 );         % number clusters
k2 = max( assig2 );
n = length( assig1 );          % nuber data points

confusion = zeros( k1, k2 );     % build confusion matrix

for i1 = 1:k1;
for i2 = 1:k2;
  confusion( i1, i2 ) = length( find( (assig1 == i1) & (assig2 == i2 )));
end;
end;

confusion = confusion/n;
p1 = sum( confusion, 1 );
p2 = sum( confusion, 2 );

idummy = find( confusion > 0 );
pdummy = p2*p1;
cdummy = confusion;
cdummy( idummy ) = confusion( idummy )./pdummy( idummy );
I12 = sum( sum( confusion.* locallog( cdummy )));
H1given2 = -sum( p1.*locallog( p1 ))-I12;
H2given1 =  -sum( p2.*locallog( p2 ))-I12;
distance = H1given2 + H2given1;

%distance = -sum( p1.*locallog( p1 ))-sum( p2.*locallog( p2 ))-2*sum( sum( confusion.* locallog( cdummy )));
distance = distance/log(2);
I12 = I12/log(2);
H1given2 = H1given2/log(2);
H2given1 = H2given1/log(2);



