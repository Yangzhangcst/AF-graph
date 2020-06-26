function y = locallog( x )

% function y = locallog( x )
%
% returns log( x ) or zero if x == 0

y = x;
ii = find( y );
y( ii ) = log( y( ii ));
