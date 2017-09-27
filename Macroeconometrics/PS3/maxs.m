function [ xop ] = maxs( D,x0 )
%This is the function that being maximized delivers the maximal shock
%effect.
x1=x0(1,1);
x2=x0(2,1);

xop=-(D(2,2)*x1+D(2,3)*x2+D(2,4)*sqrt(1-(x1^2+x2^2)));


end

