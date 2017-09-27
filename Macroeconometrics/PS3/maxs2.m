function [ xop ] = maxs2( D,x0 )
%This is the function that being maximized delivers the maximal shock
%effect.
x1=x0(1,1);
x2=x0(2,1);
x3=x0(3,1);
x4=x0(4,1);
x5=x0(5,1);
x6=x0(6,1);

xop=-(D(2,2)*x1+D(2,3)*x2+D(2,4)*x3+D(2,5)*x4+D(2,6)*x5+ ...
    D(2,7)*x6+D(2,8)*sqrt(1-(x1^2+x2^2+x3^2+x4^2+x5^2+x6^2)));


end

