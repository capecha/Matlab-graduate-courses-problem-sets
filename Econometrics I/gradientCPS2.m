% Econometrics I
% IDEA 2012/2013
% Computer Problem Set 2 - Maximum Likelihood Estimation
%
%--------------------------------------------------------------------------
% This function evaluates the gradient of the Likelihood Function.
%--------------------------------------------------------------------------
%
function g=gradientCPS2(beta,y,x)
% Set up function to be filled up (the gradient vector):
g=zeros(length(beta),1);
% Compute the gradient:
stop_val=length(beta);
for i=1:1:stop_val;
    g=x'*(y-x*beta');
end
%%

