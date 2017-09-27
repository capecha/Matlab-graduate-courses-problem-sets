% Econometrics I
% IDEA 2012/2013
% Computer Problem Set 2 - Maximum Likelihood Estimation
%
%--------------------------------------------------------------------------
% This function calculates the hessian numerically.
%--------------------------------------------------------------------------
%
function H=hessian_programmed(beta,y,x)
%
h=0.00001;
[k,a]=size(beta');	
k = length(beta)
% k is the number of regressors
%
H=zeros(k,k); % matrix to be filled up
I=eye(k,k)*h;
% We will use this matrix to obtain the numerical derivatives
% It has the increments of each of the variables by columns, so we can
% add it to beta and see how the likelihood function changes.
%
% The procedure is the following:
%
% 1) Evaluate the gradient at beta0, using the m-file "gradientCPS2";
% 2) Evaluate the gradient at beta0+h_j, an infinitesimal increment of the 
%first-order derivative, using again the m-file "gradientCPS2";
% 3) Compute the numerical derivative of beta_j: (step2-step1)/h

for i=1:1:k	% for all rows of the Hessian
	for j=1:1:k	% for all columns
	g=gradientCPS2(beta,y,x);
	% For the row i of the hessian we take i-th element of the gradient:
	step1=g(i);	
    %
	% Now, introduce a small change in the parameter j by
	% adding column j of matrix M (zeros except for j-th element):	
	beta_h=beta'+I(:,j);
    %
	% Recalculate the gradient again and pick i-th element of the gradient
	% matrix again:
	g=gradientCPS2(beta_h',y,x);
	step2=g(i);
	% Last step, build the Hessian, one cell at a time:
	H(i,j)=(step2-step1)/h;
	end
end
%%
