% Econometrics I
% IDEA 2012/2013
% Computer Problem Set 2 - Maximum Likelihood Estimation
%
%--------------------------------------------------------------------------
% This version uses the command FMINSEARCH.
% The name "simple" just means that sigma is assumed to be 1, as given in
% the exercise.
% FMINSEARCH finds the minimum of a scalar function of several variables, 
% starting at an initial estimate. This is generally referred to as 
% unconstrained nonlinear optimization.
%--------------------------------------------------------------------------
%
% Load data:
%
clear all;
clc;
load ComputerPS2_MLE_dataset.mat;
%
% 4 - Estimate parameters by ML: 
% 
% Create experience squared:
%
exper2=exper.^2;
%
% Declare variables of the model:
%
N=size(logwk,1);
x=cat(2,ones(N,1),educ,exper,exper2);
y=logwk;
%
% Estimate parameters using OLS:
%
b_OLS=inv(x'*x)*x'*y;
disp(b_OLS)
%
% When moving to ML estimation, take the OLS estimators as starting values:
%
beta0=b_OLS';
[beta,fval,exitflag,output] = fminsearch(@(beta) LogLik(N,y,x,beta),beta0);
BaseLLF=-fval; % we will use this in question 5, for the LR (denominator)
%
%
%5 - Hypothesis Testing
%5.1 - Likelihood Ratio Test
%
% We first need to specify the Likelihood when the null-hypothesis is true:
%
beta0=b_OLS';
A=[];
b=[];
Aeq=eye(4,4);
beq=[beta0(1) 0 beta0(3) beta0(4)];
[beta_null,fval,exitflag,output] = fmincon(@(beta) LogLik(N,y,x,beta),beta0,A,b,Aeq,beq);
NullLLF=-fval; % numerator of the Likelihood Ratio
%
% Using the two values of the Likelihoods:
%
LR = -2*(BaseLLF - NullLLF);
p_value=chi2cdf(LR,1); % because degrees of freedom are the difference between the "free"
% parameters under the alternative and under the null, i.e. 4-3=1.
%
% 5.2 - Wald Test
%
% For the Wald test, we need to compute matrices h_prime and v_hat, where:
% (i)  h_prime is the gradient of the null hypothesis (wrt the parameters);
% (ii) v_hat is the estimated asymptotic variance.
%
% How to estimate these?
%
% (i)  First, let us compute the gradient. For that we need to formalise
% the null hypotesis, R*beta=r.
%
R=[0 1 0 0];
r=0;
%
% In this simple null hypotesis, we know that:
%
h_prime=R;
%
% (ii) Second, we must compute the hessian matrix:
%
H=hessian_programmed(beta,y,x);
%
% Now, we have all the information to compute the Wald test:
%
W=(R*beta_null'-r)'*inv(R*H*R')*(R*beta_null'-r);
alpha=0.05;
compare=chi2inv(1-alpha,1);
if compare > W
    fprintf('We reject H0 that the coefficient on education is not statistically significant.')
end
%
% We can also compute the p-values:
%
p_value2=chi2cdf(W,1);
%
%%







