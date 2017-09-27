% Econometrics I
% IDEA 2012/2013
% Computer Problem Set 2 - Maximum Likelihood Estimation
%
%--------------------------------------------------------------------------
% This version uses the command MAXLIK.
% MAXLIK takes the loglikelihood function and minimizes it in a loop form.
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
sigma_OLS = (1/(N-1))*((y-x*b_OLS)'*(y-x*b_OLS));
beta0=b_OLS;
sigma0=sigma_OLS;
disp(b_OLS)
disp(sigma_OLS)
%
% When moving to ML estimation, take the OLS estimators as starting values:
%
b = cat(1,beta0,sigma0);
result = maxlik(@(beta) LogLik(N,y,x,beta),b,info,y,x);
%
% Declare the results:
%
b_MLE = result.b;
v_MLE = inv(result.hess);
sterr_MLE = diag(v_MLE);
%%




