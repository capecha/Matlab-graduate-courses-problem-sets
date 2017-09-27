function [ mse ] = fitt(y,alpha, m, d )
%This function uses equation 5.11 to obtain predicted values for log(1-F)
%and with these compute the mean squeare error with respect observational
%data

yhat=(2/(1-alpha))*log(m+(2*alpha*m/(1-alpha)))...
    -(2/(1-alpha))*log(d+(2*alpha*m/(1-alpha)));
error=(y-yhat).^2;
mse=-mean(error);
end

