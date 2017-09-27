function LogLik = LogLik(N,y,x,beta);
LogLik = (N/2)*log(2*pi)+(1/2)*(sum((y-x*beta').^2));