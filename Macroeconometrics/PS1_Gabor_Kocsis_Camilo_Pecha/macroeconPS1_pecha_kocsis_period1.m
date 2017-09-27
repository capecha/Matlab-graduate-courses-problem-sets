%% MACROECONOMETRICS PS#1
% Camilo José Pecha and Gábor Kocsis

% USING DATA FOR THE PERIOD 1970-1983

clc
clear all
load SWdata
obs=length(gdpdef);
pi=zeros(obs,1);

for i=54:obs
    pi(i,1)=400*log(gdpdef(i,1)/gdpdef(i-1,1));
end

%benchmark AR(2), that is an OLS 
%vars


%h=1
dpi=zeros(obs,1);
for i=2:obs
   dpi(i,1)= pi(i,1)-pi(i-1,1);
end

XLAG = lagmatrix(dpi,[0 1 2 3 4]);
eb1=zeros(200,1);

for i=93:148
y=[XLAG(53:i,1)];
cons=ones(length(y),1);
%for 1 lag
X1=[cons XLAG(53:i,2) XLAG(53:i,3)];
%ols 

%two lags
phi1=inv(X1'*X1)*X1'*y;
yhat=X1*phi1;
yobs=y(i-53,1);
yfor=yhat(i-53,1);
eb1(i,1)=(yobs-yfor)^2;

end

epsAR2h1=eb1(93:148);
mseAR2h1=mean(epsAR2h1);

%h=2
dpi2=zeros(obs,1);
for i=3:obs
   dpi2(i,1)= (1/2)*(pi(i,1)+pi(i-1,1))-pi(i-2,1);
end

XLAG2 = lagmatrix(dpi2,[0 1 2 3 4]);
eb2=zeros(200,1);

for i=93:148
y=[XLAG2(53:i,1)];
cons=ones(length(y),1);
%for 2 lags
X2=[cons XLAG2(53:i,2) XLAG2(53:i,3)];
%ols 

%two lags
phi2=inv(X2'*X2)*X2'*y;
yhat=X2*phi2;
yobs=y(i-53,1);
yfor=yhat(i-53,1);
eb2(i,1)=(yobs-yfor)^2;
end
epsAR2h2=eb2(93:148);
mseAR2h2=mean(epsAR2h2);

%h=4
dpi4=zeros(obs,1);
for i=5:obs
   dpi4(i,1)= (1/4)*(pi(i,1)+pi(i-1,1)+pi(i-2,1)+pi(i-3,1))-pi(i-4,1);
end

XLAG4 = lagmatrix(dpi4,[0 1 2 3 4]);
eb1=zeros(200,1);
eb2=zeros(200,1);
eb3=zeros(200,1);
eb4=zeros(200,1);

for i=93:148
y=[XLAG4(53:i,1)];
cons=ones(length(y),1);
%for 1 lag
X1=[cons XLAG4(53:i,2)];
%for 2 lags
X2=[cons XLAG4(53:i,2) XLAG4(53:i,3)];
%for 3 lags
X3=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAG4(53:i,4)];
%for 4 lags
X4=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAG4(53:i,4) XLAG4(53:i,5)];
%ols 

%two lags
phi1h4=inv(X1'*X1)*X1'*y;
yhat1h4=X1*phi1h4;
yobs=y(i-53,1);
yfor1h4=yhat1h4(i-53,1);
eb1h4(i,1)=(yobs-yfor1h4)^2;

%two lags
phi2h4=inv(X2'*X2)*X2'*y;
yhat2h4=X2*phi2h4;
yobs=y(i-53,1);
yfor2h4=yhat2h4(i-53,1);
eb2h4(i,1)=(yobs-yfor)^2;

%four lags
phi4h4=inv(X4'*X4)*X4'*y;
yhat4h4=X4*phi4h4;
yobs=y(i-53,1);
yfor4h4=yhat4h4(i-53,1);
eb4h4(i,1)=(yobs-yfor4h4)^2;
end

epsAR1h4=eb1h4(93:148);
mseAR1h4=mean(epsAR1h4);

epsAR2h4=eb2h4(93:148);
mseAR2h4=mean(epsAR2h4);

epsAR4h4=eb4h4(93:148);
mseAR4h4=mean(epsAR4h4);

%we should display as benchmark only mseAR2h1 mseAR2h2 mseAR2h4

%Computation of AO forcast
%h=1
AO1=zeros(obs,1);
for i=2:obs
    AO1(i,1)=pi(i,1);
end

%h=2
AO2=zeros(obs,1);
for i=3:obs
    AO2(i,1)=(1/2)*(pi(i,1)+pi(i-1,1));
end

%h=4
AO4=zeros(obs,1);
for i=5:obs
    AO4(i,1)=(1/4)*(pi(i,1)+pi(i-1,1)+pi(i-2,1)+pi(i-3,1));
end

%compute errors

epsilonh1=zeros(length(cons),1);
for i=93:148
epsilonh1(i,1)=pi(i,1)-AO1(i-1,1);
end
epsilonh1sq=epsilonh1.^2;
mseh1=mean(epsilonh1sq);

AOh1=mseh1/mseAR2h1

epsilonh2=zeros(length(cons),1);
for i=93:148
epsilonh2(i,1)=pi(i,1)-AO2(i-1,1);
end
epsilonh2sq=epsilonh2.^2;
mseh2=mean(epsilonh2sq);

AOh2=mseh2/mseAR2h2

epsilonh4=zeros(length(cons),1);
for i=93:148
epsilonh4(i,1)=pi(i,1)-AO4(i-1,1);
end
epsilonh4sq=epsilonh4.^2;
mseh4=mean(epsilonh4sq);

AOh4=mseh4/mseAR2h4

%Phillips curve
%first, un employment lags to calcylate PC-u

u=UNRATE;

%h=1
du=zeros(obs,1);
for i=6:obs
   du(i,1)= u(i,1)-u(i-1,1);
end

XLAGu = lagmatrix(du,[0 1 2 3 4]);
eb1=zeros(200,1);

for i=93:148
y=[XLAG(53:i,1)];
cons=ones(length(y),1);
%for 1 lag
X1=[cons XLAG(53:i,2) XLAG(53:i,3) u(53:i,1) XLAGu(53:i,2) XLAGu(53:i,3)];
%ols 

%two lags
beta1=inv(X1'*X1)*X1'*y;
yhat=X1*beta1;
yobs=y(i-53,1);
yfor=yhat(i-53,1);
eb1(i,1)=(yobs-yfor)^2;

end

epsPC_u1=eb1(93:148);
msePC_u1=mean(epsPC_u1);

%h=2
du2=zeros(obs,1);
for i=7:obs
   du2(i,1)= (1/2)*(u(i,1)+u(i-1,1))-u(i-2,1);
end

XLAGu2 = lagmatrix(du2,[0 1 2 3 4]);
eb2=zeros(200,1);

for i=93:148
y=[XLAG2(53:i,1)];
cons=ones(length(y),1);
%for 2 lags
X2=[cons XLAG2(53:i,2) XLAG2(53:i,3) u(53:i,1) XLAGu2(53:i,2) XLAGu2(53:i,3)];
%ols 


%two lags
beta2=inv(X2'*X2)*X2'*y;
yhat=X2*beta2;
yobs=y(i-53,1);
yfor=yhat(i-53,1);
eb2(i,1)=(yobs-yfor)^2;
end
epsPC_u2=eb2(93:148);
msePC_u2=mean(epsPC_u2);


%h=4
du4=zeros(obs,1);
for i=9:obs
   du4(i,1)= (1/4)*(u(i,1)+u(i-1,1)+u(i-2,1)+u(i-3,1))-u(i-4,1);
end

XLAGu4 = lagmatrix(du4,[0 1 2 3 4]);
eb1=zeros(200,1);
eb2=zeros(200,1);
eb3=zeros(200,1);
eb4=zeros(200,1);

for i=93:148
y=[XLAG4(53:i,1)];
cons=ones(length(y),1);
%for 1 lag
X1=[cons XLAG4(53:i,2) u(53:i,1) XLAGu4(53:i,2)];
%for 2 lags
X2=[cons XLAG4(53:i,2) XLAG4(53:i,3) u(53:i,1) XLAGu4(53:i,2) XLAGu4(53:i,3)];
%for 3 lags
X3=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAG4(53:i,4) u(53:i,1) XLAGu4(53:i,2) XLAGu4(53:i,3) XLAGu4(53:i,4)];
%for 4 lags
X4=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAG4(53:i,4) XLAG4(53:i,5) u(53:i,1) XLAGu4(53:i,2) XLAGu4(53:i,3) XLAGu4(53:i,4) XLAGu4(53:i,5)];
%ols 

%two lags
beta1h4=inv(X1'*X1)*X1'*y;
yhat1h4=X1*beta1h4;
yobs=y(i-53,1);
yfor1h4=yhat1h4(i-53,1);
eb1h4(i,1)=(yobs-yfor1h4)^2;

%two lags
beta2h4=inv(X2'*X2)*X2'*y;
yhat2h4=X2*beta2h4;
yobs=y(i-53,1);
yfor2h4=yhat2h4(i-53,1);
eb2h4(i,1)=(yobs-yfor)^2;

%four lags
beta4h4=inv(X4'*X4)*X4'*y;
yhat4h4=X4*beta4h4;
yobs=y(i-53,1);
yfor4h4=yhat4h4(i-53,1);
eb4h4(i,1)=(yobs-yfor4h4)^2;
end

epsPC_u1h4=eb1h4(93:148);
msePC_u1h4=mean(epsPC_u1h4);

epsPC_u2h4=eb2h4(93:148);
msePC_u2h4=mean(epsPC_u2h4);

epsPC_u4h4=eb4h4(93:148);
msePC_u4h4=mean(epsPC_u4h4);

%rations with benchmark AR2

PC_u1=msePC_u1/mseAR2h1
PC_u2=msePC_u2/mseAR2h2
PC_u4=msePC_u2h4/mseAR2h4


%Second, un employment lags to calcylate PC-dy

y=log(GDP);

%h=1
dy=zeros(obs,1);
for i=2:obs
   dy(i,1)= y(i,1)-y(i-1,1);
end

XLAGy = lagmatrix(dy,[0 1 2 3 4]);
eb1=zeros(200,1);

for i=93:148
Y=[XLAG(53:i,1)];
cons=ones(length(Y),1);
%for 1 lag
X1=[cons XLAG(53:i,2) XLAG(53:i,3) XLAGy(53:i,2) XLAGy(53:i,3)];
%ols 

%two lags
beta1=inv(X1'*X1)*X1'*Y;
yhat=X1*beta1;
yobs=Y(i-53,1);
yfor=yhat(i-53,1);
eb1(i,1)=(yobs-yfor)^2;

end

epsPC_dy1=eb1(93:148);
msePC_dy1=mean(epsPC_dy1);

%h=2
dy2=zeros(obs,1);
for i=3:obs
   dy2(i,1)= (1/2)*(y(i,1)+y(i-1,1))-y(i-2,1);
end

XLAGy2 = lagmatrix(dy2,[0 1 2 3 4]);
eb2=zeros(200,1);

for i=93:148
Y=[XLAG2(53:i,1)];
cons=ones(length(Y),1);
%for 2 lags
X2=[cons XLAG2(53:i,2) XLAG2(53:i,3) XLAGy2(53:i,2) XLAGy2(53:i,3)];
%ols 


%two lags
beta2=inv(X2'*X2)*X2'*Y;
yhat=X2*beta2;
yobs=Y(i-53,1);
yfor=yhat(i-53,1);
eb2(i,1)=(yobs-yfor)^2;
end
epsPC_dy2=eb2(93:148);
msePC_dy2=mean(epsPC_dy2);


%h=4
dy4=zeros(obs,1);
for i=5:obs
   dy4(i,1)= (1/4)*(y(i,1)+y(i-1,1)+y(i-2,1)+y(i-3,1))-y(i-4,1);
end

XLAGy4 = lagmatrix(dy4,[0 1 2 3 4]);
eb1=zeros(200,1);
eb2=zeros(200,1);
eb3=zeros(200,1);
eb4=zeros(200,1);

for i=93:148
Y=[XLAG4(53:i,1)];
cons=ones(length(Y),1);
%for 1 lag
X1=[cons XLAG4(53:i,2) XLAGy4(53:i,2)];
%for 2 lags
X2=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAGy4(53:i,2) XLAGy4(53:i,3)];
%for 3 lags
X3=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAG4(53:i,4) XLAGy4(53:i,2) XLAGy4(53:i,3) XLAGy4(53:i,4)];
%for 4 lags
X4=[cons XLAG4(53:i,2) XLAG4(53:i,3) XLAG4(53:i,4) XLAG4(53:i,5) XLAGy4(53:i,2) XLAGy4(53:i,3) XLAGy4(53:i,4) XLAGy4(53:i,5)];
%ols 

%two lags
beta1h4=inv(X1'*X1)*X1'*Y;
yhat1h4=X1*beta1h4;
yobs=Y(i-53,1);
yfor1h4=yhat1h4(i-53,1);
eb1h4(i,1)=(yobs-yfor1h4)^2;

%two lags
beta2h4=inv(X2'*X2)*X2'*Y;
yhat2h4=X2*beta2h4;
yobs=Y(i-53,1);
yfor2h4=yhat2h4(i-53,1);
eb2h4(i,1)=(yobs-yfor)^2;

%foyr lags
beta4h4=inv(X4'*X4)*X4'*Y;
yhat4h4=X4*beta4h4;
yobs=Y(i-53,1);
yfor4h4=yhat4h4(i-53,1);
eb4h4(i,1)=(yobs-yfor4h4)^2;
end

epsPC_dy1h4=eb1h4(93:148);
msePC_dy1h4=mean(epsPC_dy1h4);

epsPC_dy2h4=eb2h4(93:148);
msePC_dy2h4=mean(epsPC_dy2h4);

epsPC_dy4h4=eb4h4(93:148);
msePC_dy4h4=mean(epsPC_dy4h4);

%rations with benchmark AR2

PC_dy1=msePC_dy1/mseAR2h1
PC_dy2=msePC_dy2/mseAR2h2
PC_dy4=msePC_dy2h4/mseAR2h4


