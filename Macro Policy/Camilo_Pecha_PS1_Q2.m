% Macroeconomic Policy
% Camilo Pecha
% Problem set 1, question 2.1
% value function iteration for
% simple non-stochastic growth model

clear all
clc
% Functional Forms

% preferences: u(c)=log(ci+cbar)
% beta: discount factor
% delta is the rate of capital depreciation
% the production function is: k^alpha

% Basic Parameters

     beta=.99;
     alpha=.36;
     delta=0.025;
     theta=0.5;
     cbar=0.0001;
% Spaces
% here is specifyed the space that the capital stock lives in
% also is made an educated guess here so that the capital stock space
% includes the steady state and the hypothetical value ki=theta*Kss

% N is the size of the state space k and k' (i.e. the number of points in
% the grid)
N=1000;

% determine the steady state capital stock that we found analiticaly

Kss=(((1/(alpha*beta))-((1-delta)/alpha)))^(1/(alpha-1));

% put in bounds here for the capital stock space
% bounds are adjusted to let ki live in the state space 

Kl=Kss*.49; %the lower bound of K
Ku=Kss*1.1; %the upper bound of K
step=(Ku-Kl)/N;

K=Kl:step:Ku; % state space as a row vector
% n is the true length of the state space
n=length(K);
%steady state prices
wss=(1-alpha)*Kss^alpha;
rss=alpha*Kss^(alpha-1);
%initial capital for i
Kssi=theta*Kss;
%consumption of i given k0=theta*Kss
c=wss+(rss-delta)*(Kssi);
%Value function evaluated at the ki with steady state prices
v_ki=(1/(1-beta))*(log(c+cbar));

%Vector of production
Kalpha=K.^alpha;%the production function, a vector of capital values
% ytot is then total output available
ytot=Kalpha+(1-delta)*K; %output has size 1x1000

% VALUE FUNCTION ITERATION
% first we obtain an initial guess of the value
% function from the single period problem. 



% Matrix I plays the role of the future capital stock
colones=ones(n,1);
rowones= colones';
% I here is matrix where each column is K
I = K'*rowones;
% J is a matrix where each row is K
J= colones*K;
%consumption at each point in this space
Cbar=cbar*ones(n,n);
C = (J.^alpha)+(1-delta)*J -I ;

% current utility, current state as columns, future capital as rows
U=log(max((C+Cbar),Cbar)); %that is this is the initial guess on V


% Second we start the value function iteration routine.
V=zeros(1,n); %n is each value of k

% Iterations: using this first guess at V iterate 

T=100; % here T is the maximal number of iterations

%here we set our tolerance for convergence; 

toler=.0001; %the difference between value functions



for j=1:T;
w=ones(n,1)*V;
q=w';
r = U + beta*q;
v=max(r);
diff=sum((abs(V-v)));
if abs(diff) <= toler
   break
else
     V=v;
 end
end

%CHECK how v looks like
figure(1)
plot(K,V)
xlabel('Capital')
ylabel('value function')

% Having the vector v, I want to know the policy vector using the max to figure out which row was chosen for each
% level of the capital stock

[R,m]=max(r);

% now build a vector for future capital using the m-vector.

Kprime=[];

% now loop through the m-vector using the I matrix

for i=1:n
inv=I(m(i),i);
Kprime=[Kprime inv];
end

% here we plot the policy function against the 45 degree line

figure(2)
plot(K,Kprime,'r');
hold on;
plot(K,K);% so now we have a nice 45 degree line
xlabel('current capital')
ylabel('next period capital')
legend('policy function','current capital',0)

% use this plot to look for steady state and investment patterns
figure(3)
plot(K,Kprime-K)
xlabel('current capital')
ylabel('net investment')
% this command will plot the index of the policy function
%plot(m); % 

% simulation of transition dynamics

   P=130; % arbitrary length for transition path
   capind=ones(1,P);
   captran=ones(1,P);
   captran(1)= Kssi(capind(1)); %assuming that the economy starts in the ki=thetaKss
   
   for t=2:P
   capind(t)=(m(capind(t-1)));% so follow evolution in index space   
   captran(t)=K(capind(t)); % follow evoluation in capital space
	end
figure(4)    
plot(captran)
xlabel('time period')
ylabel('capital stock')
  
  
%Variables using the transition values 
ytr=captran.^alpha;
yss=Kss.^alpha;
Kptran=lagmatrix(captran,-1)';
%consumption with transition values

Ctran=ytr+(1-delta).*captran-Kptran;
Css=yss-delta.*Kss;
Cbartr=cbar*ones(1,P);
%Utility with trnasition values of k
Utran=log(Ctran+Cbartr);
Uss=log(Css+Cbartr);

figure(5)
plot(Ctran)
xlabel('time')
ylabel('consumption')


%Lifetime welfare


%Utility at the steady state is defined by v_ki found before

Utran=log(Ctran+Cbartr); %this is the utility given by the transition values

cond=captran-Kptran
cond(P)=0
stop=find(cond,1,'last');

V1=[];
V1(1)=Utran(1);
for i=2:stop
    V1(i)=(beta^i)*Utran(i);
end

Vnew=sum(V1)+(beta^stop/(1-beta))*Uss;

Udiff=v_ki-Vnew


%As can be seen, moving the value of theta and open the state space the
%decision of being in the economy or leave to an island change.

if Udiff>toler;
    disp('Do not go to the island since the welfare form leaving the economy is less than that if staying in.')
end





  
