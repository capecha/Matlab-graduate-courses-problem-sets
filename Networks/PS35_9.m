clc
clear all
%setting the variables
d=[0
    1
    2
    3
    4
    5
    6
    7
    8]; %degrees
p=[7
    17
    11
    9
    12
    3
    4
    3
    1]; %prisioners

m=mean(d)/2; %this is the value of m that is roughly half of the mean degree

pt=sum(p); %this is the total number of prosioners 
X=[d p];   
fr=zeros(size(d));
F=zeros(size(d));

for i=1:9
    fr(i,1)=p(i,1)./pt;
    F=fr;
end

for i=2:9
    F(i,1)=F(i-1,1)+fr(i,1); %The coummulative distribution function from data
end
F(9,1)=0.99; %setting the last value below 1, other wise no solution
y=zeros(size(d));

for i=1:9
    y(i,1)=log(1-F(i,1))    ; %this is  the left hand side of eqn 5.11 in Jackson's book
end
plot(F)
a=0.5;
fitt(y,a,m,d)
[a]= fminsearch(@(a) fitt(y,a,m,d), a ) %the program that seeks for the minimum of the mean square error wich result is a
