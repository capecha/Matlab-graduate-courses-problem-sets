%Camilo Pecha
%Macroeconomic Policy
%Problem set 2, question 4


close all;
clc 
 
%Parameters: 
e_l = 0.99;                 % Endowment when unemployed
e_h = 1.1;                  % Endowment when employed
b   = .99;                  % Beta           
gs   = 200;                  % Grid Size         
pi   = [0.9 0.1;0.2 0.8];   % Transition probability matrix
tol = 0.0000001;             % Tolerance
one = ones(gs,1);                    %


%Assets and asset price
a_grid  = linspace(-2.5,5,gs);
qhi     = 3;
qlo     = b/4;
conv    = 100;
t       = 1;
 
 
while (conv > tol) & (t < 30)
   
    q   = (qhi + qlo) / 2;
    for i = 1:gs
        asset   = a_grid(1,i);
        for j = 1:gs
            asset_np  = a_grid(1,j); %assets in t+1
           
            c_h(i,j) = asset + e_h - q*asset_np; %consumption when employed
            c_l(i,j) = asset + e_l - q*asset_np; %consumption when unemployed 
           
            u_h(i,j) = log(c_h(i,j)); %utility when employed
            u_l(i,j) = log(c_l(i,j)); %utility when unemployed
           
            if c_h(i,j) <= 0
                u_h(i,j) = -10000000000;
            end
            if c_l(i,j) <= 0
                u_l(i,j) = -10000000000;
            end
           
        end
    end
   %Initialization of value function
    v1 = zeros(1,gs); 
    v2 = v1;
    pos     = zeros(gs,2);
    conv1   = 100;
    conv2   = 100;
   
    while (conv1 > .001) | (conv2 ~= 0)
       
        v_h = (u_h + b*one*( pi(1,1)*v1 + pi(1,2)*v2 ) ); %value function for employed
        v_l = (u_l + b*one*( pi(2,1)*v1 + pi(2,2)*v2 ) ); %value function for unemployed
       
        [V1,pos1]  = max( v_h' ); %detecting the maximum values for each value function
        [V2,pos2]  = max( v_l' );
       
        V           = [V1,V2];
        v           = [v1,v2];
        Pos         = [pos1',pos2']; 
       
        conv1       = abs(V - v);
        conv1       = max( (max(conv1))' );
        conv2       = max(any(Pos-pos));
       
        v1 = V1;    
        v2 = V2; 
        pos         = Pos;
       
    end
 
    posa = zeros(gs,gs); 
    posb = posa;
 
    for i = 1:gs
        posa(i,pos1(1,i)) = 1; 
        posb(i,pos2(1,i)) = 1;
    end
   
    t_matx  = [pi(1,1)*posa,pi(1,2)*posa;
        pi(2,1)*posb,pi(2,2)*posb];
       
    dist = (ones(2*gs,1))/(2*gs);
   
    conv1   = 100;
    while conv1 > tol
        dist_p     = t_matx'*dist;
        conv1   = max(abs(dist_p-dist));
        dist       = dist_p;
    end
   
    ExDD    = [a_grid,a_grid]*dist; 
    ExcessDemand(t) = ExDD;
    Q(t)    = q;
    if (t>1) & (abs(Q(t) - Q(t-1))<tol)
        break
    end
    if ExDD < 0
        qhi = (qhi + qlo) / 2;
    else
        qlo = (qhi + qlo) / 2;
    end
    Qhi(t)  = qhi;
    Qlo(t)  = qlo;
       
    t       = 1 + t;
     
end
 
pol_1       = zeros(1,gs);
pol_2       = pol_1;
for i = 1:gs
    pol_1(1,i) = a_grid(1,pos(i,1));
    pol_2(1,i) = a_grid(1,pos(i,2));
end
 
Dist        = ( dist(1:gs,1) )' + ( dist(gs+1:2*gs,1) )';
 
close all
figure(1)
plot(a_grid,pol_1,'-',a_grid,pol_2,'-',a_grid,a_grid,'--')
title('Policy Function for Credit Holding: Economy with sigma = 1.5')
xlabel('Credit Stock Today')
ylabel('Credit Stock Tomorrow')
 
qstar=Q(1,25)
 
NBL=1/(1-qstar) %natural borrowing limit

