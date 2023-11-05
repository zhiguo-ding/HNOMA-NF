%courtesy to Zhaolin Wang @ QMUL
clear all
%close all

N =513; % number of antennas
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % wavelength
d = lambda / 2; % antenna spacing
D = (N-1)*d; %apture or size
dis_Ray = 2*D^2/lambda;
dis_Fre = sqrt(D^3/lambda)/2;
noise = 10^(-100/10); %noise power
P =0.01;% near field users' power
Roma = 5;
Rt=2;
     
%generate the channel vector  for near field users
M = 3;
r =[ 40  10 5 ];%[  5+20*[0:1:M-1]]; % distance
r = r(1:M); 
%%%%%%%%%%%%%%%%%%%%%%%%
%r(end)=10000000;
%%%%%%%%%%%%%%%%%%%%%%%

theta = ones(M,1)*45 / 180 * pi;%[0:90/M:90] / 180 * pi;%45* ones(M,1)/180 * pi;%%90*rand(nnf,1) / 180 * pi; % direction
for i = 1 : M
    w(:,i) = beamfocusing(r(i), theta(i), N, d, f)/sqrt(N); % beamforming vector
    alpham(i) = lambda/4/pi/r(i);
    h(:,i) = sqrt(N)*alpham(i)*w(:,i);
end


 
%generate the channel vector  for the single far field user
rff = 100;%dis_Ray+10; % distance
thetaff = 45 / 180 * pi;%90*rand(nnf,1) / 180 * pi; % direction
w0 = beamfocusing(rff, thetaff, N, d, f)/sqrt(N); % beamforming vector
alpha0 = lambda/4/pi/rff;
h0 = sqrt(N)*alpha0*w0;  

%simplified channel gains
gmx = (abs(h0'*w))'; % M gm
g0 = abs(h0'*w0)^2; 
hmx = diag(abs(h'*w).^2);% M diagoal ones for hm
temp1 = abs(h'*w).^2;
amx = 1./(P*sum(temp1,2) + noise);

%order the near field users according to hm
[z,hm_indx]=sort(hmx,'descend');
gm = gmx(hm_indx,1);
h = h(:,hm_indx);
w = w(:,hm_indx);
hm = hmx(hm_indx);
am =amx(hm_indx);

b0 = g0/noise;
a0 =  1./(P*sum(abs(h0'*w).^2) + noise);

A = []; % No other constraints
b = [];
Aeq = [];%
beq = [];%zeros(M+1,1); 
lb = [];
ub = [];
x0 =  ones(M+1,1); %initilization 
options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
x1 = fmincon(@(x) sum(M*x(1:M))+x(M+1),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,a0,gm,Roma,b0,P,noise,M,Rt,hm),options)
           

x2 = fmincon(@(x) sum(M*x(1:M))+x(M+1),x0,A,b,Aeq,beq,lb,ub,@(x) mycons2(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,am),options)
   
% Mx = min(M-1,M);
% x0 =  ones(Mx+1,1); %initilization 
% x3 = fmincon(@(x) sum(Mx*x(1:Mx))+x(Mx+1),x0,A,b,Aeq,beq,lb,ub,@(x) mycons3(x,a0,gm(1:Mx),Roma,b0,P,noise,Mx,Rt,hm(1:Mx),am(1:Mx)),options)

Poma = (exp(Roma)-1)/b0

%analytical results
c = min(1/P, a0*sum(gm)^2);
lamx = (M^M*exp(Roma)/(c^M*b0))^(1/(M+1));
xana = ones(M,1)*(lamx/M-1/c);
xana = [xana; lamx-1/b0]

if M==1
    beta = min(am*hm,a0*gm^2);
    P1o = sqrt(exp(Roma)/beta/b0)-1/beta;
    P0o = sqrt(exp(Roma)/beta/b0)-1/b0;
    [[P1o; P0o] x2]
end
 
%total energy consumption 
[sum(M*x1(1:M))+x1(M+1) sum(M*x2(1:M))+x2(M+1) Poma]

 
%actual data rates achieved after taking the interference into
%consideration
for m = 1 : M
    temp1 = abs(h0'*w); % |h_0^Hw_m|
    fm = h0'*w./(temp1);
    R1m(m) = log(1+abs(sum(h(:,m)'*w.*fm.*sqrt(x2(1:M)')))^2.*am(m));
end
R1m = [R1m log(1+a0*( sum(gm.*sqrt(x2(1:M))) )^2)]';
R1mapp = [log(1+am.*hm.*x2(1:M)) ;log(1+a0*( sum(gm.*sqrt(x2(1:M))) )^2)];
 [(R1m) * M + log(1+x2(end)*b0)   R1mapp*M+ log(1+x2(end)*b0)]

function [c,ceq] = mycons(x,a0,gm,Roma,b0,P,noise,M,Rt,hm)

temp1 = log(1+a0*(sum(gm.* sqrt(x(1:M)) ))^2);
c(1,1) = Roma - M*temp1-log(1+b0*x(M+1)); 
c(2:M+1,1) = x(1:M)-P/(exp(Rt)-1)+noise./hm;%Rt-log(1+P*hm./(hm.*x(1:M)+noise));
c(M+2:2*M+2,1) = -x;
ceq = [];
end

function [c,ceq] = mycons2(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,am)

temp1 = log(1+a0*(sum(gm.* sqrt(x(1:M)) ))^2);
c(1,1) = Roma - M*temp1-log(1+b0*x(M+1)); 
c(2:M+1,1) = Roma - M*log(1+am.*hm.*x(1:M)) - log(1+b0*x(M+1)); 
c(M+2:2*M+2,1) = -x;
ceq = [];
end

function [c,ceq] = mycons3(x,a0,gm,Roma,b0,P,noise,Mx,Rt,hm,am)

temp1 = log(1+a0*(sum(gm.* sqrt(x(1:Mx)) ))^2);
c(1,1) = Roma - Mx*temp1-log(1+b0*x(Mx+1)); 
c(2:Mx+1,1) = Roma - Mx*log(1+am.*hm.*x(1:Mx)) - log(1+b0*x(Mx+1)); 
c(Mx+2:2*Mx+2,1) = -x;
ceq = [];
end

 

