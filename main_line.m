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
%Roma = 8;
Rt=4;
     
%generate the channel vector  for near field users
M =3;
Mx = 3;
r =[40 10 5];%[  5+20*[0:1:M-1]]; % distance
r = r(1:M); 

theta = ones(M,1)*45 / 180 * pi;%[0:90/M:90] / 180 * pi;%45* ones(M,1)/180 * pi;%%90*rand(nnf,1) / 180 * pi; % direction
for i = 1 : M
    w(:,i) = beamfocusing(r(i), theta(i), N, d, f)/sqrt(N); % beamforming vector
    alpham(i) = lambda/4/pi/r(i);
    h(:,i) = sqrt(N)*alpham(i)*w(:,i);
end

%second choice for beamforming based on ZF
Dno = diag((diag(inv(h'*h)).^(-1/2)));
%w = h*inv(h'*h)*Dno;%/sqrt(M); 

%generate the channel vector  for the single far field user
rff = 200;%dis_Ray+10; % distance
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

Nbit =  [1:1:10];
for i = 1 : length(Nbit)
    Roma = Nbit(i);
    A = []; % No other constraints
    b = [];
    Aeq = [];%
    beq = [];%zeros(M+1,1); 
    lb = [];
    ub = [];
    x0 =  ones(M+1,1); %initilization 
    options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
    %approach 1
    %approach 1, we need to check whether a beam is qualified
        test_hm = hm-noise*(exp(Rt)-1)/P;
        [index_bad] = find(test_hm<0);
        Mnew = M-length(index_bad);        x0new =  ones(Mnew+1,1);
        x1 = fmincon(@(x) sum(M*x(1:Mnew))+x(Mnew+1),x0new,A,b,Aeq,beq,lb,ub,@(x) mycons(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,Mnew),options);
        scheme1(i) = sum(M*x1(1:Mnew))+x1(Mnew+1);

        %approach 2        
        x0 =  ones(Mx+1,1); %initilization 
        x2 = fmincon(@(x) sum(M*x(1:Mx))+x(Mx+1),x0,A,b,Aeq,beq,lb,ub,@(x) mycons2(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,am,Mx),options);
        scheme2(i) = sum(M*x2(1:Mx))+x2(Mx+1);

    %OMA
    Poma = (exp(Roma)-1)/b0;
    oma(i) = Poma;

    %analytical results
    c = min(1/P, a0*sum(gm(1:Mx))^2);
    lamx = (Mx^M*exp(Roma)/(c^M*b0))^(1/(M+1));
    xana = ones(Mx,1)*(lamx/Mx-1/c);
    if lamx/Mx-1/c<0
        xana = zeros(Mx,1);
        xana = [xana; Poma]; %simple OMA
    else
        xana = [xana; lamx-1/b0];
    end
    approx(i) = sum(M*xana(1:Mx))+xana(Mx+1);
    
    %analytical results for the 2-user special case
    if M==1
        beta = min(am*hm,a0*gm^2);
        P1o = sqrt(exp(Roma)/beta/b0)-1/beta;
        if P1o<0
            P1o = 0;
            P0o = Poma;%OMA
        else
            P0o = sqrt(exp(Roma)/beta/b0)-1/b0;
        end
        %[[P1o; P0o] x2];
        spe(i) = M*P1o+P0o;
    end

    %actual data rates achieved after taking the interference into
    %consideration, if beamfocusing is used as beamforming
    R1m = [];
    for m = 1 : Mx
        temp1x = abs(h0'*w); % |h_0^Hw_m|
        fm = h0'*w./(temp1x);
        R1m(m) = log(1+abs(sum(h(:,m)'*w(:,1:Mx).*fm(1:Mx).*sqrt(x2(1:Mx)')))^2.*am(m));
    end
    R1m = [R1m log(1+a0*( sum(gm(1:Mx).*sqrt(x2(1:Mx))) )^2)]';
    rateactx(i) = min((R1m) * M + log(1+x2(end)*b0));
    R1mapp = [log(1+am(1:Mx).*hm(1:Mx).*x2(1:Mx)) ;log(1+a0*( sum(gm(1:Mx).*sqrt(x2(1:Mx))) )^2)];
    rateappx(i)= min(R1mapp)*M+ log(1+x2(end)*b0);
end
%plot(Nbit,  (oma), Nbit,  (scheme1), Nbit,  (scheme2), Nbit, approx, Nbit,spe) 
plot(Nbit,  (oma), Nbit,  (scheme1), Nbit,  (scheme2))%, Nbit, approx) 
%plot(  Nbit,  (scheme2))
plot(  Nbit,  (scheme1), Nbit,  (scheme2))%, Nbit, approx) 

function [c,ceq] = mycons(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,Mnew)

temp1 = log(1+a0*(sum(gm(1:Mnew).* sqrt(x(1:Mnew)) ))^2);
c(1,1) = Roma - M*temp1-log(1+b0*x(Mnew+1)); %still M, not Mnew, as M time slots used
c(2:Mnew+1,1) = x(1:Mnew)-P/(exp(Rt)-1)+noise./hm(1:Mnew);%Rt-log(1+P*hm./(hm.*x(1:M)+noise));
c(Mnew+2:2*Mnew+2,1) = -x(1:Mnew+1);
ceq = [];
end

function [c,ceq] = mycons2(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,am,Mx)

temp1 = log(1+a0*(sum(gm(1:Mx).* sqrt(x(1:Mx)) ))^2);
c(1,1) = Roma - M*temp1-log(1+b0*x(Mx+1)); 
c(2:Mx+1,1) = Roma - M*log(1+am(1:Mx).*hm(1:Mx).*x(1:Mx)) - log(1+b0*x(Mx+1)); 
c(Mx+2:2*Mx+2,1) = -x;
ceq = [];
end

 

