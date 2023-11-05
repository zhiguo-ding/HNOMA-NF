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
radius2 = 200;
radius1 = 150;% the radius of the ring for the far field users
radius3 = 10;
%Roma = 8;
Rt=4;
ct=1000;
     
M = 40;
Mx = 10;
% r =10*ones(M,1);%[  5+20*[0:1:M-1]]; % distance
% 
% theta =  [-90:180/M:90] / 180 * pi;%45* ones(M,1)/180 * pi;%%90*rand(nnf,1) / 180 * pi; % direction


Nbit =  [1:1:10];
for in = 1 : length(Nbit)
    Roma = Nbit(in);
    temp1 =0;temp2 =0;temp3 =0;
    for ic = 1 : ct
        %channel
        %the location of the near-field users [5 70]
        NF_loc=[];
        while size(NF_loc,1)<=M
            x_loc = [radius1*rand(1,1) sign(randn)*radius1*rand(1,1)];
            if sqrt(x_loc*x_loc')<radius3 & sqrt(x_loc*x_loc')>5
                NF_loc = [NF_loc; x_loc];
            end
        end           
        [theta, r] = cart2pol(NF_loc(:,1), NF_loc(:,2));
        for i = 1 : M
            wch(:,i) = beamfocusing(r(i), theta(i), N, d, f)/sqrt(N); % beamforming vector
            alpham(i) = lambda/4/pi/r(i);
            h(:,i) = sqrt(N)*alpham(i)*wch(:,i);
        end

        Dno = diag((diag(inv(h'*h)).^(-1/2)));
        w = h*inv(h'*h)*Dno;%/sqrt(M);         

        %the location of the far-field user [100 150]        
        FF_loc=[];
        while size(FF_loc,1)<=0
            x_loc = [radius2*rand(1,1) sign(randn)*radius2*rand(1,1)];
            if sqrt(x_loc*x_loc')<radius2 & sqrt(x_loc*x_loc')>radius1
                FF_loc = [FF_loc; x_loc];
            end
        end
        [thetaff, rff] = cart2pol(FF_loc(:,1), FF_loc(:,2));
        
        %generate the channel vector  for the single far field user
%         rff = 100;%dis_Ray+10; % distance
%         thetaff = 45 / 180 * pi;%90*rand(nnf,1) / 180 * pi; % direction
        w0 = beamfocusing(rff, thetaff, N, d, f)/sqrt(N); % beamforming vector
        alpha0 = lambda/4/pi/rff;
        h0 = sqrt(N)*alpha0*w0;  

        %simplified channel gains
        gmx = (abs(h0'*w))'; % M gm
        g0 = abs(h0'*w0)^2; 
        hmx = diag(abs(h'*w).^2);% M diagoal ones for hm
        temp1x = abs(h'*w).^2;
        amx = 1./(P*sum(temp1x,2) + noise);

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
        options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        %approach 1
        %approach 1, we need to check whether a beam is qualified
        test_hm = hm-noise*(exp(Rt)-1)/P;
        [index_bad] = find(test_hm<0);
        Mnew = M-length(index_bad);        x0new =  ones(Mnew+1,1);
        x1 = fmincon(@(x) sum(M*x(1:Mnew))+x(Mnew+1),x0new,A,b,Aeq,beq,lb,ub,@(x) mycons(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,Mnew),options);
        temp1 = temp1 + sum(M*x1(1:Mnew))+x1(Mnew+1);

        %approach 2        
        x0 =  ones(Mx+1,1); %initilization 
        x2 = fmincon(@(x) sum(M*x(1:Mx))+x(Mx+1),x0,A,b,Aeq,beq,lb,ub,@(x) mycons2(x,a0,gm,Roma,b0,P,noise,M,Rt,hm,am,Mx),options);
        temp2 = temp2+ sum(M*x2(1:Mx))+x2(Mx+1);

        %OMA
        Poma = (exp(Roma)-1)/b0;
        temp3 = temp3+Poma;
        
        %actual data rates achieved after taking the interference into
        %consideration
        R1m = [];
        for m = 1 : Mx
            temp1x = abs(h0'*w); % |h_0^Hw_m|
            fm = h0'*w./(temp1x);
            R1m(m) = log(1+abs(sum(h(:,m)'*w(:,1:Mx).*fm(1:Mx).*sqrt(x2(1:Mx)')))^2.*am(m));
        end
        R1m = [R1m log(1+a0*( sum(gm(1:Mx).*sqrt(x2(1:Mx))) )^2)]';
        rateactx(ic) = min((R1m) * M + log(1+x2(end)*b0));
        R1mapp = [log(1+am(1:Mx).*hm(1:Mx).*x2(1:Mx)) ;log(1+a0*( sum(gm(1:Mx).*sqrt(x2(1:Mx))) )^2)];
        rateappx(ic)= min(R1mapp)*M+ log(1+x2(end)*b0);

    end
    scheme1(in) = temp1/ct;
    scheme2(in) = temp2/ct;
    oma(in) = temp3/ct;
    rateact(in) = min(rateactx);
    rateapp(in) = min(rateappx);
    
    % c = min(1/P, a0*sum(gm)^2);
    % lamx = (M^M*exp(Roma)/(c^M*b0))^(1/(M+1));
    % xana = ones(M,1)*(lamx/M-1/c);
    % if lamx/M-1/c<0
    %     xana = zeros(M,1);
    %     xana = [xana; Poma]; %simple OMA
    % else
    %     xana = [xana; lamx-1/b0];
    % end
    % approx(in) = sum(M*xana(1:M))+xana(M+1);

end
plot(Nbit,  (oma), Nbit,  (scheme1), Nbit,  (scheme2))%,Nbit, approx)

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

 
 

