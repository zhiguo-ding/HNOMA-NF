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
    
% r = 10;%[5:5:20]; % distance
% theta = [25 30 35 40 45 50 55 60 65] / 180 * pi; % direction
% nfx = length(theta);
% for i = 1 : nfx
%     w(:,i) = beamfocusing(r, theta(i), N, d, f)/sqrt(N); % beamforming vector
% end
  
nnf = 3;
r =[5 10 40];%[  ceil(dis_Fre)+10*[0:1:nnf-1]]; % distance
theta = 45 / 180 * pi;%90*rand(nnf,1) / 180 * pi; % direction
%nfx = length(r);
for i = 1 : nnf
    w(:,i) = beamfocusing(r(i), theta, N, d, f)/sqrt(N); % beamforming vector
end

% calculate beampattern
m = 500;
area =  50;%dis_Ray; 
area2=  0;%dis_Ray/2;
X = linspace(area2, area, m);
Y = linspace(area2, area, m);
[X, Y] = meshgrid(X, Y);
[theta_all, r_all] = cart2pol(X, Y);

P = zeros(length(r_all), length(theta_all));
P_ff = zeros(length(r_all), length(theta_all));
parfor i = 1:m
    for j = 1:m
        a = beamfocusing(r_all(i,j), theta_all(i,j), N, d, f)/sqrt(N);
        P(i,j) = sum(abs(a' * w))^2;
        if P(i,j)>1.2
             zzzz=a;
            %break
        end
    end
    
end
P(1,1)=1;
P_ff(1,1)=1;

% plot beampattern
[X, Y] = pol2cart(theta_all, r_all);
figure; hold on; colormap jet; colorbar;
mesh(Y,X,P); 
view([90,-90]);
xlim([area2,area]); ylim([area2,area]);
xlabel('x (m)'); ylabel('y (m)');
title(['NF',' ,M =',num2str(N), ' ,r = ',num2str(r),' ,Rayleigh distance is ',num2str(dis_Ray), ])

% % plot beampattern FF
% [X, Y] = pol2cart(theta_all, r_all);
% figure; hold on; colormap jet; colorbar;
% mesh(X,Y,P_ff); 
% view([90,-90]);
% xlim([area2,area]); ylim([area2,area]);
% xlabel('x (m)'); ylabel('y (m)');
% title(['FF',' ,M =',num2str(N),' ,r = ',num2str(r),' ,Rayleigh distance is ',num2str(dis_Ray), ])

