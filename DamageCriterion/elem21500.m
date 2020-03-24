%==========================================================================
%
%           This script plot the Mohr-Coulomb cercle on the plan (sigma-tau)
%
%==========================================================================

clear
clc
%% plot the failure criterion

D          = 0;
beta       = 0.1;
beta_angle = 1;
hr         = 0.33;            
c          = 2.8; %MPa
sigma_R    = 1.2; %MPa
phi        = 31/180*pi; % 

gD = (1-D).*(1-beta*log(1-D));
hD = hr+(1-D).^beta_angle*(1-hr);
XD = gD./hD;

tau_c = ( c^2+(sigma_R*tan(phi))^2 ) / (2*sigma_R*tan(phi));
b = ( c^2 - (sigma_R*tan(phi))^2 ) / (2*sigma_R*tan(phi));

tau_axe = [-5:0.01:5]*tau_c;
sigma_axe = [-5:0.01:1]*tau_c/tan(phi);
sigma_n = 1/hD/tan(phi)*(gD*tau_c-sqrt(tau_axe.^2+gD^2*b^2));

tau_asymp_positive = hD*tan(phi)*sigma_axe-gD*tau_c;
tau_asymp_negative = -hD*tan(phi)*sigma_axe+gD*tau_c;
sigma_asymp_negative = 1/hD/tan(phi)*(gD*tau_c-tau_axe);
figure(4);
clf;
hold on;
plot(sigma_n,tau_axe)
plot(sigma_axe,tau_asymp_positive,'--k',...
    sigma_axe,tau_asymp_negative,'--k');
%  xlim([-1,1]*tau_c*gD/hD/tan(phi))
daspect([1,1,1])
box on; grid on;
%% plot the Mohr Cirle

sigma1 = -20;  % MPa
sigma2 = -40;  % MPa

R = (sigma1 - sigma2)/2;
SigmaO = (sigma2 + sigma1)/2;

for i = 1:101
    alpha = pi/100*(i-1);
    sn(i) = SigmaO + R*cos(2*alpha);
    tn(i) = R*sin(2*alpha);
end
plot(sn,tn,'k')

%% plot the Mohr circle in effective stress
p_initial = 19;
p_max_approx = gD*tau_c/hD/tan(phi) - R/sin(phi)- SigmaO;
biot = 0.86;

for i = 1:101
    alpha = pi/100*(i-1);
    sn_eff_initial(i) = SigmaO + R*cos(2*alpha) + biot*p_initial;
    tn_eff_initial(i) = R*sin(2*alpha);
end

plot(sn_eff_initial,tn_eff_initial,'-g')



%% Analytical solution
k = hD*tan(phi);
x0 = 1/k*( gD*tau_c -sqrt( (R^2+gD^2*b^2)/(1+k^2) ) );
SO = (1+k^2)* x0 - k*gD*tau_c;
y0 = sqrt(  ( k*x0-gD*tau_c)^2 - gD^2*b^2  );
k0  = -(x0-SO)/y0;
k0  = k*(k*x0-gD*tau_c)/y0;
max_Pressure = (SO - SigmaO)/biot

y0^2 - (k*x0-gD*tau_c)^2 +gD^2*b^2;
y0^2 + (x0 -SO)^2;
plot(x0,y0,'*')   % tangent point
plot(SO,0,'*')    % center of the cercle tangent to the failure criteria
grid on;

for i = 1:1001
    alpha = pi/1000*(i-1);
    sn_eff_exact(i) = SigmaO + R*cos(2*alpha) + biot*max_Pressure;
    tn_eff_exact(i) = R*sin(2*alpha);
end

plot(sn_eff_exact,tn_eff_exact,'-')
plot(sigma_axe,k0*(sigma_axe-x0) + y0 )
 
 
 alpha_best = 1/2* atan(y0/(x0-SO))/pi*180;
 
%% stress on the fault

% Initial Total stress
Sn_fault = SigmaO + R*cos(2*alpha_best/180*pi)
Tn_fault = R*sin(2*alpha_best/180*pi)
plot(Sn_fault,Tn_fault,'o')
 
% Initial effective stress
Sn_eff_fault = Sn_fault + biot*p_initial;

plot(Sn_eff_fault,Tn_fault,'o')
 
% Activated effective stress
Sn_eff_fault = Sn_fault + biot*max_Pressure;

plot(Sn_eff_fault,Tn_fault,'o')
 