function out = vop_ode_ave(t,y)
%% DESCRIPTION
% This is the set of six 1st order odes which descibe the perturbations of
% the orbital elements.  

% In this case the time averaged rates of change of the the orbital 
% elements is considered.

%% DEFINE CONSTANTS
global mu R_e omega_e J_2 mass C_d Area rho_0

a = y(1);
i = y(2);
OMEGA = y(3);
e = y(4);
omega = y(5);
M = y(6);

%% COMPUTE NEEDED PARAMETERS
% Use The COE's to Compute the True Anomaly f
% Uses Vallado (2007) Algorithm 2, 6
f = True_Anomaly(M,e);

% mean motion
n = sqrt(mu/a^3);

% Calculate the Position Vector of The S/C
% X_pos = COE2RV(a,e,i,OMEGA,omega,f);

% Calculate local atm. density rho above the ellipsoidal Earth
rho = rho_0; %density_altitude_model(X_pos);

% Compute the dimensionless ballistic coefficient
beta_inv = (rho*1000*C_d*Area*a/mass);

%% COMPUTE DERIVATIVES
% a
out(1,1) = -(beta_inv)*a*n*(1-2*(omega_e/n)*cos(i));
% i
out(2,1) = -.25*beta_inv*omega_e*sin(i);
% OMEGA
out(3,1) = -1.5*J_2*n*(R_e/a)^2*cos(i);
% e
out(4,1) = -0.5*beta_inv*e*n;
% w
out(5,1) = 0.75*J_2*n*(R_e/a)^2*(5*cos(i)^2-1);
% M
out(6,1) = n+0.75*J_2*n*(R_e/a)^2*(3*cos(i)^2-1);