function out = vop_ode(t,y)
%% DESCRIPTION
% This is the set of six 1st order odes which descibe the perturbations of
% the orbital elements.

% In this case we have included Full J2 and Full Dyanamic Drag.

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

% compute r
r = a*(1-e^2)/(1+e*cos(f));

% define arg of lat u,p,& h
u = omega+f;

p = a*(1-e^2);

h = sqrt(mu*a*(1-e^2));

n = sqrt(mu/(a^3));

%% COMPUTE J2 ACCELERATIONS

j2con = -1.5*J_2*mu*R_e^2/r^4;
a_j2 = j2con*[(1-3*sin(i)^2*sin(u)^2) sin(i)^2*sin(2*u) sin(2*i)*sin(u)];

%% COMPUTE DRAG ACCELERATIONS

% Calculate the Position Vector of The S/C
% X_pos = COE2RV(a,e,i,OMEGA,omega,f);

% Calculate local atm. density rho above the ellipsoidal Earth
rho = rho_0; %density_altitude_model(X_pos);

% ballistic coeff
B = mass/(C_d*Area);

v_rel(1) = (h/p)*e*sin(f);
v_rel(2) = r*((n*((1+e*cos(f))^2)/(1-e^2)^(3/2))-omega_e*cos(i));
v_rel(3) = r*omega_e*cos(u)*sin(i);

a_aero = -0.5*(rho*1000/B)*norm(v_rel)*v_rel;

%% COMBINE PERTURBATION ACCELERATIONS
a_pert = a_aero + a_j2;

%% COMPUTE DERIVATIVES

% a
out(1,1)=(2*a^2/h)*(e*sin(f)*a_pert(1) + (p/r)*a_pert(2));
% i
out(2,1) = r*cos(u)*a_pert(3)/h;
% OMEGA
out(3,1) = r*sin(u)*a_pert(3)/(h*sin(i));
% e
out(4,1) = (p*sin(f)*a_pert(1) + ((p+r)*cos(f)+r*e)*a_pert(2))/h;
% w
out(5,1) = (-p*cos(f)*a_pert(1) + (p+r)*sin(f)*a_pert(2))/(h*e) - r*sin(u)*cos(i)*a_pert(3)/(h*sin(i));
% M
out(6,1) = n + sqrt(1-e^2)*((p*cos(f)-2*r*e)*a_pert(1) - (p+r)*sin(f)*a_pert(2))/(h*e);