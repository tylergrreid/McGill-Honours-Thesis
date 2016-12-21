function out = Equations_of_Motion(t,y)
%% DESCRIPTION
%% %% DEFINE CONSTANTS
global mu R_e omega_e J_2 B_inv_CHIEF B_inv_DEPUTY

% pass orbital elements (for ease of review)
a = y(1);
e = y(2);
i = y(3);
OMEGA = y(4);
omega = y(5);
M = y(6);

% pass relative position parameters (again, for ease of review)
x_deputy = y(7);
y_deputy = y(8);
z_deputy = y(9);
xdot_deputy = y(10);
ydot_deputy = y(11);
zdot_deputy = y(12);

% define the relative position parameters
xrel_DEPUTY = [y(7) y(8) y(9)]';

% do not allow eccentricity to go below 0 (this makes no sense physically)
if e<0
    e=0;
end

%% COMPUTE NEEDED PARAMETERS

% Use The input COE's to Compute the True Anomaly f
% Uses Vallado (2007) Algorithm 2, 6
f = True_Anomaly(M,e);

% mean orbital rate
n = sqrt(mu/a^3);

% compute the argument of latitude
theta = omega + f;

%% CALCULATIONS FOR THE CHIEF

% Calculate the Position Vector of The CHIEF s/c
X_CHIEF = COE2RV(a,e,i,OMEGA,omega,f);

% Calculate local atm. density rho above the ellipsoidal Earth for the CHIEF
rho_CHIEF = density_altitude_model(X_CHIEF);

% Compute the dimensionless ballistic coefficient for the CHIEF
beta_inv_CHIEF = (rho_CHIEF*a*B_inv_CHIEF*1000);

%% CALCULATIONS FOR THE DEPUTY

% Calculate the Position Vector of The DEPUTY s/c
% {This is based on the Linear Transformation of Appendix C}
X_DEPUTY = Hill2ECI_Linear(X_CHIEF,xrel_DEPUTY,i,OMEGA,theta);

% Calculate local atm. density rho above the ellipsoidal Earth for the DEPUTY
rho_DEPUTY = density_altitude_model(X_DEPUTY);

% Compute the dimensionless ballistic coefficient for the DEPUTY
beta_inv_DEPUTY = (rho_DEPUTY*a*B_inv_DEPUTY*1000);

%% OUTPUT FOR THE CHIEF

% a
out(1,1) = -(beta_inv_CHIEF)*a*(n-2*omega_e*cos(i));
% e
out(2,1) = -0.5*beta_inv_CHIEF*e*n;
% i
out(3,1) = -.25*beta_inv_CHIEF*omega_e*sin(i);
% OMEGA
out(4,1) = -1.5*J_2*n*(R_e/a)^2*cos(i);
% omega
out(5,1) = 0.75*J_2*n*(R_e/a)^2*(5*(cos(i)^2)-1);
% M
out(6,1) = n + 0.75*J_2*n*(R_e/a)^2*(3*(cos(i)^2)-1);

%% OUTPUT FOR THE DEPUTY
% [refer to state-space form Xdot = AX + W

% define variables
s = (3/8)*J_2*(R_e/a)^2*(1+3*cos(2*i));
stilda = (15/4)*J_2*(R_e/a)^2*sin(2*i);
sigma1 = 1-(omega_e/n)*cos(i)+(2/3)*e*cos(f);
sigma2 = 1-2*(omega_e/n)*cos(i)+(4/3)*e*cos(f);

% divided by n^2 already:
thetadotdot = (1.5*beta_inv_CHIEF*sigma2 - 2*e*sin(f)*(1+2*e*cos(f)-0.5*s));

% x-equation
out(7,1) = xdot_deputy;
% y-equation
out(8,1) = ydot_deputy;
% z-equation
out(9,1) = zdot_deputy;
% xdot-equation
out(10,1) = n*(-0.5*beta_inv_DEPUTY*sigma1*xdot_deputy + 2*(1+2*e*cos(f)+0.5*s)*ydot_deputy)...
    +(n^2)*((3+10*e*cos(f)+5*s)*x_deputy + (thetadotdot+0.5*beta_inv_DEPUTY*sigma2)*y_deputy...
    +(4*e*stilda*sin(omega)*z_deputy));
% ydot-equation
out(11,1) = n*(-2*(1+2*e*cos(f)+0.5*s)*xdot_deputy - beta_inv_DEPUTY*sigma1*ydot_deputy)...
    +(n^2)*((thetadotdot-beta_inv_DEPUTY*sigma2)*x_deputy + (e*cos(f))*y_deputy...
    -e*stilda*cos(omega)*z_deputy) - 0.5*(beta_inv_DEPUTY-beta_inv_CHIEF)*a*(n^2)*sigma2;
% zdot-equation
out(12,1) = n*(-0.5*beta_inv_DEPUTY*sigma1*zdot_deputy) + (n^2)*(4*e*stilda*sin(omega)*x_deputy...
    -e*stilda*cos(omega)*y_deputy-(1+3*e*cos(f)+3*s)*z_deputy);
