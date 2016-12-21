function out = Equations_of_Motion(t,y)
%% DESCIPTION 

% Given a time (t) and a vector (y) this function will output the necessary
% vector needed to implement the ODE45 - adaptive Runge-Kutta 4th/5th order
% numerical integrator.

% Equations of motion are based on those presented by Schweighart and
% Sedwick in the article 'High-Fidelity Linearized J_2 Model for Satellite
% Formation Flight' - Journal of Guidance, Control, and Dynamics Vol. 25,
% No. 6, November-December 2002.

% The equations used here are modified to accept dimensionless position and
% velocity parameters for numerical stability.  

%% INPUT

% t     - Current dimensionless time tau in the integration [-]

% y     - Current dimensionless positions and velocity parameters of the 
%         chief and the deputy s/c

%% OUPUT

% out   - The output in this case is a dimensionless vector of
%         differentials/slopes for the ode45 algorithm to use in the next
%         of numerical integration.

%% DECLARE GLOBAL VARIABLES USED
global R_e J_2 omega_e c n_0 n k r_ref i_ref e_ref q l phi DRAG Area_Chief Area_Deputy d tau_0

%% APPLY ORBIT RE-SET FILTER

% Check that individual motion from reference orbit Hill frame has not 
% drifted further than allowable tolerance.  This will re-set the initial
% conditions and reference orbit parameters and constants.

% y = ref_orbit_reset_filter(t,y);

%% COMPUTE THE CURRENT ORBITAL PARAMETERS OF THE REFERENCE ORBIT 

% compute the orbital parameters of the ref. circ. orbit
[i_current,OMEGA_current,theta_current] = ref_orbit_classical_elements(t);

%% FOR SIMPLICITY RE-DEFINE THE VARIABLES THAT THE INPUT VECTOR CORRESPONDS TO

xh_chief = y(1);
xh_prime_chief = y(2);
yh_chief = y(3);
yh_prime_chief = y(4);
zh_chief = y(5);
zh_prime_chief = y(6);

xh_deputy = y(7);
xh_prime_deputy = y(8);
yh_deputy = y(9);
yh_prime_deputy = y(10);
zh_deputy = y(11);
zh_prime_deputy = y(12);

%% COMPUTE DRAG EFFECTS
%% STEP 1 - DETERMINE THE POS. AND VEL. OF THE REF. CIR. Hill FRAME IN ECI FRAME
% THESE VARIABLES WILL HAVE UNITS OF [m] & [m/s]

[R_Hill,V_Hill] = ref_orbit_classical_elements2ECI(r_ref,e_ref,i_current,OMEGA_current,0,theta_current);

%% STEP 2 - DETERMINE THE FINAL DIMESIONLESS DRAG EFFECTS FOR THE CHIEF

% Determine the dimensionless velocity of the chief relative to the 
% rotating atmosphere expressed in the Hill frame

v_rel_chief(1) = xh_prime_chief-yh_chief*(n*c/n_0-omega_e*cos(i_current)/n_0)-...
    zh_chief*omega_e*cos(theta_current)*sin(i_current)/n_0;
v_rel_chief(2) = yh_prime_chief+(xh_chief+r_ref/d)*(n*c/n_0-omega_e*cos(i_current)/n_0)+...
    zh_chief*sin(i_current)*omega_e*sin(theta_current)/n_0;
v_rel_chief(3) = zh_prime_chief-yh_chief*sin(i_current)*omega_e*sin(theta_current)/n_0+...
    (r_ref/d+xh_chief)*omega_e*cos(theta_current)*sin(i_current)/n_0;

% Define the position and velocity of the CHIEF in the Hill frame [m]&[m/s]
R_chief_Hill = d*[xh_chief yh_chief zh_chief];
V_chief_Hill = d*n_0*[xh_prime_chief yh_prime_chief zh_prime_chief];

% Determine the position and velocity of the CHIEF in the ECI frame
[R_chief_ECI,V_chief_ECI] = Hill2ECI_New(R_Hill',V_Hill',R_chief_Hill',V_chief_Hill',0,i_current,OMEGA_current,theta_current);

% Determine drag effects for the CHIEF
f_drag_chief = relative_drag_perturbation(R_chief_ECI,v_rel_chief,Area_Chief);

%% STEP 2A - DETERMINE THE FINAL DIMESIONLESS DRAG EFFECTS FOR THE DEPUTY 

% Determine the dimensionless velocity of the deputy relative to the 
% rotating atmosphere expressed in the Hill frame
v_rel_deputy(1) = xh_prime_deputy-yh_deputy*(n*c/n_0-omega_e*cos(i_current)/n_0)-...
    zh_deputy*omega_e*cos(theta_current)*sin(i_current)/n_0;
v_rel_deputy(2) = yh_prime_deputy+(xh_deputy+r_ref/d)*(n*c/n_0-omega_e*cos(i_current)/n_0)+...
    zh_deputy*sin(i_current)*sin(theta_current)*omega_e/n_0;
v_rel_deputy(3) = zh_prime_deputy-yh_deputy*sin(i_current)*sin(theta_current)*omega_e/n_0+...
    (r_ref/d+xh_deputy)*omega_e*cos(theta_current)*sin(i_current)/n_0;

% Define the position and velocity of the DEPUTY in the Hill frame [m]&[m/s]
R_deputy_Hill = d*[xh_deputy yh_deputy zh_deputy];
V_deputy_Hill = d*n_0*[xh_prime_deputy yh_prime_deputy zh_prime_deputy];

% Determine the position and velocity of the DEPUTY in the ECI frame
[R_deputy_ECI,V_deputy_ECI] = Hill2ECI_New(R_Hill',V_Hill',R_deputy_Hill',V_deputy_Hill',0,i_current,OMEGA_current,theta_current);

% Determine drag effects for the DEPUTY
f_drag_deputy = relative_drag_perturbation(R_deputy_ECI,v_rel_deputy,Area_Deputy);

%% EQUATIONS OF MOTION FOR THE CHIEF

% x_chief component
out(1,1) = xh_prime_chief;
out(2,1) = 2*(n/n_0)*c*yh_prime_chief+(5*c^2-2)*(n/n_0)^2*xh_chief...
    -3*(n/n_0)^2*J_2*(R_e^2/(r_ref*d))*(0.5-(1.5*sin(i_ref)^2*sin(k*(t-tau_0)/n_0)^2)-...
    (1+3*cos(2*i_ref))/8) + DRAG*f_drag_chief(1);

% y_chief component
out(3,1) = yh_prime_chief;
out(4,1) = -2*(n/n_0)*c*xh_prime_chief-...
    3*(n/n_0)^2*J_2*(R_e^2/(r_ref*d))*sin(i_ref)^2*sin(k*(t-tau_0)/n_0)*cos(k*(t-tau_0)/n_0)...
    + DRAG*f_drag_chief(2);

% z_chief component
out(5,1) = zh_prime_chief;
out(6,1) = -(q/n_0)^2*zh_chief + 2*(l*q/(d*n_0^2))*cos(q*(t-tau_0)/n_0+phi) + DRAG*f_drag_chief(3);

%% EQUATIONS OF MOTION FOR THE DEPUTY

% x_deputy component
out(7,1) = xh_prime_deputy;
out(8,1) = 2*(n/n_0)*c*yh_prime_deputy+(5*c^2-2)*(n/n_0)^2*xh_deputy...
    -3*(n/n_0)^2*J_2*(R_e^2/(r_ref*d))*(0.5-(1.5*sin(i_ref)^2*sin(k*(t-tau_0)/n_0)^2)-...
    (1+3*cos(2*i_ref))/8) + DRAG*f_drag_deputy(1);

% y_deputy component
out(9,1) = yh_prime_deputy;
out(10,1) = -2*(n/n_0)*c*xh_prime_deputy-...
    3*(n/n_0)^2*J_2*(R_e^2/(r_ref*d))*(sin(i_ref)^2)*sin(k*(t-tau_0)/n_0)*cos(k*(t-tau_0)/n_0)...
    + DRAG*f_drag_deputy(2);

% z_deputy component
out(11,1) = zh_prime_deputy;
out(12,1) = -(q/n_0)^2*zh_deputy + 2*(l*q/(d*n_0^2))*cos(q*(t-tau_0)/n_0+phi) + DRAG*f_drag_deputy(3);
