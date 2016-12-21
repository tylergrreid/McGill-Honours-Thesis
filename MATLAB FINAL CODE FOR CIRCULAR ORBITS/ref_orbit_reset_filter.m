function out = ref_orbit_reset_filter(t,y,i_current,OMEGA_current,theta_current)
%% DESCRIPTION

%% INPUT

%% OUPUT

%% NOTES

%% IMPLEMENTATION:
%% DEFINE GLOBAL VARIABLES TO BE USED

% this parameters maybe re-defined based on new reference orbit if
% reference orbit re-set criteria is met
global mu J_2 R_e d r_ref e_ref i_ref OMEGA_ref theta_ref s c n n_0 k l q tau_0 ref_orbit_tol

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

%% APPLY RE-SET CRITERION
pos = d*[xh_chief yh_chief zh_chief];
dist = norm(pos);

if dist/r_ref < ref_orbit_tol
    % no change, do not alter global variables
    % do not alter IC's or ref. orbit
    out = y;
    
elseif dist >= ref_orbit_tol
    %% compute current reference orbital parameters (classical elements)
    [i_current,OMEGA_current,theta_current] = ref_orbit_classical_elements(t);
    
    %% compute reference orbit ECI position and velocity vectors
    [R_Hill,V_Hill] = ref_orbit_classical_elements2ECI(r_ref,e_ref,i_current,OMEGA_current,0,theta_current);
    
    %% re-set global tau_0
    tau_0 = t;
    % this essentially re-sets the time in the integratation for the
    % equations of motion and in computing the classical orbital elements
    % of the reference circular orbit
    
    %% compute the ECI position and velocity vectors of the CHIEF
    
    % Define the position and velocity of the CHIEF in the Hill frame [m]&[m/s]
    R_chief_Hill = d*[xh_chief yh_chief zh_chief];
    V_chief_Hill = d*n_0*[xh_prime_chief yh_prime_chief zh_prime_chief];
    
    % Compute the position and velocity of the CHIEF in the ECI frame [m]&[m/s]
    [R_chief_ECI,V_chief_ECI] = Hill2ECI(R_Hill',V_Hill',R_chief_Hill',V_chief_Hill',1);
    
    %% compute the ECI position and velocity vectors of the DEPUTY
    
    % Define the position and velocity of the DEPUTY in the Hill frame [m]&[m/s]
    R_deputy_Hill = d*[xh_deputy yh_deputy zh_deputy];
    V_deputy_Hill = d*n_0*[xh_prime_deputy yh_prime_deputy zh_prime_deputy];
    
    % Determine the position and velocity of the DEPUTY in the ECI frame
    [R_deputy_ECI,V_deputy_ECI] = Hill2ECI(R_Hill',V_Hill',R_deputy_Hill',V_deputy_Hill',1);
    
    %% re-define reference orbit parameters
    
    % radius of the reference orbit
    r_ref = norm(R_chief_ECI);
    
    % inclination of the reference orbit
    i_ref = i_current;
    
    % RAAN of the reference orbit
    OMEGA_ref = OMEGA_current;
    
    % determine the new argument of latitude based on ECI position of the
    % Chief spacecraft [this is the new centering point]
    theta_ref = theta_reset(r_ref,i_ref,OMEGA_ref,R_chief_ECI);
    
    % based on these elements - determine the ECI position and velocity of
    % the newly defined Hill frame 
    [R_Hill_new,V_Hill_new] = ref_orbit_classical_elements2ECI(r_ref,e_ref,i_current,OMEGA_current,0,theta_current);
    
    %% update s,n,c,k
    
    s = 3*J_2*(R_e^2)*(1+3*cos(2*i_ref))/(8*r_ref^2);
    c = sqrt(1+s);
    n = sqrt(mu/(r_ref^3));
    k = n*c + 3*n*J_2*(R_e^2)*(cos(i_ref)^2)/(2*(r_ref^2));
    
    %% redefine ICs of CHIEF
    
    % given the ECI state vectors of the chief and the reference Hill frame
    % output the state vector of the chief expressed in the newly defined
    % Hill frame
    [R_chief_Hill_New,V_chief_Hill_New] = ECI2Hill(R_Hill_new',V_Hill_new',R_chief_ECI,V_chief_ECI);
    
    %% redefine ICs of DEPUTY
    
    % given the ECI state vectors of the chief and the reference Hill frame
    % output the state vector of the chief expressed in the newly defined
    % Hill frame
    [R_deputy_Hill_New,V_deputy_Hill_New] = ECI2Hill(R_Hill_new',V_Hill_new',R_deputy_ECI,V_deputy_ECI);
    
    %% for simplicity - specify which variable is which
    
    x_chief = R_chief_Hill_New(1);
    y_chief = R_chief_Hill_New(2);
    z_chief = R_chief_Hill_New(3);
    x_dot_chief = V_chief_Hill_New(1);
    y_dot_chief = V_chief_Hill_New(2);
    z_dot_chief = V_chief_Hill_New(3);

    x_deputy = R_deputy_Hill_New(1);
    y_deputy = R_deputy_Hill_New(2);
    z_deputy = R_deputy_Hill_New(3);
    x_dot_deputy = V_deputy_Hill_New(1);
    y_dot_deputy = V_deputy_Hill_New(2);
    z_dot_deputy = V_deputy_Hill_New(3);

    %% update remaining constants used in EOM l,q
       
    i_sat2 = i_ref;
    i_sat1 = i_sat2 + (z_dot_deputy - z_dot_chief)/(k*r_ref);
    delta_OMEGA_0 = (z_deputy - z_chief)/(r_ref*sin(i_ref));
    gamma_0 = acot((cot(i_sat2)*sin(i_sat1)-...
        cos(i_sat1)*cos(delta_OMEGA_0))/sin(delta_OMEGA_0));
    PHI_0 = acos(cos(i_sat1)*cos(i_sat2) + sin(i_sat1)*sin(i_sat2)*cos(delta_OMEGA_0));
    OMEGA_sat1 = -3*n*J_2*(R_e^2)*cos(i_sat1)/(2*(r_ref^2));
    OMEGA_sat2 = -3*n*J_2*(R_e^2)*cos(i_sat2)/(2*(r_ref^2));
    q = n*c - (cos(gamma_0)*sin(gamma_0)*cot(delta_OMEGA_0)-...
        (sin(gamma_0)^2)*cos(i_sat1))*(OMEGA_sat1-OMEGA_sat2) - OMEGA_sat1*cos(i_sat1);
    l = -r_ref*sin(i_sat1)*sin(i_sat2)*sin(delta_OMEGA_0)*(OMEGA_sat1-OMEGA_sat2)/sin(PHI_0);
    
    %% Non-dimensionalize the ICs of CHIEF & DEPUTY and OUTPUT to integrator
    
    % non-dimensionalze the chief & deputy IC's
    xh_chief = x_chief/d;
    yh_chief = y_chief/d;
    zh_chief = y_chief/d;
    xh_prime_chief = x_dot_chief/(d*n_0);
    yh_prime_chief = y_dot_chief/(d*n_0);
    zh_prime_chief = z_dot_chief/(d*n_0);
    
    xh_deputy = x_deputy/d;
    yh_deputy = y_deputy/d;
    zh_deputy = z_deputy/d;
    xh_prime_deputy = x_dot_deputy/(d*n_0);
    yh_prime_deputy = y_dot_deputy/(d*n_0);
    zh_prime_deputy = z_dot_deputy/(d*n_0);
    
    % output in the correct order for integrator
    out(1) = xh_chief;
    out(2) = xh_prime_chief;
    out(3) = yh_chief;
    out(4) = yh_prime_chief;
    out(5) = zh_chief;
    out(6) = zh_prime_chief;
    out(7) = xh_deputy;
    out(8) = xh_prime_deputy;
    out(9) = yh_deputy;
    out(10) = yh_prime_deputy;
    out(11) = zh_deputy;
    out(12) = zh_prime_deputy;
    
else
    printf('There is an ERROR in ref_orbit_reset_filter.m\n');
end
