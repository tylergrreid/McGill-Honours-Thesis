function rho = density_altitude_model(R)
%% DESCRIPTION

% Given the s/c position above the ellipsoidal earth in the ECI frame,
% output the density of the atmosphere at this position based on an
% exponential model and this altitude.

% Algorithm is based on:

% Position above elliposoidal earth:
% Based on Algorithm 12 of 'Fundamentals of Astrodynamics
% and Applications 3rd Ed. (2007)' by David A. Vallado


% Density altitude model:
% Refers to values taken from table 8.4 of 'Fundamentals of Astrodynamics
% and Applications 3rd Ed. (2007)' by David A. Vallado

%% INPUT

% R       - Position vector of the s/c in ECI expressed in ECI Frame.  [m]

%% OUPUT

% rho     - The density of the atmosphere at this inertial point.   [kg/m^3]

%% NOTES

% (1) - S/C IS SHORT FOR SPACECRAFT.
% (2) - PLEASE BE SURE THAT VARIABLES BEING PASSED ARE IN THE CORRECT
%       UNITS.

%% IMPLEMENTATION:
%% DEFINE GLOBAL VARIABLES TO BE USED

global R_e

%% DETERMINE THE ALTITUDE ABOVE THE ELLIPSOIDAL EARTH (h_ellp)

% Define r_I, r_J, r_K
r_I = R(1);
r_J = R(2);
r_K = R(3);

r_delta = sqrt((r_I^2)+(r_J^2));
 
% Eccentricity of ellipsoidal earth - eccentricity of wgs84 ellipsoid
e_E = 0.081819221456;

delta = atan(r_K/r_delta);

% Determine the Values of phi_gd and C_E:
phi = 1e5; %dummy value
phi_old = delta; %inital guess
tol = 1e-16; %tolerence to be used in the iteration

% Perform the iteration:
while abs(phi_old-phi)>tol
    C_E = R_e/sqrt(1-(e_E^2)*(sin(phi_old))^2);
    
    phi = atan((r_K+C_E*(e_E^2)*sin(phi_old))/r_delta);

    phi_old=phi;
end

% Using the values of phi_gd and C_E determined above, compute h_ellp [km]
h_ellp = (r_delta/cos(phi)-C_E);

%% DETERMINE THE CONSTANTS TO USE IN THE EXPONENTIAL MODEL [h_0,H,rho_0]

% if we select orbits in the range 150+ km altitudee, we need to
% select the values from Table A1: 

% NOTE: the values determined from this table are given in [km] this is why
% h_ellp had to be converted to [km] the values of density given are given
% in units of [kg/m^3] and the units of [km] cancel out in the final
% equation.

if 150<=h_ellp && h_ellp<180
    h_0 = 150;
    rho_0 = 2.070e-9;
    H = 22.523;
elseif 180<=h_ellp && h_ellp<200
    h_0 = 180;
    rho_0 = 5.464e-10;
    H = 29.740;
elseif 200<=h_ellp && h_ellp<250
    h_0 = 200;
    rho_0 = 2.789e-10;
    H = 37.105;
elseif 250<=h_ellp && h_ellp<300
    h_0 = 250;
    rho_0 = 7.248e-11;
    H = 45.546;
elseif 300<=h_ellp && h_ellp<350
    h_0 = 300;
    rho_0 = 2.418e-11;
    H = 53.628;
elseif 350<=h_ellp && h_ellp<400
    h_0 = 350;
    rho_0 = 9.158e-12;
    H = 53.298;
elseif 400<=h_ellp && h_ellp<450
    h_0 = 400;
    rho_0 = 3.725e-12;
    H = 58.515;
elseif 450<=h_ellp && h_ellp<500
    h_0 = 450;
    rho_0 = 1.585e-12;
    H = 60.828;
elseif 500<=h_ellp && h_ellp<600
    h_0 = 500;
    rho_0 = 6.967e-13;
    H = 63.822;
elseif 600<=h_ellp && h_ellp<700
    h_0 = 600;
    rho_0 = 1.454e-13;
    H = 71.835;
elseif 700<=h_ellp && h_ellp<800
    h_0 = 700;
    rho_0 = 3.614e-14;
    H = 88.667;
elseif 800<=h_ellp && h_ellp<900
    h_0 = 900;
    rho_0 = 1.170e-14;
    H = 124.64;
elseif 900<=h_ellp && h_ellp<1000
    h_0 = 900;
    rho_0 = 5.245e-15;
    H = 181.05;
elseif h_ellp>=1000
    h_0 = 1000;
    rho_0 = 3.019e-15;
    H = 268.00;
end
%% DETERMINE ATMOSPHERIC DENSITY RHO

rho = rho_0*exp(-(h_ellp-h_0)/H);
