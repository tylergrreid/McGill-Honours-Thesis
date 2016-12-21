function [R_chaser_Hill,V_chaser_Hill] = ECI2Hill_New(R_target_ECI,V_target_ECI,R_chaser_ECI,V_chaser_ECI,i,OMEGA,theta)
%% DESCRIPTION

% Given the position and velocity vectors of the target and chaser s/c's
% expressed in the ECI frame, output the postion and velocity of the chaser
% s/c in the Hill frame which is centered at the target

% Algorithm is based on the eci2hill algorithm found in 'Fundamentals of
% Astrodynamics and Applications 3rd ed.' by David A. Vallado (2007) 
% p.413-414

%% INPUT

% POSITION VECTORS HAVE UNITS OF [m]
% VELOCITY VECTORS HAVE UNITS OF [m/s]
% ANGULAR MEASUREMENTS HAVE UNITS OF [rad]

% R_target_ECI     - Position vector of the target s/c in the ECI frame
% V_target_ECI     - Velocity vector of the target s/c in the ECI frame
% R_chaser_ECI     - Position vector of the chaser s/c in the ECI frame
% V_chaser_ECI     - Velocity vector of the chaser s/c in the ECI frame
% i                - Inclination of the reference orbit
% OMEGA            - longitude of the ascending node of the reference orbit
% theta            - argument of latitude of the reference orbit


%% OUTPUT

% POSITION VECTORS HAVE UNITS OF [m]
% VELOCITY VECTORS HAVE UNITS OF [m/s]

% R_chaser_Hill    - Position vector of the chaser s/c in the Hill frame
% V_chaser_HIll    - Velocity vector of the chaser s/c in the Hill frame

%% NOTES

% (1) - ALL INPUT AND OUTPUT ARGUMENTS ARE COLUMN VECTORS 

%% Define Global Variables Used
global n

%% IMPLEMENTATION:
%% Determine the transformation matrix from Hill to ECI frame
trans_matrix(1,1) = cos(OMEGA)*cos(theta)-sin(OMEGA)*sin(theta)*cos(i);
trans_matrix(2,1) = sin(OMEGA)*cos(theta)+cos(OMEGA)*sin(theta)*cos(i);
trans_matrix(3,1) = sin(theta)*sin(i);
trans_matrix(1,2) = -cos(OMEGA)*sin(theta)-sin(OMEGA)*cos(theta)*cos(i);
trans_matrix(2,2) = -sin(OMEGA)*sin(theta)+cos(OMEGA)*cos(theta)*cos(i);
trans_matrix(3,2) = cos(theta)*sin(i);
trans_matrix(1,3) = sin(OMEGA)*sin(i);
trans_matrix(2,3) = -cos(OMEGA)*sin(i);
trans_matrix(3,3) = cos(i);

%% Compute Hill Position Vector of The Spacecraft wrt the Reference Orbit

R_chaser_Hill = trans_matrix'*(R_chaser_ECI - R_target_ECI);

%% Compute Hill Velocity Vector of The Spacecraft wrt the Reference Orbit

vect1(1,1) = n*R_chaser_Hill(2);
vect1(2,1) = -n*R_chaser_Hill(1);
vect1(3,1) = 0;

V_chaser_Hill = trans_matrix'*V_chaser_ECI - trans_matrix'*V_target_ECI + vect1;
