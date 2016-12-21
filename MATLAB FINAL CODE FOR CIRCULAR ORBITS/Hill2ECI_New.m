function [R_chaser_ECI,V_chaser_ECI] = Hill2ECI_New(R_target_ECI,V_target_ECI,R_chaser_Hill,V_chaser_Hill,CALC,i,OMEGA,theta)
%% DESCRIPTION

% Given the position vector of the target s/c in the ECI frame and the
% position vector of the chaser s/c in the Hill frame, this function will
% return the position of the chaser s/c in the ECI frame.

% Algorithm is based on the hill2eci algorithm found in 'Fundamentals of
% Astrodynamics and Applications 3rd ed.' by David A. Vallado (2007) 
% p.413-414

% refer to the prepared coordinate transformation document

%% INPUT

% POSITION VECTORS HAVE UNITS OF [m]
% VELOCITY VECTORS HAVE UNITS OF [m/s]
% ANGULAR MEASUREMENTS HAVE UNITS OF [rad]

% R_target_ECI     - Position vector of the target s/c in the ECI frame
% V_target_ECI     - Velocity vector of the target s/c in the ECI frame
% R_chaser_Hill    - Position vector of the chaser s/c in the Hill frame
% V_chaser_Hill    - Velocity vector of the chaser s/c in the Hill frame
% CALC             - This is a logic variable and can have an integer value
%                    of 0 or 1, if 0 - do not calculate the ECI Velocity 
%                    vector, if 1 - calculate the ECI Velocity vector, this
%                    is done in order to optimize run time and # of flops.
%                    It is a type of check bit.
% i                - Inclination of the Reference orbit
% OMEGA            - longitude of the ascending node of the reference orbit
% theta            - argument of latitude of the reference orbit

%% OUTPUT

% POSITION VECTORS HAVE UNITS OF [m]
% VELOCITY VECTORS HAVE UNITS OF [m/s]

% R_chaser_ECI     - Position vector of the chaser s/c in the ECI frame
% V_chaser_ECI     - Velocity vector of the chaser s/c in the ECI frame

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

%% Compute ECI Position Vector of the Spacecraft
% (i.e. chief or deputy w.r.t. the reference orbit)

R_chaser_ECI = R_target_ECI + trans_matrix*R_chaser_Hill;

%% Compute ECI Velocity Vector of the Spacecraft

% determine first if it needs to be calculated
if CALC == 1
    % define additional vector to be transformed
    vect(1,1) = V_chaser_Hill(1) - n*R_chaser_Hill(2);
    vect(2,1) = V_chaser_Hill(2) + n*R_chaser_Hill(1);
    vect(3,1) = V_chaser_Hill(3);
    
    % determine ECI velocity of the spacecraft
    V_chaser_ECI = V_target_ECI + trans_matrix*vect;
else
    % if the velocity does not need to be calculated simply set it to 0
    % as it will not be needed or used in subsequent calculations
    V_chaser_ECI = 0;
end
