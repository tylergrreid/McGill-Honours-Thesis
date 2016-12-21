function [R_chaser_Hill,V_chaser_Hill] = ECI2Hill(R_target_ECI,V_target_ECI,R_chaser_ECI,V_chaser_ECI)
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

% R_target_ECI     - Position vector of the target s/c in the ECI frame
% V_target_ECI     - Velocity vector of the target s/c in the ECI frame
% R_chaser_ECI     - Position vector of the chaser s/c in the ECI frame
% V_chaser_ECI     - Velocity vector of the chaser s/c in the ECI frame


%% OUTPUT

% POSITION VECTORS HAVE UNITS OF [m]
% VELOCITY VECTORS HAVE UNITS OF [m/s]

% R_chaser_Hill    - Position vector of the chaser s/c in the Hill frame
% V_chaser_HIll    - Velocity vector of the chaser s/c in the Hill frame

%% NOTES

% (1) - ALL INPUT AND OUTPUT ARGUMENTS ARE COLUMN VECTORS 
% (2) - THIS ROUTINE MAKES USE OF THE r-s-w NOTATION USED IN VALLADO (2007)
%       THIS FRAME IS EQUIVALENT TO THE HILL FRAME WHERE r=x, s=y, w=z
% (3) - 'TARGET' IS SYNONOMOUS WITH 'CHIEF' AND 'CHASER' IS SYNONOMOUS 
%       WITH 'DEPUTY'

%% IMPLEMENTATION:

%% Determine the orientation of the Hill (r-s-w) frame & Transformation Matrix

% r-vector
r_vect = R_target_ECI/norm(R_target_ECI);

% w-vector (same direction as angular momentum vector h)
w_vect = cross(R_target_ECI,V_target_ECI);
w_vect = w_vect/norm(w_vect);

% s-vector (complete the triad)
s_vect = cross(w_vect,r_vect);
s_vect = s_vect/norm(s_vect);

% Define the transformation matrix from ECI(X-Y-Z) to r-s-w frame
trans_matrix(1,1) = r_vect(1);
trans_matrix(2,1) = r_vect(2);
trans_matrix(3,1) = r_vect(3);
trans_matrix(1,2) = s_vect(1);
trans_matrix(2,2) = s_vect(2);
trans_matrix(3,2) = s_vect(3);
trans_matrix(1,3) = w_vect(1);
trans_matrix(2,3) = w_vect(2);
trans_matrix(3,3) = w_vect(3);

%% Compute The Position of The Chaser wrt The Target in The Hill Frame

% -----transform position vector of the target to the r-s-w frame-----
R_target_RSW = trans_matrix'*R_target_ECI;

% -----transform postion vector of the chaser to the r-s-w frane-----
R_chaser_RSW = trans_matrix'*R_chaser_ECI;

% -----compute y-offset to correct vector-----
angle_y = atan(R_chaser_RSW(2)/norm(R_target_RSW));
% determine rotation matrix about the 3rd axis ROT3
rot3 = [cos(angle_y) sin(angle_y) 0; -sin(angle_y) cos(angle_y) 0; 0 0 1];
% apply rotation matrix to the chaser
r_chaser_temp = rot3*R_chaser_RSW;

% -----compute z-offset to correct vector-----
angle_z = atan(R_chaser_RSW(3)/norm(R_target_RSW));
% determine rotation matrix
rot2 = [cos(-angle_z) 0 -sin(-angle_z); 0 1 0; sin(-angle_z) 0 cos(-angle_z)];
% apply rotation matrix to the chaser
r_chaser_temp = rot2*r_chaser_temp;

% -----output the position of the chaser in the Hill frame-----
R_chaser_Hill(1,1) = r_chaser_temp(1) - norm(R_target_RSW);
R_chaser_Hill(2,1) = angle_y*norm(R_target_RSW);
R_chaser_Hill(3,1) = -angle_z*norm(R_target_RSW);

%% Compute The Velocity of The Chaser wrt The Target in the Hill Frame

% -----transform the velocity vector of the target to the r-s-w frame-----
V_target_RSW = trans_matrix'*V_target_ECI;

% -----transform the velocity vector of the chaser to the r-s-w frame-----
V_chaser_RSW = trans_matrix'*V_chaser_ECI;

% -----output the velocity of the chaser wrt target in the Hill frame -----
temp1 = cross(R_target_RSW,V_target_RSW)/norm(R_target_RSW)^2;
V_chaser_Hill = V_chaser_RSW - V_target_RSW - cross(temp1,R_chaser_Hill);
