function [R_DEPUTY] = Hill2ECI(R_CHIEF,r_rel_DUPUTY,i,OMEGA,theta)
%% DESCRIPTION

%% IMPLEMENTATION:
%% Determine the Direction Vectors of the Hill (r-s-w) Frame in the ECI
%% Frame (X-Y-Z) and the Transformation Between Them

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

%% Compute ECI Position Vector of the Chaser s/c 
% (i.e. chief or deputy w.r.t. the reference orbit)

% transform target s/c position vector to r-s-w (Hill) frame
R_target_RSW = trans_matrix'*R_target_ECI;

% fix x
r_chaser_temp = R_target_RSW;
r_chaser_temp(1) = r_chaser_temp(1)+R_chaser_Hill(1);

% Perform rotation to fix y
angle_y = R_chaser_Hill(2)/norm(R_target_ECI);
rot3 = [cos(-angle_y) sin(-angle_y) 0; -sin(-angle_y) cos(-angle_y) 0; 0 0 1];
r_chaser_temp = rot3*r_chaser_temp;

% Perform rotation to fix z
angle_z = R_chaser_Hill(3)/norm(R_target_ECI);
rot2 = [cos(angle_z) 0 -sin(angle_z); 0 1 0; sin(angle_z) 0 cos(angle_z)];
r_chaser_temp = rot2*r_chaser_temp;

% Transform to ECI Coordinates to Obtain Final Answer
R_chaser_ECI = trans_matrix*r_chaser_temp;

%% Compute ECI Velocity Vector of the Chaser s/c

if CALC == 1
    % transform target s/c velocity vector to r-s-w (Hill) frame
    V_target_RSW = trans_matrix'*V_target_ECI;
    
    % compute third term needed for velocity (refer to algorithm by Vallado)
    temp1 = cross(R_target_RSW,V_target_RSW)/(norm(R_target_RSW)^2);
    temp2 = cross(temp1,R_chaser_Hill);
    
    v_chaser_temp = V_target_RSW + V_chaser_Hill + temp2;
    
    V_chaser_ECI = trans_matrix*v_chaser_temp;
    
else
    % if the velocity does not need to be calculated simply set it to 0
    % as it will not be needed or used in subsequent calculations
    V_chaser_ECI = 0;
end
