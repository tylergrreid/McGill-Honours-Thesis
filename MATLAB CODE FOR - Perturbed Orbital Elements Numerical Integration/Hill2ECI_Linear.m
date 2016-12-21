function [R_DEPUTY] = Hill2ECI_Linear(R_CHIEF,r_rel_DUPUTY,i,OMEGA,theta)
%% DESCRIPTION

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

R_DEPUTY = R_CHIEF + trans_matrix*r_rel_DUPUTY;

