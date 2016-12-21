function theta = theta_reset(r,i,OMEGA,R_ECI)
%% DESCRIPTION

% This algorithm is specific to this simulation and is not intended to be a
% general agorithm.  

% Given the radius of an INCLINED CIRCULAR orbit as well as the
% inclination, RAAN, and the ECI position vector, output the value of the
% argument of latitude 'u' or 'theta'.

% Based in part by algorithm 10 'COE2RV' by Vallado (2007)
% from 'Fundamentals of Astrodynamics and Applications 3rd Ed.'

%% INPUT

% r           - Radius of the orbit               [m]
% i           - Inclination of the orbit          [rad]
% OMEGA       - RAAN of the orbit                 [rad]
% R_ECI       - ECI Position Vector of the s/c    [m]


%% OUTPUT

% theta       - Argument of latitude              [rad]

%% NOTES

% R_ECI is assumed to be a column vector

%% IMPLEMENTATION
%% Compute Transformation Matrix
% This transformation matrix transforms from the ECI I-J-K coordinate frame
% to the P-Q-W coordinate frame [refer to Vallado]

% define the argument of peri-apsis (0 for circular orbits)
w = 0;

ROT_ijk2pqw(1,1) = cos(OMEGA)*cos(w)-sin(OMEGA)*sin(w)*cos(i);
ROT_ijk2pqw(1,2) = -cos(OMEGA)*sin(w)-sin(OMEGA)*cos(w)*cos(i);
ROT_ijk2pqw(1,3) = sin(OMEGA)*sin(i);
ROT_ijk2pqw(2,1) = sin(OMEGA)*cos(w)+cos(OMEGA)*sin(w)*cos(i);
ROT_ijk2pqw(2,2) = -sin(OMEGA)*sin(w)+cos(OMEGA)*cos(w)*cos(i);
ROT_ijk2pqw(2,3) = -cos(OMEGA)*sin(i);
ROT_ijk2pqw(3,1) = sin(w)*sin(i);
ROT_ijk2pqw(3,2) = cos(w)*sin(i);
ROT_ijk2pqw(3,3) = cos(i);

%% Transform from ECI to PQR coordinates

R_PQW = ROT_ijk2pqw'*R_ECI;

%% Determine Argument of Latitude

% refer to page 126 of Vallado (2007) for logic

% NOTE that R_PQW(1)/r = cos(theta) && R_PQW(2)/r = sin(theta)
cosu = R_PQW(1)/r;
sinu = R_PQW(2)/r;

% there are four conditions to consider

% first quadrant
if cosu>0 && sinu>0
    theta = acos(cosu);
    
% second quadrant
elseif cosu<0 && sinu>0
    theta = acos(cosu);
    
% third quadrant
elseif cosu<0 && sinu<0
    theta = 2*pi - acos(cosu);
    
% fourth quadrant
elseif cosu>0 && sinu<0
    theta = 2*pi - acos(cosu);
    
end