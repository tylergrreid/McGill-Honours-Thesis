function [R_ECI,V_ECI] = ref_orbit_classical_elements2ECI(r,e,i,OMEGA,w,u)
%% DESCRIPTION

% Given the classical orbital elements of an INCLINED CIRCULAR ORBIT:
%   - radius(r)
%   - eccentricity (e = 0)
%   - inlination (i)
%   - Right ascension of the ascending node (OMEGA)
%   - Argument of perigee (w)
%   - Argument of Latitude (u) - a.k.a. theta
% Ouput the velocity of the s/c in the ECI frame

% Algorithm is based on Algorithm 10 (COE2RV) of 'Fundamentals of
% Astrodynamics and Applications 3ed.' by David A. Vallado (2007)
% p.126-127 and the thesis of Mathieu Landry

%% INPUT

% R      - radius of circular reference orbit            - [m]
% e      - Eccentricity = 0                              - [n/a]
% i      - Inclination                                   - [rad]
% OMEGA  - Right ascension of the ascending node (RAAN)  - [rad]
% w      - Argument of perigee = 0                       - [rad]
% u      - Argument of latitude [a.k.a. theta]           - [rad]

%% OUTPUT

% R_ECI  - Position vector of the s/c in the ECI frame   - [m]
% V_ECI  - Velocity vector of the s/c in the ECI frame   - [m/s]

%% NOTES
 
% (1) - a,e,i,OMEGA,w,u ARE SCALAR QUANTITIES
% (2) - V & R ARE ROW VECTOR QUANTITIES
% (3) - THIS ALGORITHM ONLY GIVES THE VELOCITY IN THE ECI FRAME OF AN
%       INCLINED CIRCULAR ORBIT

%% IMPLEMENTATION:
%% DEFINE GLOBAL VARIABLES

global mu
%% DETERMINE ECI POSITION

R_ECI(1) = r*(cos(OMEGA)*cos(u)-sin(OMEGA)*cos(i)*sin(u));
R_ECI(2) = r*(sin(OMEGA)*cos(u)+cos(OMEGA)*cos(i)*sin(u));
R_ECI(3) = r*(sin(i)*sin(u));

%% DETERMINE ECI VELOCITY

const1 = sqrt(mu/r);

V_ECI(1) = const1*(-cos(OMEGA)*sin(u)-sin(OMEGA)*cos(i)*cos(u));
V_ECI(2) = const1*(-sin(OMEGA)*sin(u)+cos(OMEGA)*cos(i)*cos(u));
V_ECI(3) = const1*(sin(i)*cos(u));
