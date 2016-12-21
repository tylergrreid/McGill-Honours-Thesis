function X = COE2RV(a,e,i,OMEGA,omega,f)
%% DESCRIPTION:
% This is based on Vallado (2007) Algorithm 10

% This algorithm will only compute the Position Vector Based on the Given
% Classiscal Orbital Elements.

% All input angles are in radians [rad]
% Semi major axis may be in [km] or [m] as long as R_e is consistant

%% DEFINE NECESSARY CONSTANTS

global R_e

%% MAIN ALGORITHM

% normalize semi-major axis wrt R_e for calc purposes
a = a/R_e;

p = a*(1-e^2);

% Define the position vectpr in PQW coordinates
X(1,1) = p*cos(f)/(1+e*cos(f));
X(2,1) = p*sin(f)/(1+e*cos(f));
X(3,1) = 0;

% Define Transformation Matrix To IJK
Trans(1,1) = cos(OMEGA)*cos(omega)-sin(OMEGA)*sin(omega)*cos(i);
Trans(1,2) = -cos(OMEGA)*sin(omega)-sin(OMEGA)*cos(omega)*cos(i);
Trans(1,3) = sin(OMEGA)*sin(i);

Trans(2,1) = sin(OMEGA)*cos(omega)+cos(OMEGA)*sin(omega)*cos(i);
Trans(2,2) = -sin(OMEGA)*sin(omega)+cos(OMEGA)*cos(omega)*cos(i);
Trans(2,3) = -cos(OMEGA)*sin(i);

Trans(3,1) = sin(omega)*sin(i);
Trans(3,2) = cos(omega)*sin(i);
Trans(3,3) = cos(i);

X = Trans*X;

% Transform Back to [km]
X = X*R_e;