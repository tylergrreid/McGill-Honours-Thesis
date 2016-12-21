function f_drag = relative_drag_perturbation(X_ECI,V_rel,Area_Spacecraft)
%% DESCRIPTION:
% To compute the relative atmospheric drag effects on the s/c under
% consideration.

%% INPUT

% X_ECI            - The position of the s/c in the ECI frame.    [m]
% V_rel            - The velocity of the s/c in the Hill frame    [-]
%                    w.r.t. to the rotating atmosphere.
% Area_Spacecraft  - s/c projected area normal to the direction   [m^2]  
%                    of flight.

%% OUPUT

% f_drag           - Dimensionless drag perturbation effects on the    [-]
%                    s/c in the Hill frame (vector expressed 
%                    in the hill frame).

%% IMPLEMENTATION:
%% DEFINE GLOBAL VARIABLES TO BE USED 
global Mass C_d d r_ref

%% ATMOSPHERIC DRAG EFFECTS ON  S/C

% Compute the density [kg/m^3]
rho = density_altitude_model(X_ECI);

% Compute the normalized ballistic coefficient
beta = (rho*C_d*Area_Spacecraft*r_ref/Mass)^-1;

% Compute the f_drag perturbations [m/s^2]
f_drag = -(1/2)*(1/beta)*(d/r_ref)*norm(V_rel)*V_rel;

