function [i,OMEGA,theta] = ref_orbit_classical_elements(tau)
%% DESCRIPTION:
% To compute the orbital elements of the reference circular orbit at a
% given time tau

% Based on the equations given by Schweighart & Sedwick (2002)

%% INPUT

% tau              - Dimensionless time tau                          [-] 

%% OUPUT

% i                - Inclination of the reference orbit at time tau  [rad]
% OMEGA            - RAAN of the reference orbit at time tau         [rad]
% theta            - Argument of latitude of the reference orbit     [rad]

%% IMPLEMENTATION:
%% DEFINE GLOBAL VARIABLES TO BE USED 
global mu J_2 R_e k r_ref i_ref OMEGA_ref theta_ref n_0 tau_0

%% COMPUTE ELEMENTS

% Current inclination of the reference orbit
i = i_ref-3*sqrt(mu)*J_2*R_e^2*cos(i_ref)*sin(i_ref)*(sin(k*(tau-tau_0)/n_0)^2)/(2*k*(r_ref^(7/2)));

% Current longitude of the ascending node of the reference orbit
OMEGA = OMEGA_ref-3*sqrt(mu)*J_2*R_e^2*cos(i_ref)*((tau-tau_0)/n_0)/(2*(r_ref^(7/2)));

% Current argument of latitude of the reference orbit
theta = theta_ref+k*(tau-tau_0)/n_0;