%% NUMERICAL INTEGRATION OF PERTURBED ORBITAL ELEMENTS

%% DESCRIPTION
% This is the main configuration 

%% SETUP AND DEFINITION OF UNIVERSAL CONSTANTS
clear;
clc;
format long;

global mu R_e omega_e J_2 mass C_d Area rho_0

mu = 3.986005e5; %[km^3s^-2]
R_e = 6.378137e3; %[km]
J_2 = 1.086e-3;
omega_e = 7.292e-5; %[rad/s]
OMEGA_SS = 1.991063853e-7; %[rad/s]

%% Define Spacecraft Physical Parameters To Be Used
mass = 175; %[kg]
C_d = 2.3;
Area = 2.22; %[m^2]

%% Define Initial Reference Orbital Elements

e_0 = 0.02;
a_0 = (R_e+300)/(1-e_0) %[km]
i_0 = 97.13523*pi/180 %acos(-(2*OMEGA_SS*(a_0^(7/2))*(1-e_0^2)^2)/(3*R_e^2*J_2*sqrt(mu))); %[rad]
OMEGA_0 = 0; %[rad]
omega_0 = 0; %[rad]
M_0 = 0; %[rad]

i_ss = i_0*180/pi

% Calculate the IJK position of the s/c at t_initial
X_pos = COE2RV(a_0,e_0,i_0,OMEGA_0,omega_0,M_0);

% Calculate the density at this starting position
rho_0 = density_altitude_model(X_pos)


Y_0 = [a_0 i_0 OMEGA_0 e_0 omega_0 M_0]';


%% Numerical Integration

% Define the Time of The Integration
T_initial = 0;
T_final = 3600*60;
tspan = [T_initial T_final];

% Integrate Numerically The Actual Gauss VOP Equations
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,y] = ode45(@vop_ode,tspan,Y_0,options);


Y_0 = [6804.4 97.1404*pi/180 OMEGA_0 0.01947 omega_0 M_0]';
% Integrate Numerically The Simplified Equations
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_ave,y_ave] = ode45(@vop_ode_ave,tspan,Y_0,options);


% Smooth the actual data 
% for i = 1:6
%     y_fltd(:,i) = smooth(y(:,i),length(y));
% end

%% Plot The Results

t = t/3600;
t_ave = t_ave/3600;

figure;
subplot(2,3,1)
plot(t,y(:,1),'Color',[0.6 0.6 0.6]);hold on;
%plot(t,y_fltd(:,1),'m','LineWidth',2.5);hold on
plot(t_ave,y_ave(:,1),'k--','LineWidth',1.5)
title('Semi-Major Axis a')
xlabel('Time [hrs]')
ylabel('a [km]')


subplot(2,3,2);
plot(t,y(:,2)*180/pi,'Color',[0.6 0.6 0.6]);hold on;
%plot(t,y_fltd(:,2)*180/pi,'m','LineWidth',2.5);hold on
plot(t_ave,y_ave(:,2)*180/pi,'k--','LineWidth',1.5)
title('Inclination i')
xlabel('Time [hrs]')
ylabel('i [deg]')


subplot(2,3,3)
plot(t,y(:,3)*180/pi,'Color',[0.6 0.6 0.6]);hold on;
%plot(t,y_fltd(:,3)*180/pi,'m','LineWidth',2.5);hold on
plot(t_ave,y_ave(:,3)*180/pi,'k--','LineWidth',1.5)
title('RAAN \Omega')
xlabel('Time [hrs]')
ylabel('\Omega [deg]')
subplot(2,3,4)


plot(t,y(:,4),'Color',[0.6 0.6 0.6])
title('Eccentricity e');hold on;
%plot(t,y_fltd(:,4),'m','LineWidth',2.5);hold on
plot(t_ave,y_ave(:,4),'k--','LineWidth',1.5)
xlabel('Time [hrs]')
ylabel('e')


subplot(2,3,5)
plot(t,y(:,5)*180/pi,'Color',[0.6 0.6 0.6]);hold on;
%plot(t,y_fltd(:,5)*180/pi,'m','LineWidth',2.5);hold on
plot(t_ave,y_ave(:,5)*180/pi,'k--','LineWidth',1.5)
title('Argument of Perigee \omega')
xlabel('Time [hrs]')
ylabel('\omega [deg]')


subplot(2,3,6)
plot(t,y(:,6)*180/pi,'Color',[0.6 0.6 0.6]);hold on;
%plot(t,y_fltd(:,6)*180/pi,'m','LineWidth',2.5);hold on
plot(t_ave,y_ave(:,6)*180/pi,'k--','LineWidth',1.5)
title('Mean Anomaly')
xlabel('Time [hrs]')
ylabel('M [deg]')
hold off
% legend1 = legend('Actual Drift','Mean Secular Drift');
% set(legend1,'Position',[0.721 0.255 0.1611 0.06325]);

% input = y(:,4);
% output = smooth(input,1000001);
% 
% figure;
% plot(t,output);
% hold on;
% plot(t_ave,y_ave(:,4),'r','LineWidth',2);
% title('Eccentricity - Smoothed Actual Drift vs. Mean Simplified Secular Drift');
% xlabel('Time [hrs]');
% ylabel('e');
% legend1 = legend('Smoothed Actual Drift','Mean Secular Drift');
% set(legend1,'Position',[0.6868 0.7836 0.1611 0.06325]);