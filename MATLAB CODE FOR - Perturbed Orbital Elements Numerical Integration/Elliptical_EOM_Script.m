%% Elliptical Equations of Motion Script
%  By: Tyler Reid (2009/2010)
%% DESCRIPTION

%% SETUP AND DEFINITION OF UNIVERSAL CONSTANTS
clear;
clc;
format long;

global mu R_e omega_e J_2 

mu = 3.986005e5; %[km^3*s^-2]
R_e = 6.378137e3; %[km]
J_2 = 1.0826267e-3;
omega_e = 7.292e-5; %[rad/s]
OMEGA_SS = 1.991063853e-7; %[rad/s]

%% Define Spacecraft Physical Parameters

global B_inv_CHIEF B_inv_DEPUTY

Mass_Chief = 175; %[kg]
C_d_Chief = 2.3; % [-]
Area_Chief = 2.22; %[m^2]
B_inv_CHIEF = Area_Chief*C_d_Chief/Mass_Chief; %[m^2/kg]

Mass_Deputy = 175; %[kg]
C_d_Deputy = 2.3; % [-]
Area_Deputy = 2.22; %[m^2]
B_inv_DEPUTY = Area_Deputy*C_d_Deputy/Mass_Deputy; %[m^2/kg]

%% Define Initial Reference Orbital Elements

e_0 = 0;
a_0 = (R_e+300)/(1-e_0); % [km]
i_0 = 78*pi/180;%acos(-(2*OMEGA_SS*(a_0^(7/2))*(1-e_0^2)^2)/(3*R_e^2*J_2*sqrt(mu))); %[rad]
OMEGA_0 = 320*pi/180; % [rad]
omega_0 = 0*pi/180; % [rad]
M_0 = 0; % [rad]

%% Define Relative Initial Conditions

n_0 = sqrt(mu/a_0^3);

d = 100/1000; %[km]

s_0 = (3/8)*J_2*(R_e/a_0)^2*(1+3*cos(2*i_0));
c = sqrt(1+s_0);

% alpha = pi/2;

% Initial Conditions
x_0 = 0; %[km]
y_0 = d; %[km]
z_0 = 0; %[km]
xdot_0 = n_0*d/2; %[km/s]
ydot_0 = 0; %[km/s]
zdot_0 = n_0*d; %[km/s]


% % from: Sabol et al. (2001)
% x_0 = d*cos(alpha)/2 %[km]
% xdot_0 = -d*n_0*sin(alpha)/2 %[km/s]
% y_0 = 2*xdot_0/n_0 %[km]
% z_0 = 2*x_0 %[km]
% ydot_0 = -2*n_0*x_0 %[km/s]
% zdot_0 = 2*xdot_0 %[km/s]

% % Initial Conditions
% x_0 = d/2 % [km]
% y_0 = 0 % [km]
% z_0 = d % [km]
% xdot_0 = 0 % [km/s]
% ydot_0 = -2*x_0*n_0*sqrt(1+s_0) % [km/s]
% zdot_0 = 0 % [km/s]
% 
% x_0 = 0; %[km]
% y_0 = d; %[km]
% z_0 = -(omega_e/(n_0*c))*d*sin(i_0); %[km]
% xdot_0 = 0; %[km/s]
% ydot_0 = 0; %[km/s]
% zdot_0 = 0; %[km/s]

%% Numerical Integration Setup

% Define the Time of The Integration
T_initial = 0;
T_final = 3600*24;
tspan = [T_initial T_final];

% Define Intial Conditions
Y_0 = [a_0 e_0 i_0 OMEGA_0 omega_0 M_0 x_0 y_0 z_0 xdot_0 ydot_0 zdot_0];

% Integrate Numerically 
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_out,y_out] = ode45(@Equations_of_Motion,tspan,Y_0,options);

%% PLOT RESULTS
%% Plot Orbital Elements

t_out = t_out/3600;

figure;
subplot(2,3,1)
plot(t_out,y_out(:,1),'r')
title('Semi-Major Axis a')
xlabel('Time [hrs]')
ylabel('a [km]')

subplot(2,3,2);
plot(t_out,y_out(:,2),'r')
title('Eccentricity e')
xlabel('Time [hrs]')
ylabel('e')

subplot(2,3,3)
plot(t_out,y_out(:,3)*180/pi,'r')
title('Inclination i')
xlabel('Time [hrs]')
ylabel('i [deg]')

subplot(2,3,4)
plot(t_out,y_out(:,4)*180/pi,'r')
title('RAAN \Omega')
xlabel('Time [hrs]')
ylabel('\Omega [deg]')

subplot(2,3,5)
plot(t_out,y_out(:,5)*180/pi,'r')
title('Argument of Perigee \omega')
xlabel('Time [hrs]')
ylabel('\omega [deg]')

subplot(2,3,6)
plot(t_out,y_out(:,6)*180/pi,'r')
title('Mean Anomaly')
xlabel('Time [hrs]')
ylabel('M [deg]')

%% Plot Relative Motion

Y_x = y_out(:,7)*1000;
Y_y = y_out(:,8)*1000;
Y_z = y_out(:,9)*1000;

figure;
plot3(Y_x,Y_y,Y_z);
axis square

figure;
title('Relative Motion of the Deputy With Respect to the Chief');
subplot(2,2,1);plot3(Y_z,Y_y,Y_x,'b');
zlabel({'Radial';'[m]'});
ylabel({'In-Track';'[m]'});
xlabel({'Cross-Track';'[m]'});
axis([-100 100 -100 100 -60 60]);
set(gca,'XTick',-100:50:100)
set(gca,'YTick',-100:50:100)
set(gca,'ZTick',-60:20:60)
grid on
hold on;
subplot(2,2,2);plot(Y_x,Y_y,'b');
xlabel({'Radial';'[m]'});
ylabel({'In-Track';'[m]'});
axis([-60 60 -120 120]);
set(gca,'XTick',-60:20:60)
set(gca,'YTick',-120:20:120)
grid on;
hold on;
subplot(2,2,3);plot(Y_z,Y_x,'b');
xlabel({'Cross-Track';'[m]'});
ylabel({'Radial';'[m]'});
axis([-120 120 -60 60]);
set(gca,'XTick',-120:20:120)
set(gca,'YTick',-60:20:60)
grid on
hold on;
subplot(2,2,4);plot(Y_z,Y_y,'b');
xlabel({'Cross-Track';'[m]'});
ylabel({'In-Track';'[m]'});
axis([-120 120 -120 120]);
set(gca,'XTick',-120:20:120)
set(gca,'YTick',-120:20:120)
grid on;
hold off;
%Create textbox
% annotation(figure(1),'textbox',[0.1865 0.9244 0.6943 0.0836],...
%         'String',{'1 Day Simulation - J_2 and Aerodynamic Drag Perturbations - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - Updated Drag Equations (#0000111)'},...
%     'HorizontalAlignment','center',...
%     'LineStyle','none');

%% Plot Axes As a Function of Time

figure;
subplot(3,1,1);plot(t_out,Y_x,'b');
xlabel('Time [hrs]');
ylabel('Radial [m]');
axis tight;
hold on;
subplot(3,1,2);plot(t_out,Y_y,'r');
xlabel('Time [hrs]');
ylabel('In-Track [m]');
axis tight;
hold on;
subplot(3,1,3);plot(t_out,Y_z,'g');
xlabel('Time [hrs]');
ylabel('Cross-Track [m]');
axis tight;
hold off
% annotation(figure(2),'textbox',[0.1865 0.9244 0.6943 0.0836],...
%         'String',{'1 Day Simulation - J_2 and Aerodynamic Drag Perturbations - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - Updated Drag Equations (#0000111)'},...
%     'HorizontalAlignment','center',...
%     'LineStyle','none');