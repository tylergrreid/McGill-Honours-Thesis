%% BLOCK 1-C CODE 
%% SIMULATION INITIALIZATION SCRIPT
%% DETAILS

%   - Includes normalized equations of motion (non-dimensionalized)
%   - Valid for a circular reference orbit ONLY
%   - Input Parameters are in SI units (meters, seconds,radians)

%% SETUP

% Clear Memory
clear all;
% Clear Output Screen
clc;
% Ouput Format
format long;
% Initiate clock
tic;

%% DEFINE UNIVERSAL CONSTANTS
global mu omega_e R_e J_2 DRAG;

% Earth's Mean Equitorial Radius [m]
R_e = 6.37813649e6;

% Earth Gravitational Constant mu = G*M_earth [m^3/s^2]
mu = 3.986005e14;

% Mean Angular Velocity of the Earth [rad/s]
omega_e = 7.2921150e-5;

% J2
J_2 = 1.0826300e-3;

% Define logic variables
YES = 1;
NO = 0;

% To eliminate a perturbation set J_2 and/or DRAG = 0
J_2 = YES*J_2;
DRAG = YES;

%% DEFINE PHYSICAL SPACECRAFT PARAMETERS
global Area_Chief Area_Deputy Mass C_d Pannel_Area

% Cross-sectional area of the chief and deputy s/c [m^2]
Area_Chief = 2.22;
Area_Deputy = 2.22*1.000043;

% Assumed Additional Pannel Area of Both Spacecrafts [m^2]
Pannel_Area = 0;

% Mass of both spacecrafts [kg]
Mass = 175;

% Drag Coefficient of both spacecrafts [-]
C_d = 2.3;

%% DEFINE CIRCULAR REFERENCE ORBIT PARAMETERS
global r_ref i_ref OMEGA_ref e_ref theta_ref n n_0

% Define Reference orbit approximate altitude [m]
altitude = 250e3;

% Radius of circular reference orbit [m]
r_ref = altitude + R_e;

% Define the reference orbit period
n_0 = sqrt(mu/(r_ref^3));
n = n_0;

% Define the d(OMEGA)/dt Rate for a Sun Synchronous Orbit [rad/s]
OMEGA_sunsynch = 1.991063853e-7;

% Initial Inclination of circular reference orbit [rad]
cosi_sun_synch = -2*OMEGA_sunsynch*(r_ref/R_e)^2/(3*1.0826300e-3*n_0);

% Outout Sun Synchronous Inclination
acosd(cosi_sun_synch)

% Initial Inlination of the Reference Orbit
i_ref = 78*pi/180;%acos(cosi_sun_synch);

% Initial RAAN of circular reference orbit [rad]
OMEGA_ref = 320*pi/180;

% Initial Argument of Latitude of the reference orbit [rad[
theta_ref = 0;

% Eccentricity of the circular reference orbit [0 by definition]
e_ref = 0;

%% DEFINE CONSTANTS TO BE USED IN THE EQUATIONS OF MOTION & ICs
global s c k 

s = 3*J_2*(R_e^2)*(1+3*cos(2*i_ref))/(8*r_ref^2);
c = sqrt(1+s);
k = n*c + 3*n*J_2*(R_e^2)*(cos(i_ref)^2)/(2*(r_ref^2));

%% INITIAL CONDITIONS OF THE CHIEF AND DEPUTY S/C's FOR A CIRCULAR FORMATION
global d alpha

% Define the Characteristic length parameter associated with relative
% motion geometry - this could be the diameter of a projected circular 
% formation or the constant offset value in an in-track formation [m]
d = 100;

% Define alpha (needed for projected circular formation only) [rad]
alpha = 0;

% Define ICs for the chief s/c
x_chief = 0
y_chief = 0
z_chief = 0
x_dot_chief = 0
y_dot_chief = 3*J_2*R_e^2*n^2*sin(i_ref)^2/(4*k*r_ref)
z_dot_chief = 0
 
% % Define ICs for the deputy s/c
x_deputy = (d/2)*cos(alpha)
y_deputy = d*sin(alpha)
z_deputy = d*cos(alpha)
x_dot_deputy = -(y_deputy*n*(1-s)/sqrt(1+s))/2
y_dot_deputy = -2*n*x_deputy*sqrt(1+s)+3*J_2*R_e^2*n^2*sin(i_ref)^2/(4*k*r_ref)
z_dot_deputy = 2*(y_deputy*n*(1-s)/sqrt(1+s))/2

% x_deputy = 0
% y_deputy = d
% z_deputy = d
% x_dot_deputy = n_0*d/2%+(y_deputy*n*(1-s)/sqrt(1+s))/2
% y_dot_deputy = 0%-2*n*x_deputy*sqrt(1+s)+3*J_2*R_e^2*n^2*sin(i_ref)^2/(4*k*r_ref)
% z_dot_deputy = n_0*d

% x_deputy = 0
% y_deputy = d
% z_deputy = -(omega_e/(n*c))*d*sin(i_ref)
% x_dot_deputy = 0
% y_dot_deputy = 3*J_2*R_e^2*n^2*sin(i_ref)^2/(4*k*r_ref)
% z_dot_deputy = 0

%% DEFINE THE REMAINING CONSTANTS BASED ON NON-NORMALIZED ICs
global q l phi

% select a value for phi [rad]
phi = 0;

i_sat2 = i_ref;
i_sat1 = i_sat2 + (z_dot_deputy - z_dot_chief)/(k*r_ref);
delta_OMEGA_0 = (z_deputy - z_chief)/(r_ref*sin(i_ref));
gamma_0 = acot((cot(i_sat2)*sin(i_sat1)-...
    cos(i_sat1)*cos(delta_OMEGA_0))/sin(delta_OMEGA_0));
PHI_0 = acos(cos(i_sat1)*cos(i_sat2) + sin(i_sat1)*sin(i_sat2)*cos(delta_OMEGA_0));
OMEGA_sat1 = -3*n*J_2*(R_e^2)*cos(i_sat1)/(2*(r_ref^2));
OMEGA_sat2 = -3*n*J_2*(R_e^2)*cos(i_sat2)/(2*(r_ref^2));
q = n*c - (cos(gamma_0)*sin(gamma_0)*cot(delta_OMEGA_0)-...
    (sin(gamma_0)^2)*cos(i_sat1))*(OMEGA_sat1-OMEGA_sat2) - OMEGA_sat1*cos(i_sat1);
l = -r_ref*sin(i_sat1)*sin(i_sat2)*sin(delta_OMEGA_0)*(OMEGA_sat1-OMEGA_sat2)/sin(PHI_0);

%% DEFINE THE TIME BOUNDS OF THE INTEGRATION [seconds]

T_0 = 0;
T_FINAL = 3600*24;

%% NORMALIZE THE INITIAL CONDITIONS AND TIME PARAMETERS

% refer to non-dimensionalized document xh = x/d, tau = n_0*t 
% where h denotes hat

% non-dimensionalize time
tau_initial = n_0*T_0;
tau_final = n_0*T_FINAL;

% initial time offset
global tau_0
tau_0 = tau_initial;

% non-dimensionalze the chief/deputy IC's
xh_chief = x_chief/d;
yh_chief = y_chief/d;
zh_chief = y_chief/d;
xh_prime_chief = x_dot_chief/(d*n_0);
yh_prime_chief = y_dot_chief/(d*n_0);
zh_prime_chief = z_dot_chief/(d*n_0);

xh_deputy = x_deputy/d;
yh_deputy = y_deputy/d;
zh_deputy = z_deputy/d;
xh_prime_deputy = x_dot_deputy/(d*n_0);
yh_prime_deputy = y_dot_deputy/(d*n_0);
zh_prime_deputy = z_dot_deputy/(d*n_0);

%% DEFINE REFERENCE ORBIT RESET TOLERANCE
global ref_orbit_tol

ref_orbit_tol = 1e-4;

%% SOLVE THE DIFFERENTIAL EQUATIONS OF MOTION NUMERICALLY

% Initial Conditions Vector
Y0 = [xh_chief xh_prime_chief yh_chief yh_prime_chief zh_chief zh_prime_chief ...
    xh_deputy xh_prime_deputy yh_deputy yh_prime_deputy zh_deputy zh_prime_deputy ];

% Set up time bounds
tspan = [tau_initial tau_final];

% Define solver integration tolerances and options
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Function
odefun_1 = @Equations_of_Motion;

% Solve Differential Equations Numericallly
[TAU_OUT,YOUT] = ode45(odefun_1,tspan,Y0,options);


fprintf('Solution Obtained - Plotting Results\n');

%% PLOT THE RESULTS

% convert back to original units and plot
Y_x = (YOUT(:,7)-YOUT(:,1))*d;
Y_y = (YOUT(:,9)-YOUT(:,3))*d;
Y_z = (YOUT(:,11)-YOUT(:,5))*d;

% convert time back to original units
TOUT = TAU_OUT/n_0;
TOUT = TOUT/3600;
% figure;
% plot3(Y_x,Y_y,Y_z);
% xlabel({'In-Track';'[m]'});
% ylabel({'Radial';'[m]'});
% zlabel({'Cross-Track';'[m]'});
% grid on;
% axis square;

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
% 
% 
figure;
subplot(3,1,1);plot(TOUT,Y_x,'b');
xlabel('Time [hrs]');
ylabel('Radial [m]');
axis tight;
hold on;
subplot(3,1,2);plot(TOUT,Y_y,'r');
xlabel('Time [hrs]');
ylabel('In-Track [m]');
axis tight;
hold on;
subplot(3,1,3);plot(TOUT,Y_z,'g');
xlabel('Time [hrs]');
ylabel('Cross-Track [m]');
axis tight;
hold off
% % annotation(figure(2),'textbox',[0.1865 0.9244 0.6943 0.0836],...
% %         'String',{'1 Day Simulation - J_2 and Aerodynamic Drag Perturbations - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - Updated Drag Equations (#0000111)'},...
% %     'HorizontalAlignment','center',...
% %     'LineStyle','none');
% 
% figure;
% subplot(3,1,1);plot(TOUT/3600,YOUT(:,1)*d,'b');
% hold on;
% plot(TOUT/3600,YOUT(:,7)*d,'--r');
% hold off;
% xlabel('Time [hrs]');
% ylabel('Radial [m]');
% axis tight;
% subplot(3,1,2);plot(TOUT/3600,YOUT(:,3)*d,'b');
% hold on;
% plot(TOUT/3600,YOUT(:,9)*d,'--r');
% hold off;
% xlabel('Time [hrs]');
% ylabel('In-Track [m]');
% axis tight;
% subplot(3,1,3);plot(TOUT/3600,YOUT(:,5)*d,'b');
% hold on;
% plot(TOUT/3600,YOUT(:,11)*d,'--r');
% xlabel('Time [hrs]');
% ylabel('Cross-Track [m]');
% axis tight;
% legend('Chief','Deputy');
% set(legend,'Position',[0.806 0.02157 0.09668 0.06325]);
% hold off
% annotation(figure(3),'textbox',[0.1865 0.9244 0.6943 0.0836],...
%         'String',{'1 Day Simulation - J_2 and Aerodynamic Drag Perturbations - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - Updated Drag Equations (#0000111)'},...
%     'HorizontalAlignment','center',...
%     'LineStyle','none');
% 
% figure;
% subplot(3,1,1);plot(TOUT/3600,YOUT(:,2)*n_0*d,'k');
% hold on;
% plot(TOUT/3600,YOUT(:,8)*n_0*d,'--g');
% hold off;
% xlabel('Time [hrs]');
% ylabel('Radial Vel. [m/s]');
% axis tight;
% subplot(3,1,2);plot(TOUT/3600,YOUT(:,4)*n_0*d,'k');
% hold on;
% plot(TOUT/3600,YOUT(:,10)*n_0*d,'--g');
% hold off;
% xlabel('Time [hrs]');
% ylabel('In-Track Vel. [m/s]');
% axis tight;
% subplot(3,1,3);plot(TOUT/3600,YOUT(:,6)*n_0*d,'k');
% hold on;
% plot(TOUT/3600,YOUT(:,12)*n_0*d,'--g');
% xlabel('Time [hrs]');
% ylabel('Cross-Track Vel. [m/s]');
% axis tight;
% legend('Chief','Deputy');
% set(legend,'Position',[0.806 0.02157 0.09668 0.06325]);
% hold off
% annotation(figure(4),'textbox',[0.1865 0.9244 0.6943 0.0836],...
%         'String',{'1 Day Simulation - J_2 and Aerodynamic Drag Perturbations - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - Updated Drag Equations (#0000111)'},...
%     'HorizontalAlignment','center',...
%     'LineStyle','none');
% 
% % UNCOMMENT FOR A COMET PLOT
% % figure;
% % plot(0,0);
% % xlabel({'Radial';'[m]'});
% % ylabel({'In-Track';'[m]'});
% % axis([-60 60 -120 120]);
% % set(gca,'XTick',-60:20:60)
% % set(gca,'YTick',-120:20:120)
% % grid on;
% % hold on;
% % comet(Y_x,Y_y);
% 
% 
% % UNCOMMENT FOR A PLOT OF NON-DIMENSIONALIZED PARAMETERS
% figure;
% subplot(3,1,1);plot(TAU_OUT,YOUT(:,1),'b');
% hold on;
% plot(TAU_OUT,YOUT(:,7),'--r');
% hold off;
% xlabel('Time [-]');
% ylabel('Radial [-]');
% axis tight;
% subplot(3,1,2);plot(TAU_OUT,YOUT(:,3),'b');
% hold on;
% plot(TAU_OUT,YOUT(:,9),'--r');
% hold off;
% xlabel('Time [-]');
% ylabel('In-Track [-]');
% axis tight;
% subplot(3,1,3);plot(TAU_OUT,YOUT(:,5),'b');
% hold on;
% plot(TAU_OUT,YOUT(:,11),'--r');
% xlabel('Time [-]');
% ylabel('Cross-Track [-]');
% axis tight;
% legend('Chief','Deputy');
% set(legend,'Position',[0.806 0.02157 0.09668 0.06325]);
% hold off
% % annotation(figure(5),'textbox',[0.1865 0.9244 0.6943 0.0836],...
% %         'String',{'1 Day Simulation - Aerodynamic Drag Perturbations Only - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - ICs Set 2(#000090)'},...
% %     'HorizontalAlignment','center',...
% %     'LineStyle','none');
% % 
% figure;
% subplot(3,1,1);plot(TAU_OUT,YOUT(:,2),'k');
% hold on;
% plot(TAU_OUT,YOUT(:,8),'--g');
% hold off;
% xlabel('Time [-]');
% ylabel('Radial Vel. [-]');
% axis tight;
% subplot(3,1,2);plot(TAU_OUT,YOUT(:,4),'k');
% hold on;
% plot(TAU_OUT,YOUT(:,10),'--g');
% hold off;
% xlabel('Time [-]');
% ylabel('In-Track Vel. [-]');
% axis tight;
% subplot(3,1,3);plot(TAU_OUT,YOUT(:,6),'k');
% hold on;
% plot(TAU_OUT,YOUT(:,12),'--g');
% xlabel('Time [-]');
% ylabel('Cross-Track Vel. [-]');
% axis tight;
% legend('Chief','Deputy');
% set(legend,'Position',[0.806 0.02157 0.09668 0.06325]);
% hold off
% % annotation(figure(6),'textbox',[0.1865 0.9244 0.6943 0.0836],...
% %         'String',{'1 Day Simulation - Aerodynamic Drag Perturbations Only - Reference Orbit Inclination = 78.1 deg, Altitude = 250 km','TEST 01 Physical Parameters - ICs Set 2(#000090)'},...
% %     'HorizontalAlignment','center',...
% %     'LineStyle','none');

% Output runtime
toc
