clear;
clc;

global R_e
R_e =  6.378137e3;
e = 0;
i = 45*pi/180;

for i = 1:500
   rp(i) = R_e+250+i;
   a(i) = rp(i)/(1-e);
   X = [rp(i) 0 0]';
   rho = density_altitude_model(X);
   B_inv(i) =  rho*a(i)*1000;
end

figure;
plot(rp-R_e,B_inv,'r','LineWidth',2);
title('\beta_{C}^{-1} as a Function of Orbit Altitude');
ylabel('\beta_{C}^{-1}');
xlabel('Altitude [km]');
grid on

