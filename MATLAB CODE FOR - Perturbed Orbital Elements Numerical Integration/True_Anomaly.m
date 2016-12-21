function f = True_Anomaly(M,e)
%% DESCRIPTION:
% Given the mean anomaly M and the eccentricity e, output the true anomaly
% f.  This algormithm is based on Algorithm 2 and 6 from Fundamentals of
% Astrodynamics and Applications (Vallado 2007).

% This algorithm is valid only in the supposed context, i.e. closed
% elliptical orbits.

%% Find Eccentric Anomaly

if (M>-pi && M<0) || M>pi
    E0 = M-e;
else
    E0 = M+e;
end

TRUE = 1;
FALSE = 0;

DONE = FALSE;

while DONE == FALSE
    E1 = E0+(M-E0+e*sin(E0))/(1-e*cos(E0));
    
    if abs(E1-E0)<10^-10
        DONE = TRUE;
    else
        E0 = E1;
    end
end

%% Find True Anomaly

Y = sin(E1)*sqrt(1-e^2)/(1-e*cos(E1));

X = (cos(E1)-e)/(1-e*cos(E1));

f = atan2(Y,X);
