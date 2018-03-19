function [satPositions, satClkCorr] = satpos(transmitTime, prnList, ...
    eph, settings)
%SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
%given ephemeris EPH. Coordinates are calculated for each satellite in the
%list PRNLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph, settings);
%
%   Inputs:
%       transmitTime  - transmission time
%       prnList       - list of PRN-s to be processed
%       eph           - ephemeridies of satellites
%       settings      - receiver settings
%
%   Outputs:
%       satPositions  - positions of satellites (in ECEF system [X; Y; Z;])
%       satClkCorr    - correction of satellites clocks

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================
numOfSatellites = size(prnList, 2);

% GPS constatns

bdsPi          = 3.1415926535898;  % Pi used in the GPS coordinate
% system

%--- Constants for satellite position calculation -------------------------
Omegae_e     = 7.2921151467e-5;     % Earth rotation rate, [rad/s]
GM             = 3.986004418e14;   % Earth's universal
% gravitational parameter,
% [m^3/s^2]
F              = -4.44280731e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    
    %% Find initial satellite clock correction --------------------------------
    
    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime - eph(prn).t_oc);
    
    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph(prn).a2 * dt + eph(prn).a1) * dt + eph(prn).a0 - eph(prn).T_GD;
    
    time = transmitTime - satClkCorr(satNr);
    
    %% Find satellite's position ----------------------------------------------
    
    %Restore semi-major axis
    a   = eph(prn).sqrtA * eph(prn).sqrtA;
    
    %Time correction
    tk  = check_t(time - eph(prn).t_oe);
    
    %Initial mean motion
    n0  = sqrt(GM / a^3);
    %Mean motion
    n   = n0 + eph(prn).deltan;
    
    %Mean anomaly
    M   = eph(prn).M_0 + n * tk;
    %Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*bdsPi, 2*bdsPi);
    
    %Initial guess of eccentric anomaly
    E   = M;
    
    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*bdsPi);
        
        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end
    
    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*bdsPi, 2*bdsPi);
    
    %Compute relativistic correction term
    dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);
    
    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);
    
    %Compute angle phi
    phi = nu + eph(prn).omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*bdsPi);
    
    %Correct argument of latitude
    u = phi + eph(prn).C_uc * cos(2*phi) + eph(prn).C_us * sin(2*phi);
    %Correct radius
    r = a * (1 - eph(prn).e*cos(E)) + eph(prn).C_rc * cos(2*phi) + eph(prn).C_rs * sin(2*phi);
    %Correct inclination
    i = eph(prn).i_0 + eph(prn).IDOT * tk + eph(prn).C_ic * cos(2*phi) + eph(prn).C_is * sin(2*phi);
    if prn > 5 && prn ~= 17
        %Compute the angle between the ascending node and the Greenwich meridian
        Omega = eph(prn).omega_0 + (eph(prn).omega_dot - Omegae_e)*tk - Omegae_e * eph(prn).t_oe;
        %Reduce to between 0 and 360 deg
        Omega = rem(Omega + 2*bdsPi, 2*bdsPi);
        
        %--- Compute satellite coordinates ------------------------------------
        satPositions(1, satNr) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
        satPositions(2, satNr) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
        satPositions(3, satNr) = sin(u)*r * sin(i);
    else
        Omega = eph(prn).omega_0 + eph(prn).omega_dot * tk - Omegae_e * eph(prn).t_oe;
        Omega = rem(Omega + 2*bdsPi, 2*bdsPi);
        Xk = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
        Yk = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
        Zk = sin(u)*r * sin(i);
        Rx = [1,      0,       0;
            0, cosd(-5),sind(-5);
            0,-sind(-5),cosd(-5)];
        Rz = [cos(Omegae_e * tk),sin(Omegae_e * tk),0;
            -sin(Omegae_e * tk) ,cos(Omegae_e * tk),0;
              0                   ,                   0,1];
        satPositions(:, satNr) = Rz * Rx *[Xk;Yk;Zk];
    end
    
    %% Include relativistic correction in clock correction --------------------
    satClkCorr(satNr) = (eph(prn).a2 * dt + eph(prn).a1) * dt + eph(prn).a0 - eph(prn).T_GD + dtr;
    
end % for satNr = 1 : numOfSatellites
