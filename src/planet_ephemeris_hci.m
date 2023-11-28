function [r,v] = planet_ephemeris_hci(jd,planet)
% Computes the ephemerides of a given planet at a specified date in
% heliocentric coordinates
% Inputs:
%   jd - date in Julian days
%   planet - string name of planet (e.g. 'Mars')
% Outputs:
%   r - position in HCI
%   v - velocity in HCI
    [r_eq,v_eq] = planetEphemeris(jd,'Sun',planet);
    
    % Transform from equatorial coordinates to HCI (ecliptic)
    theta = deg2rad(23.45); % Earth-ecliptic tilt angle
    
    % Transformation matrix (clockwise rotation about x)
    T_eci2hci = [1 0 0;
                 0 cos(theta) sin(theta);
                 0 -sin(theta) cos(theta)];
    r = (T_eci2hci*r_eq');
    v = (T_eci2hci*v_eq');
end