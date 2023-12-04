% Incoming Hyperbola Flight Path Angle
% Given the orbital elements of an incoming trajectory, computes the
% re-entry parameters position, velocity, and flight path angle
% Inputs:
%   [a,e,i,Om,w] - orbital elements of incoming hyperbola
%   mu - gravitational parameter of Earth [km^3/s^2]
% Outputs:
%   r_e - position at entry interface [m]
%   v_e - velocity at entry [m/s]
%   fpa_e - flight path angle with local horizontal at entry

function [r_e,v_e,fpa_e] = hyperbola_edl(a,e,i,Om,w,mu)

    h_atmos = 100e3; % [m]; altitude of the sensible atmosphere for EDL analysis
    R_Earth = 6378e3; % [m]; radius of Earth
    r_atmos = h_atmos + R_Earth;

    % Compute true anomaly at the entry interface
    nu_e = -acos(1/e*( (a*(1-e^2))/r_atmos -1) ); % negative to be on the correct side of the hyperbola
    p = a*(1-e^2);
    [r_e,v_e] = orb2rv(p/1000,e,i,Om,w,nu_e,0,0,0,mu/1e9);
    r_e = r_e*1000; % convert to [m]
    v_e = v_e*1000;

    % Calculate flight path angle at atmospheric interface
    % Local horizontal to the velocity vector    
    fpa_e = rad2deg( acos(norm(cross(r_e,v_e)) / (norm(r_e)*norm(v_e))) );

end
