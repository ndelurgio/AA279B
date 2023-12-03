function r_eci = hci2eci(r_hci)

    theta = deg2rad(23.45); % Earth-ecliptic angle

    % Transformation matrix (clockwise rotation about x)
    T_eci2hci = [1 0 0;
                 0 cos(theta) sin(theta);
                 0 -sin(theta) cos(theta)];
    T_hci2eci = T_eci2hci';
    r_eci = (T_hci2eci*r_hci);
end