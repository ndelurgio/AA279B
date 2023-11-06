clear
close all;
%% CONFIGURABLE PARAMETERS

% Initial Time
ti_utc = datetime(2028, 1, 1, 1, 1, 1);
% Initial Position
capsule_pos_ti_j2000 = [5.0e+06,-2.0e+06,5.0e+06];
% Landing Time
tf_utc = datetime(2028, 1, 1, 1, 10, 1);
% Landing Position: Utah test & training range
range_lat_deg = 40.489831374;
range_lon_deg = -113.635330792; 
range_alt = 0;

%% CONSTANTS
mu_earth = 3.986004415E14;
dt_sec = 6;

%% MODEL
range_pos_tf_j2000 = lla2eci([range_lat_deg,range_lon_deg,range_alt],datetime2vec(tf_utc));
t_lambert = seconds(tf_utc - ti_utc);
[capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'pro', 0, t_lambert);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_traj,capsule_traj] = ode113(@fode,0:dt_sec:t_lambert,[capsule_pos_ti_j2000,capsule_vel_ti_j2000],options,mu_earth);

LLA_traj = eci2lla_datetime(capsule_traj(:,1:3),ti_utc + seconds(t_traj));

%% Plotting
uif = uifigure;
g = geoglobe(uif);
geoplot3(g,LLA_traj(:,1),LLA_traj(:,2),LLA_traj(:,3),"c")