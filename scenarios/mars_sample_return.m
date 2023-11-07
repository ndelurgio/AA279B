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
range_LLA = [range_lat_deg,range_lon_deg,range_alt];

%% CONSTANTS
mu_earth = 3.986004415E14;
dt_sec = 1;

% Space Weather
space_weather_data = aeroReadSpaceWeatherData('scenarios/data/SW-Last5Years.csv');

% Capsule Config
capsule = struct();
capsule.Cd = 2.2;
capsule.A = 3;
capsule.m = 200;
capsule.B = capsule.Cd * capsule.A / capsule.m;

%% MODEL
range_pos_tf_j2000 = lla2eci([range_lat_deg,range_lon_deg,range_alt],datetime2vec(tf_utc));
t_lambert = seconds(tf_utc - ti_utc);
[capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'pro', 0, t_lambert);

% FODE
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_traj,capsule_traj] = ode113(...
    @fode,...
    0:dt_sec:t_lambert,...
    [capsule_pos_ti_j2000,capsule_vel_ti_j2000],...
    options,mu_earth);
LLA_traj = eci2lla_datetime(capsule_traj(:,1:3),ti_utc + seconds(t_traj));

% With Drag
[t_traj_drag,capsule_traj_drag] = ode113(...
    @fode_drag,...
    0:dt_sec:t_lambert,...
    [capsule_pos_ti_j2000,capsule_vel_ti_j2000],...
    options,mu_earth,capsule,ti_utc,space_weather_data);
LLA_traj_drag = eci2lla_datetime(capsule_traj_drag(:,1:3),ti_utc + seconds(t_traj_drag));

% Reentry Shooter
% [vel_ti_shooter,t_duration_shooter,iter_shooter] = reentry_shooter_dv(capsule_vel_ti_j2000,t_lambert,capsule_pos_ti_j2000,range_pos_tf_j2000,ti_utc,tf_utc,capsule,mu_earth);
[pos_ti_shooter,t_duration_shooter,iter_shooter] = reentry_shooter_dp(capsule_pos_ti_j2000,t_lambert,capsule_vel_ti_j2000,range_pos_tf_j2000,ti_utc,tf_utc,capsule,mu_earth);
ti_shooter = tf_utc - seconds(t_duration_shooter);
[t_traj_shooter,capsule_traj_shooter] = ode113(...
    @fode_drag,...
    0:dt_sec:t_duration_shooter,...
    ...% [capsule_pos_ti_j2000,vel_ti_shooter],...
    [pos_ti_shooter,capsule_vel_ti_j2000],...
    options,mu_earth,capsule,ti_shooter,space_weather_data);
LLA_traj_shooter = eci2lla_datetime(capsule_traj_shooter(:,1:3),ti_shooter + seconds(t_traj_shooter));
%% Plotting
uif = uifigure;
g = geoglobe(uif,'NextPlot','add');
geoplot3(g,LLA_traj(:,1),LLA_traj(:,2),LLA_traj(:,3),'LineWidth',1,'Color','red')
hold(g,'on');
geoplot3(g,LLA_traj_drag(:,1),LLA_traj_drag(:,2),LLA_traj_drag(:,3),'LineWidth',2,'Color','cyan')
geoplot3(g,LLA_traj_shooter(:,1),LLA_traj_shooter(:,2),LLA_traj_shooter(:,3),'LineWidth',2,'Color','yellow')
% legend(["Lambert Solution (Kepler)","Kepler+Drag"])