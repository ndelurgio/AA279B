clear
close all;
%% Run Mars to Earth
mars_to_earth_traj_design;
%% CONFIGURABLE PARAMETERS

% Initial Time
% ti_utc = datetime(2028, 1, 1, 1, 1, 1);
t_duration = days(1);
ti_utc = datetime(tf,'ConvertFrom','juliandate') - t_duration;
% Initial Position
earth_soi = 0.929E9;
v_inf = (v2_arr - v2_earth_hci)*10^3;
u_inf = v_inf/norm(v_inf);
capsule_pos_ti_j2000 = -u_inf*earth_soi;
% capsule_pos_ti_j2000 = [5.0e+06,-2.0e+06,5.0e+06];
% Landing Time
% tf_utc = datetime(2028, 1, 1, 1, 10, 1);
tf_utc = ti_utc + t_duration;
% Landing Position: Utah test & training range
range_lat_deg = 40.489831374;
range_lon_deg = -113.635330792; 
range_alt = 0;
range_LLA = [range_lat_deg,range_lon_deg,range_alt];

%% CONSTANTS
mu_earth = 3.986004415E14;
dt_sec = 6;

% Space Weather
space_weather_data = aeroReadSpaceWeatherData("scenarios/data/SW-Last5Years.csv");

% Capsule Config
capsule = struct();
capsule.Cd = 2.2;
capsule.A = 3;
capsule.m = 200;
capsule.B = capsule.Cd * capsule.A / capsule.m;

%% MODEL
range_pos_tf_j2000 = lla2eci([range_lat_deg,range_lon_deg,range_alt],datetime2vec(tf_utc));
t_lambert = seconds(tf_utc - ti_utc);
% [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
% [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_vallado_u(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 's', 0, t_lambert);
[capsule_pos_ti_j2000,t_lambert,iter] = reentry_shooter_hyperbola(capsule_pos_ti_j2000,t_lambert,v_inf,range_pos_tf_j2000,mu_earth);
[capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
t_sim = t_lambert + 60;
% [capsule_pos_ti_j2000,t_lambert,iter] = unperturbed_shooter_dp(capsule_pos_ti_j2000,t_lambert,capsule_vel_ti_j2000,range_pos_tf_j2000,mu_earth);
% [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
%% Hyperbola Orbit Elements
% [a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = rv2orb(capsule_pos_ti_j2000', capsule_vel_ti_j2000', mu_earth);
%% FODE
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% [t_traj,capsule_traj] = ode113(...
%     @fode,...
%     0:dt_sec:t_lambert,...
%     [capsule_pos_ti_j2000,capsule_vel_ti_j2000],...
%     options,mu_earth);
% LLA_traj = eci2lla_datetime(capsule_traj(:,1:3),tf_utc - seconds(t_lambert) + seconds(t_traj));
%% Kepler
% t_traj = 0:dt_sec:(t_lambert);
% t_traj = linspace(0,t_lambert,length(t_traj));
[t_traj,capsule_traj] = fast_unperturbed_sim(t_sim,capsule_pos_ti_j2000,capsule_vel_ti_j2000,mu_earth);
% capsule_traj = zeros(length(t_traj),6);
% M0 = nu2M_hyp(wrapToPi(nu),ecc);
% M_hist = zeros(length(t_traj),1);
% nu_hist = zeros(length(t_traj),1);
% n = sqrt(mu_earth/abs(a)^3);
% warning ('off','all');
% for i = 1:length(t_traj)
%     M = M0 + n*t_traj(i);
%     M_hist(i) = M;
%     nu = wrapToPi(M2nu_hyp(M,ecc));
%     nu_hist(i) = nu;
%     [r,v] = orb2rv(a*(1-ecc^2)/1000,ecc,incl,RAAN,argp,wrapTo2Pi(nu));
%     capsule_traj(i,:) = [r;v]'*1000;
% end
% warning ('on','all');
LLA_traj = eci2lla_datetime(capsule_traj(:,1:3),tf_utc - seconds(t_lambert) + seconds(t_traj));

%% With Drag
dt_drag = 6;
[t_traj_drag,capsule_traj_drag] = fast_drag_sim(t_sim,capsule_pos_ti_j2000,capsule_vel_ti_j2000,mu_earth,capsule,tf_utc - seconds(t_lambert),space_weather_data,options);
% [t_traj_drag,capsule_traj_drag] = ode113(...
%     @fode_drag,...
%     0:dt_sec:t_lambert,...
%     [capsule_pos_ti_j2000,capsule_vel_ti_j2000],...
%     options,mu_earth,capsule,ti_utc,space_weather_data);
LLA_traj_drag = eci2lla_datetime(capsule_traj_drag(:,1:3),tf_utc - seconds(t_lambert) + seconds(t_traj_drag));

%% Reentry Shooter
% [vel_ti_shooter,t_duration_shooter,iter_shooter] = reentry_shooter_dv(capsule_vel_ti_j2000,t_lambert,capsule_pos_ti_j2000,range_pos_tf_j2000,ti_utc,tf_utc,capsule,mu_earth);
[pos_ti_shooter,t_duration_shooter,iter_shooter] = reentry_shooter_dp(capsule_pos_ti_j2000,t_lambert,capsule_vel_ti_j2000,range_pos_tf_j2000,tf_utc - seconds(t_lambert),tf_utc,capsule,mu_earth);
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
legend(["Lambert Solution (Kepler)","Kepler+Drag"])
%% Hyperbolas
[t_traj_2,capsule_traj_2] = fast_unperturbed_sim(2*t_sim,capsule_pos_ti_j2000,capsule_vel_ti_j2000,mu_earth);
capsule_traj_plot = capsule_traj(vecnorm(capsule_traj(:,1:3),2,2) < 10000e3, :);
capsule_traj_2_plot = capsule_traj_2(vecnorm(capsule_traj_2(:,1:3),2,2) < 10000e3, :);

figure;
hold on;
plot3(capsule_traj_plot(:,1)/1000,capsule_traj_plot(:,2)/1000,capsule_traj_plot(:,3)/1000,'-b',LineWidth=2)
plot3(capsule_traj_2_plot(round(end/2):end,1)/1000,capsule_traj_2_plot(round(end/2):end,2)/1000,capsule_traj_2_plot(round(end/2):end,3)/1000,'--b',LineWidth=2)
% axis equal;
grid on;
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
view(30,30)
% zlim([-1e6,1e6])