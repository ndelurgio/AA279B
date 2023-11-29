% MARS TO EARTH INTERPLANETARY PATCHED CONIC
% Given initial and final absolute dates, finds a Lambert solution from
% Mars to Earth or Earth to Mars (or any 2 specified bodies)

% Requires: Aerospace Toolbox, aeroDataPackage

%% Setup
% setup_constants;

%% CONFIGURABLE PARAMETERS

start_body = 'Mars';
target_body = 'Earth';
t1_mars_departure = date2JD('2033-01-14 00:00:00'); %'2020-07-30 04:50'); %'2031-01-01 00:00'); %'2005-06-20 00:00'); %'2020-07-30 04:50');
t2_earth_arrival = date2JD('2033-08-20'); %'2021-02-18 00:00');  %'2033-08-01 00:00'); %'2021-02-18 00:00');
% Mars 2020: C3 14.49 km^2/s^2, tof 213 days (https://www.jpl.nasa.gov/news/press_kits/mars_2020/launch/mission/)

%% LAMBERT SOLVER

dt_sec = timeOfFlight(t1_mars_departure,t2_earth_arrival); % Time of transfer
[r1_mars_hci,v1_mars_hci] = planet_ephemeris_hci(t1_mars_departure,start_body); % start position
[r2_earth_hci,v2_earth_hci] = planet_ephemeris_hci(t2_earth_arrival,target_body); % end position

nrev = 0; % number of revolutions in trajectory
[v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,dt_sec);

v_inf = norm(v1_dep - v1_mars_hci);
c3 = v_inf^2;

% Propagate orbit with Lambert initial conditions from Mars
ic = [r1_mars_hci,v1_dep];

% Simulation settings
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
n_steps = 1000;
tspan = linspace(0,dt_sec,n_steps);
[t,traj] = ode113(@fode,tspan,ic,options,mu_Sun);

% Earth and Mars positions
tspan_d = linspace(t1_mars_departure,t2_earth_arrival,n_steps);
r_mars = planet_ephemeris_hci(tspan_d','Mars');
r_earth = planet_ephemeris_hci(tspan_d','Earth');

%% Trajectory visualization
earth_blue = [11,98,212]/255;
mars_red = [252,80,3]/255;
figure('Name','Trajectory'); visualize_heliocentric_3d()
scatter3(r1_mars_hci(1),r1_mars_hci(2),r1_mars_hci(3),36,mars_red,'filled','Marker','o')
scatter3(r2_earth_hci(1),r2_earth_hci(2),r2_earth_hci(3),36,earth_blue,'filled','Marker','o')
plot3(traj(:,1),traj(:,2),traj(:,3))
plot3(r_mars(:,1),r_mars(:,2),r_mars(:,3),'Color',mars_red)
plot3(r_earth(:,1),r_earth(:,2),r_earth(:,3),'Color',earth_blue)
legend('Sun','Mars','Earth','Trajectory')
title([start_body,' to ', target_body])
axis equal

tof_days = t2_earth_arrival-t1_mars_departure;
disp('Time of flight (days):'); disp(tof_days);
disp('V_inf (km/s):'); disp(v_inf);
disp('C3 (km^2/s^2:'); disp(c3);

% Display text on plot
annstr = sprintf('TOF (days): %0.3g \n $v_{inf}$ (km/s): %0.3g \n C3 (km$^2$/s$^2$): %0.3g', ...
    tof_days,v_inf,c3'); % annotation text
annpos = [0.7 0.5 0.1 0.1]; % annotation position in figure coordinates
ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
ha.HorizontalAlignment = 'left';
ha.BackgroundColor = [1 1 1];

