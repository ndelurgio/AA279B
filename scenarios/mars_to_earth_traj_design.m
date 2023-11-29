% MARS TO EARTH INTERPLANETARY TRAJECTORY DESIGN
% ========================================================================
% Explores a range of times of departure and times of flight from Mars to
% Earth given a starting departure date. Generates a porkchop plot to
% visualize opportunities within the specified times.
% Requires: Aerospace Toolbox, aeroDataPackage

%% Setup, clear data, and close windows
setup_constants;

%% CONFIGURABLE PARAMETERS
t1_mars_departure = date2JD('2032-11-01 00:00'); %'2031-01-01 00:00'); %'2005-06-20 00:00'); %'2020-07-30 04:50');

% Times of launch
tl_int = 1/3; % time interval to step through [days]
tl_max = 179; % maximum launch delay [days]
tl_range = (0:tl_int:tl_max) + t1_mars_departure; % absolute time in days
verify_tl = datetime(tl_range,'ConvertFrom','juliandate');

% Times of arrival
tf_int = 1/3; % time interval to step through [days]
ta_min = 178; % minimum duration of transfer time [days]
tf_max = 15*30; % maximum duration of transfer time [days]
ta_range = (ta_min:tf_int:tf_max) + t1_mars_departure; % [days]
verify_ta = datetime(ta_range,'ConvertFrom','juliandate');


%% Explore design space using Lambert solver
num_cases = length(tl_range)*length(ta_range); % matrix grid of possible launch/transfer time combinations
dVs = NaN*ones(length(tl_range),length(ta_range)); % indices match that of tl_range by tf_range
% Pre-compute planet ephermides
[r1_mars_hc,v1_mars_hc] = planet_ephemeris_hci(tl_range','Mars');
[r2_earth_hc,v2_earth_hc] = planet_ephemeris_hci(ta_range','Earth');

tic
for i=1:length(tl_range)
    % Absolute time of launch
    t1_depart = tl_range(i);

    % Initial position of Mars at launch time
    r1_mars_hci = r1_mars_hc(i,:);
    v1_mars_hci = v1_mars_hc(i,:);

    for j=1:length(ta_range)
        % Absolute time of arrival
        t2_arrival = ta_range(j);
        dt_days = t2_arrival - t1_depart;
        tof = seconds(days(dt_days));
        % Final position of Earth at arrival time
        r2_earth_hci = r2_earth_hc(j,:);

        % Lambert solver
        nrev = 0;
        [v1_dep,~] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,tof);
        % Calculate delta-Vs
        dV1 = norm(v1_dep-v1_mars_hci); % v_inf at t1
        dVs(i,j) = dV1;
    end % end loop through times of flight
end % end loop through times of departure/launch

toc


%% Porkchop plots for interplanetary transfer

figure('Name','Porkchop plot')
[x,y] = meshgrid(tl_range,ta_range);
c3s = (dVs.^2).';
c3_lvls = 0:2:20;
%contour(x,y,c3s,c3_lvls)
[c_cont,h_cont] = contourf(x,y,c3s,c3_lvls,'ShowText',true,"FaceAlpha",0.3);
clabel(c_cont,h_cont,2:2:20);
title('Mars to Earth')
xlabel('Mars departure date (UTC)')
ylabel('Earth arrival date (UTC)')
% col = colorbar;
% col.Label.String = 'C3 (km$^2$/s$^2$)'; col.Label.Interpreter = 'latex';
% set(col, 'TickLabelInterpreter', 'latex');

% Set tick labels as dates
xt = xticks;
xdt = datetime(xt,'ConvertFrom','juliandate','Format','MMM-dd-yyy');
xticklabels(string(xdt));

yt = yticks;
ydt = datetime(yt,'ConvertFrom','juliandate','Format','MMM-dd-yyy');
yticklabels(string(ydt));

grid on

%% Plot 3D trajectory for minimum v_inf path

min_vinf = min(min(dVs));
[i_min,j_min] = find(dVs==min_vinf);
tl_min = tl_range(i_min);
ta_min = ta_range(j_min);

dt_sec = (ta_min-tl_min)*24*3600;
[r1_mars_hci,v1_mars_hci] = planet_ephemeris_hci(tl_min,'Mars');
[r2_earth_hci,v2_earth_hci] = planet_ephemeris_hci(ta_min,'Earth'); % arrival

nrev = 0; % number of revolutions / phasing
[v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,dt_sec);

v_inf = norm(v1_dep - v1_mars_hci);
c3 = v_inf^2;

% Propagate orbit with Lambert initial conditions
ic = [r1_mars_hci,v1_dep];

% Simulation settings
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
n_steps = 1000;
tspan = linspace(0,dt_sec,n_steps);
[t,traj] = ode113(@fode,tspan,ic,options,mu_Sun);

% Earth and Mars positions
tspan_d = linspace(tl_min,ta_min,n_steps);
r_mars = planet_ephemeris_hci(tspan_d','Mars');
r_earth = planet_ephemeris_hci(tspan_d','Earth');

% Trajectory visualization
earth_blue = [11,98,212]/255;
mars_red = [252,80,3]/255;

figure('Name','Trajectory'); visualize_heliocentric_3d()
scatter3(r1_mars_hci(1),r1_mars_hci(2),r1_mars_hci(3),36,mars_red,'filled','Marker','o')
scatter3(r2_earth_hci(1),r2_earth_hci(2),r2_earth_hci(3),36,earth_blue,'filled','Marker','o')
plot3(traj(:,1),traj(:,2),traj(:,3))
plot3(r_mars(:,1),r_mars(:,2),r_mars(:,3),'Color',mars_red)
plot3(r_earth(:,1),r_earth(:,2),r_earth(:,3),'Color',earth_blue)
legend('Sun','Mars','Earth','Trajectory')
title('Minimum energy trajectory')
axis equal

% Date of min. energy departure and arrival
dt_dep = datetime(tl_min,'ConvertFrom','juliandate');
dt_arr = datetime(ta_min,'ConvertFrom','juliandate');

disp('Departure date:'); disp(dt_dep)
disp('Arrival date:'); disp(dt_arr)
disp('Time of flight (days):'); disp(ta_min-tl_min);
disp('V_inf (km/s):'); disp(v_inf);
disp('C3 (km^2/s^2:'); disp(c3);

% Display text on plot
% annstr = sprintf('TOF (days): %0.3g \n $v_{inf}$ (km/s): %0.3g \n C3 (km$^2$/s$^2$): %0.3g', ...
%     ta_min-tl_min,v_inf,c3'); % annotation text
% annpos = [0.7 0.5 0.1 0.1]; % annotation position in figure coordinates
% ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
% ha.HorizontalAlignment = 'left';
% ha.BackgroundColor = [1 1 1];


