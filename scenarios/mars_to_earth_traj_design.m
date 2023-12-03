% MARS TO EARTH INTERPLANETARY TRAJECTORY DESIGN
% ========================================================================
% Explores a range of times of departure and times of flight from Mars to
% Earth given a starting departure date. Generates a porkchop plot to
% visualize opportunities within the specified times.
% Requires: Aerospace Toolbox, aeroDataPackage

%% Setup, clear data, and close windows
setup_constants;

%% CONFIGURABLE PARAMETERS
t1_mars_departure = date2JD('2032-09-01 00:00');

% Times of launch
tl_int = 1; % time interval to step through [days]
tl_max = 179+60; % maximum launch delay [days]
tl_range = (0:tl_int:tl_max) + t1_mars_departure; % absolute time in days
verify_tl = datetime(tl_range,'ConvertFrom','juliandate');

% Times of arrival
% tf_int = 1; % time interval to step through [days]
% ta_min = 178; % minimum duration of transfer time [days]
% tf_max = 15*30; % maximum duration of transfer time [days]
%ta_range = (ta_min:tf_int:tf_max) + t1_mars_departure; % [days]
ta_start = date2JD('2033-06-01 00:00');
ta_end = date2JD('2033-12-31 00:00');
ta_range = ta_start:1:ta_end;
verify_ta = datetime(ta_range,'ConvertFrom','juliandate');


%% Explore design space using Lambert solver
n_cases = length(tl_range)*length(ta_range); % matrix grid of possible launch/transfer time combinations
% Data storage
dVs = NaN(length(tl_range),length(ta_range)); % indices match that of tl_range by tf_range
v1_data = NaN(length(tl_range),length(ta_range),3);
v2_data = NaN(length(tl_range),length(ta_range),3);

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
        [v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,tof);
        v1_data(i,j,:) = v1_dep; % store dep velocity
        v2_data(i,j,:) = v2_arr; % store arrival velocity for edl analysis

        % Calculate delta-Vs
        dV1 = norm(v1_dep-v1_mars_hci); % v_inf at t1
        dVs(i,j) = dV1;
    end % end loop through times of flight
end % end loop through times of departure/launch

toc


%% Porkchop plots for interplanetary transfer

figure('Name','Porkchop plot'); hold on
[x,y] = meshgrid(tl_range,ta_range);
c3s = (dVs.^2);
c3_lvls = 0:2:20;
%contour(x,y,c3s,c3_lvls)
[c_cont,h_cont] = contourf(x,y,c3s.',c3_lvls,'ShowText',true,"FaceAlpha",0.3);
clabel(c_cont,h_cont,2:2:20);
title('Mars to Earth')
xlabel('Mars departure date (UTC)')
ylabel('Earth arrival date (UTC)')
% col = colorbar;
% col.Label.String = 'C3 (km$^2$/s$^2$)'; col.Label.Interpreter = 'latex';
% set(col, 'TickLabelInterpreter', 'latex');

% Time of flight contours
tof = y-x;
tof_mo = tof / 30; 
gray = [.5 .5 .5];
[c,h] = contour(x,y,tof_mo,2:1:15,'ShowText',true,'EdgeColor',gray); 
h.LevelList = round(h.LevelList, 1);
clabel(c,h,'Color',gray);

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

[min_energy_v_inf,min_energy_c3,r1_mars_hci,v1_mars_hci,r2_earth_hci,v2_earth_hci] = mars2earth_traj(tl_min,ta_min,mu_Sun);
