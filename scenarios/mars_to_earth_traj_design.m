% MARS TO EARTH INTERPLANETARY TRAJECTORY DESIGN
% Explores a range of times of departure and times of flight from Mars to
% Earth given a starting departure date. Generates a porkchop plot to
% visualize opportunities within the specified times.
% Requires: Aerospace Toolbox, aeroDataPackage

%% Setup
setup_constants;

%% CONFIGURABLE PARAMETERS
t1_mars_departure = dateToJD('2018-03-07 00:00'); %'2031-01-01 00:00'); %'2005-06-20 00:00'); %'2020-07-30 04:50');

% Times of launch and transfer times
tl_int = 1/3; % [days]
tl_max = 112; % [days]
tl_range = (0:tl_int:tl_max) + t1_mars_departure; % absolute time in days
verify = datetime(tl_range,'ConvertFrom','juliandate');

tf_int = 1/3; % [days]
tf_min = 178;
tf_max = 13*30; % maximum duration of transfer time [days]
tf_range = (tf_min:tf_int:tf_max) + t1_mars_departure; % [days] % 14 months
verify_tf = datetime(tf_range,'ConvertFrom','juliandate');

%% Explore design space using Lambert solver
num_cases = length(tl_range)*length(tf_range); % matrix grid of possible launch/transfer time combinations
dVs = NaN*ones(length(tl_range),length(tf_range)); % indices match that of tl_range by tf_range

for i=1:length(tl_range)
    % Launch delay [days]
    % launch_delay = tl_range(i);
    % Absolute time of launch
    % t1_depart = t1_mars_departure + launch_delay; % [days]
    t1_depart = tl_range(i);
    % Initial position of Mars at launch time
    %%[r1_mars_hci,v1_mars_hci] = planetEphemeris(t1_depart,'Sun','Mars');
    [r1_mars_hci,v1_mars_hci] = planetEphemeris(t1_depart,'Sun','Earth');

    for j=1:length(tf_range)
        % Transfer time [days]
        % dt_days = tf_range(j);

        % Absolute time of arrival
        %t2_arrival = t1_depart + dt_days; % [days]
        t2_arrival = tf_range(j);
        dt_days = t2_arrival - t1_depart;
        tof = days2sec(dt_days);
        % Final position of Earth at arrival time
        %%[r2_earth_hci,v2_earth_hci] = planetEphemeris(t2_arrival,'Sun','Earth');
        [r2_earth_hci,v2_earth_hci] = planetEphemeris(t2_arrival,'Sun','Mars');

        % Lambert solver
        nrev = 0;
        [v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,tof);

        % Calculate delta-Vs
        dV1 = norm(v1_dep-v1_mars_hci); % v_inf at t1
        dVs(i,j) = dV1;
    end

end

%% Porkchop plots for interplanetary transfer
% 
% tl_dates = cell(size(tl_range));
% tf_dates = cell(size(tf_range));
% 
% for i=1:length(tl_range)
%     tl_dates(i) = datetime(tl_range(i),'ConvertFrom','juliandate');
% end
% 
% for j=1:length(tf_range)
%     tf_dates(i) = datetime(tf_range(j),'ConvertFrom','juliandate');
% end


[x,y] = meshgrid(tl_range,tf_range);

% launch_dates = datetime(x,'ConvertFrom','juliandate','Format','yyyy-MM-dd HH:mm:ss');
% arrival_dates = datetime(y,'ConvertFrom','juliandate','Format','yyyy-MM-dd HH:mm:ss');
%% C3 plot
figure('Name','Porkchop plot')
%contour(x,y,dVs.',[0:.1:15])
c3s = (dVs.^2).';
c3_lvls = 8:1:16;
[c_cont,h_cont] = contour(x,y,c3s,c3_lvls,'ShowText',true);
clabel(c_cont,h_cont,[8:2:16]);
title('Mars to Earth')
xlabel('Mars departure date (UTC)')
ylabel('Earth arrival date (UTC)')
col = colorbar;
col.Label.String = 'C3 (km$^2$/s$^2$)'; col.Label.Interpreter = 'latex';
set(col, 'TickLabelInterpreter', 'latex');

xt = xticks;
a = tl_range(1:round(length(tl_range)/length(xt)):end);
b = datetime(a,'ConvertFrom','juliandate','Format','MMM-dd-yyy');
xticklabels(string(b));

yt = yticks;
a = tf_range(1:round(length(tf_range)/length(yt)):end);
b = datetime(a,'ConvertFrom','juliandate','Format','MMM-dd-yyy');
yticklabels(string(b));



%% LAMBERT SOLVER

min_vinf = min(min(dVs));
[i_min,j_min] = find(dVs==min_vinf);
tl = tl_range(i_min);
tf = tf_range(j_min);

dt_sec = (tf-tl)*24*3600;
[r1_mars_hci,v1_mars_hci] = planetEphemeris(tl,'Sun','Mars');
[r2_earth_hci,v2_earth_hci] = planetEphemeris(tf,'Sun','Earth'); % arrival

nrev = 0; % number of revolutions / phasing
[v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,dt_sec);

v_inf = norm(v1_dep - v1_mars_hci);
c3 = v_inf^2;

% Propagate orbit with Lambert initial conditions
ic = [r1_mars_hci,v1_dep];

% Simulation settings
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
n_steps = 1000;

[t,traj] = ode113(@fode,linspace(0,dt_sec,n_steps),ic,options,mu_Sun);

% Trajectory visualization
earth_blue = [11,98,212]/255;
mars_red = [252,80,3]/255;

figure('Name','Trajectory'); visualize_heliocentric_3d()
scatter3(r1_mars_hci(1),r1_mars_hci(2),r1_mars_hci(3),36,mars_red,'filled','Marker','o')
scatter3(r2_earth_hci(1),r2_earth_hci(2),r2_earth_hci(3),36,earth_blue,'filled','Marker','o')
plot3(traj(:,1),traj(:,2),traj(:,3))
legend('Sun','Mars','Earth','Trajectory')

disp('Time of flight (days):'); disp(tf-tl);
disp('V_inf (km/s):'); disp(v_inf);
disp('C3 (km^2/s^2:'); disp(c3);

% Display text on plot
annstr = sprintf('TOF (days): %0.3g \n $v_{inf}$ (km/s): %0.3g \n C3 (km$^2$/s$^2$): %0.3g', ...
    tf-tl,v_inf,c3'); % annotation text
annpos = [0.7 0.5 0.1 0.1]; % annotation position in figure coordinates
ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
ha.HorizontalAlignment = 'left';
ha.BackgroundColor = [1 1 1];




%% Helper functions
function jd = dateToJD(datestring)
% Inputs:
%   datestring - date in UTC
% Outputs:
%   jd - Julian date in days
    D = datetime(datestring);
    jd = juliandate(D);
end


function visualize_heliocentric_3d()
% Visualizes a body of specified radius
% Outputs: Plot of body
    
    scatter3(0,0,0,72,[227,193,2]/255,'filled','o')
    hold on; grid on; %axis equal;
    title('Heliocentric Inertial Frame')
    xlabel('X');
    ylabel('Y');
    zlabel('Z')
    view(3);
end

function s = days2sec(d)
    s = d*24*3600;
end

