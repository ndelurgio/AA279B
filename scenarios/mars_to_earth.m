% MARS TO EARTH INTERPLANETARY PATCHED CONIC

% Requires: Aerospace Toolbox, aeroDataPackage

%% Setup
setup_constants;

%% CONFIGURABLE PARAMETERS

% Mars 2020: C3 14.49 km^2/s^2, tof 213 days (https://www.jpl.nasa.gov/news/press_kits/mars_2020/launch/mission/)

t1_mars_departure = dateToJD('2020-07-30 04:50'); %'2031-01-01 00:00'); %'2005-06-20 00:00'); %'2020-07-30 04:50');
t2_earth_arrival = dateToJD('2021-02-18 00:00');  %'2033-08-01 00:00'); %'2021-02-18 00:00');

%% LAMBERT SOLVER

dt_sec = timeOfFlight(t1_mars_departure,t2_earth_arrival); % Time of transfer
[r1_mars_hci,v1_mars_hci] = planetEphemeris(t1_mars_departure,'Sun','Mars'); % start position
[r2_earth_hci,v2_earth_hci] = planetEphemeris(t2_earth_arrival,'Sun','Earth'); % end position
% [r1_mars_hci,v1_mars_hci] = planetEphemeris(t1_mars_departure,'Sun','Earth'); % start position
% [r2_earth_hci,v2_earth_hci] = planetEphemeris(t2_earth_arrival,'Sun','Mars'); % end position

nrev = 0; % number of revolutions in trajectory
[v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,dt_sec);

v_inf = norm(v1_dep - v1_mars_hci);
c3 = v_inf^2;

% Propagate orbit with Lambert initial conditions from Mars
ic = [r1_mars_hci,v1_dep];

% Simulation settings
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
n_steps = 1000;

[t,traj] = ode113(@fode,linspace(0,dt_sec,n_steps),ic,options,mu_Sun);

%% Trajectory visualization
figure('Name','Trajectory'); visualize_heliocentric_3d()
scatter3(r1_mars_hci(1),r1_mars_hci(2),r1_mars_hci(3),36,mars_red,'filled','Marker','o')
scatter3(r2_earth_hci(1),r2_earth_hci(2),r2_earth_hci(3),36,earth_blue,'filled','Marker','o')
plot3(traj(:,1),traj(:,2),traj(:,3))
legend('Sun','Mars','Earth','Trajectory')

disp('Time of flight (days):'); disp(tf-tl);
disp('V_inf (km/s):'); disp(v_inf);
disp('C3 (km^2/s^2:'); disp(c3);

% Display text on plot
annstr = sprintf('TOF (days): %0.3g \n v_inf (km/s): %0.3g \n C3 (km^2/s^2): %0.3g', ...
    tf-tl,v_inf,c3); % annotation text
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

function tf = timeOfFlight(jd1,jd2)
% Return time of flight in seconds between 2 Julian dates, converting days to seconds 
    tf = (jd2-jd1)*24*3600;
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

