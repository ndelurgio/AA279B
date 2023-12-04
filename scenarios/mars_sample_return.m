% MARS SAMPLE RETURN
% Interplanetary approach to Earth entry, descent and landing
% ========================================================================
%% Run Mars to Earth
%clear
%close all;
%mars_to_earth_traj_design;

%% CONFIGURABLE PARAMETERS
% Flight path angle limit
fpa_limit = 25; % [deg]

% Landing Position: Utah test & training range
range_lat_deg = 40.489831374;
range_lon_deg = -113.635330792; 
range_alt = 0;
range_LLA = [range_lat_deg,range_lon_deg,range_alt];

% Scheduled time of arrival and time of launch in Julian days
load("scenarios/data/mins.mat")
% load("mins.mat")
ta = ta_min_dv; % JD
tl = tl_min_dv;
% ta = ta_min; % JD
% tl = tl_min;

%% CONSTANTS
mu_earth = 3.986004415E14;
earth_soi = 0.929E9;
dt_sec = 6;

% Space Weather
space_weather_data = aeroReadSpaceWeatherData("scenarios/data/SW-Last5Years.csv");

% Capsule Config
% Source: https://ntrs.nasa.gov/api/citations/20080023907/downloads/20080023907.pdf
capsule = struct();
capsule.hca = 60; % [deg], half-cone angle of blunt body
capsule.Cd = 2*sind(capsule.hca)^2; % drag coefficient
capsule.A = pi*(.9/2)^2; % [m^2], frontal area
capsule.m = 44; % [kg], mass
capsule.B = capsule.Cd * capsule.A / capsule.m; % ballistic coefficient

%% Compute position at landing time and arrival velocity
tof_sec = (ta-tl)*24*3600;
[r1_mars_hci,v1_mars_hci] = planet_ephemeris_hci(tl,'Mars');
[r2_earth_hci,v2_earth_hci] = planet_ephemeris_hci(ta,'Earth'); % arrival
[~,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,tof_sec);

% Initial time = time of arrival (approximate)
ti_utc = datetime(ta,'ConvertFrom','juliandate');

% Initial Position
v_inf_hci = (v2_arr - v2_earth_hci)*10^3; % [m/s]
v_inf = hci2eci(v_inf_hci')'; % excess arrival velocity in ECI expressed in Earth-centered coordinates
u_inf = v_inf/norm(v_inf);

% Landing Times - search within a 24-hour window to span range
tf_utc_range = ti_utc + hours(.1:0.1:24);

%% Loop through range of times of flight to find the best hyperbola for a desired FPA
fpa_data = NaN(length(tf_utc_range),1);
ve_data = NaN(length(tf_utc_range),1);
aer_data = NaN(length(tf_utc_range),3);
tof_range = hours(tf_utc_range-ti_utc);

for k=1:length(tf_utc_range)
    tf_utc = tf_utc_range(k);
    range_pos_tf_j2000 = lla2eci(range_LLA,datetime2vec(tf_utc));
% n = cross(v_inf,range_pos_tf_j2000)/norm(cross(v_inf,range_pos_tf_j2000));
% p = cross(n,u_inf);
% capsule_pos_ti_j2000 = -u_inf*earth_soi - p*5e8;


% if ang(range_pos_tf_j2000,v_inf) < pi/2
%     t_duration = days(1);
%     tf_utc = ti_utc + t_duration;
%     range_pos_tf_j2000 = lla2eci([range_lat_deg,range_lon_deg,range_alt],datetime2vec(tf_utc));
% end

    % New Hyperbola Design
    [a,e,i,Om,w,converge] = computeHyperbola(range_pos_tf_j2000,v_inf,mu_earth);
    if converge
        nuf = -acos(a*(1-e^2)/(norm(range_pos_tf_j2000)*e)-1/e);
        nui = -acos(a*(1-e^2)/(earth_soi*e)-1/e);
        [capsule_pos_ti_j2000, capsule_vel_ti_j2000] = orb2rv(a*(1-e^2)/1000,e,i,Om,w,nui);
        capsule_pos_ti_j2000 = capsule_pos_ti_j2000'*1000;
        capsule_vel_ti_j2000 = capsule_vel_ti_j2000'*1000;
        
        n = sqrt(mu_earth/abs(a)^3);
        t_lambert = (nu2M_hyp(nuf,e) - nu2M_hyp(nui,e))/n;
        t_sim = t_lambert + 60;
    
        % Design for flight path angle
        [r_e,v_e,fpa_e] = hyperbola_edl(a,e,i,Om,w,mu_earth);
        
        % Store data
        fpa_data(k) = fpa_e;
        ve_data(k) = norm(v_e);
        aer_data(k,:) = eci2aer(r_e',datetime2vec(tf_utc),range_LLA); % compute azimuth, elevation, range from landing site
    end
end

% Plot angle vs. tof
figure('Name','FPA vs. TOF'); hold on
yyaxis left
plot(tof_range,fpa_data)
plot(tof_range,aer_data(:,1))
ylabel('Angle (deg)')
yyaxis right
plot(tof_range,ve_data/10^3) % [km/s]
xlabel('TOF (hours)')
ylabel('Entry velocity (km/s)')
legend('FPA','Azimuth','Velocity')

%% Choose an angle at the earliest time it hits below 25 deg
idx = find(fpa_data<fpa_limit,1);
tf_utc = tf_utc_range(idx);
range_pos_tf_j2000 = lla2eci(range_LLA,datetime2vec(tf_utc));
% Design hyperbola
[a,e,i,Om,w] = computeHyperbola(range_pos_tf_j2000,v_inf,mu_earth);
nuf = -acos(a*(1-e^2)/(norm(range_pos_tf_j2000)*e)-1/e);
nui = -acos(a*(1-e^2)/(earth_soi*e)-1/e);
[capsule_pos_ti_j2000, capsule_vel_ti_j2000] = orb2rv(a*(1-e^2)/1000,e,i,Om,w,nui);
capsule_pos_ti_j2000 = capsule_pos_ti_j2000'*1000;
capsule_vel_ti_j2000 = capsule_vel_ti_j2000'*1000;

n = sqrt(mu_earth/abs(a)^3);
t_lambert = (nu2M_hyp(nuf,e) - nu2M_hyp(nui,e))/n;
t_sim = t_lambert + 60;

% Flight path angle / EDL parameters
[r_e,v_e,fpa_e] = hyperbola_edl(a,e,i,Om,w,mu_earth);
disp('Earth entry velocity (km/s):'); disp(norm(v_e)/10^3)
disp('Earth re-entry angle (deg):'); disp(fpa_e)


%% MODEL
% t_lambert = seconds(tf_utc - ti_utc);
% % [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
% % [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_vallado_u(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 's', 0, t_lambert);
% [capsule_pos_ti_j2000,t_lambert,iter] = reentry_shooter_hyperbola(capsule_pos_ti_j2000,t_lambert,v_inf,range_pos_tf_j2000,mu_earth);
% [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
% t_sim = t_lambert + 60;
%% FODE
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% [t_traj,capsule_traj] = ode113(...
%     @fode,...
%     0:dt_sec:t_lambert,...
%     [capsule_pos_ti_j2000,capsule_vel_ti_j2000],...
%     options,mu_earth);
% LLA_traj = eci2lla_datetime(capsule_traj(:,1:3),tf_utc - seconds(t_lambert) + seconds(t_traj));
%% Kepler
[t_traj,capsule_traj] = fast_unperturbed_sim(t_sim,capsule_pos_ti_j2000,capsule_vel_ti_j2000,mu_earth);
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

%% Geoglobe visualization
uif = uifigure;
g = geoglobe(uif,'NextPlot','add');
geoplot3(g,LLA_traj(:,1),LLA_traj(:,2),LLA_traj(:,3),'LineWidth',1,'Color','red')
hold(g,'on');
geoplot3(g,LLA_traj_drag(:,1),LLA_traj_drag(:,2),LLA_traj_drag(:,3),'LineWidth',2,'Color','cyan')
geoplot3(g,LLA_traj_shooter(:,1),LLA_traj_shooter(:,2),LLA_traj_shooter(:,3),'LineWidth',2,'Color','yellow')

%% Hyperbolas
[t_traj_2,capsule_traj_2] = fast_unperturbed_sim(2*t_sim,capsule_pos_ti_j2000,capsule_vel_ti_j2000,mu_earth);
capsule_traj_plot = capsule_traj(vecnorm(capsule_traj(:,1:3),2,2) < 100000e3, :);
capsule_traj_2_plot = capsule_traj_2(vecnorm(capsule_traj_2(:,1:3),2,2) < 100000e3, :);

figure;
hold on;
plot3(capsule_traj_plot(:,1)/1000,capsule_traj_plot(:,2)/1000,capsule_traj_plot(:,3)/1000,'-b',LineWidth=2)
plot3(capsule_traj_2_plot(round(end/2):end,1)/1000,capsule_traj_2_plot(round(end/2):end,2)/1000,capsule_traj_2_plot(round(end/2):end,3)/1000,'--b',LineWidth=2)
axis equal;
I = imread("earth2.jpg");
[x,y,z] = sphere;              % create a sphere 
s = surface(6.378e3*x,6.378e3*y,6.378e3*z);            % plot spherical surface

s.FaceColor = 'texturemap';    % use texture mapping
s.CData = flipud(I);                % set color data to topographic data
s.EdgeColor = 'none';          % remove edges
s.FaceLighting = 'gouraud';    % preferred lighting for curved surfaces
s.SpecularStrength = 0.4;      % change the strength of the reflected light

light('Position',[-1 0 1])     % add a light
grid on;
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
% view(30,30)
% zlim([-1e6,1e6])
% ylim([-1e6,1e6])
% xlim([-1e6,1e6])
xlim([-104148 77231])
ylim([-126572 54806])
zlim([-137721 43657])
view([-434 12])

function ang = ang(u1,u2)
ang = acos(dot(u1,u2)/(norm(u1)*norm(u2)));
end
