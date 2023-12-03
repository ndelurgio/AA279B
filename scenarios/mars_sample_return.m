%clear
%close all;
%% Run Mars to Earth
%mars_to_earth_traj_design;
%% CONFIGURABLE PARAMETERS
load("mins.mat")
ta = ta_min; % JD
% tl = date2JD('2033-01-01 12:00'); % JD
tl = tl_min;

dt_sec = (ta-tl)*24*3600;
[v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,dt_sec);

% Initial Time
% ti_utc = datetime(2028, 1, 1, 1, 1, 1);
t_duration = days(0.5);
ti_utc = datetime(ta,'ConvertFrom','juliandate') - t_duration;
% Initial Position
earth_soi = 0.929E9;
[~,v2_earth_hci] = planet_ephemeris_hci(ta,'Earth'); % arrival
v_inf = (v2_arr - v2_earth_hci)*10^3;
% theta = deg2rad(23.45);
% T_eci2hci = [1 0 0;
%              0 cos(theta) sin(theta);
%              0 -sin(theta) cos(theta)];
% v_inf = (T_eci2hci'*v_inf')'; 
u_inf = v_inf/norm(v_inf);
% Landing Time
tf_utc = ti_utc + t_duration;
% Landing Position: Utah test & training range
range_lat_deg = 40.489831374;
range_lon_deg = -113.635330792; 
range_alt = 0;
range_LLA = [range_lat_deg,range_lon_deg,range_alt];

range_pos_tf_j2000 = lla2eci([range_lat_deg,range_lon_deg,range_alt],datetime2vec(tf_utc));
% n = cross(v_inf,range_pos_tf_j2000)/norm(cross(v_inf,range_pos_tf_j2000));
% p = cross(n,u_inf);
% capsule_pos_ti_j2000 = -u_inf*earth_soi - p*5e8;


% if ang(range_pos_tf_j2000,v_inf) < pi/2
%     t_duration = days(1);
%     tf_utc = ti_utc + t_duration;
%     range_pos_tf_j2000 = lla2eci([range_lat_deg,range_lon_deg,range_alt],datetime2vec(tf_utc));
% end



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
%% New Hyperbola Design
theta = acos(dot(v_inf,range_pos_tf_j2000)/(norm(v_inf)*norm(range_pos_tf_j2000)));
rt = norm(range_pos_tf_j2000);
energy = norm(v_inf)^2/2;
a = -mu_earth/(2*energy);

e_max = 1 - rt/a;
e_vec = 1:0.001:e_max;
nu_vec = zeros(1,length(e_vec));
psi_vec = zeros(1,length(e_vec));

for i = 1:length(e_vec)
    e = e_vec(i);
    psi_vec(i) = acos(-1/e);
    nu_vec(i) = -acos(a*(1-e^2)/(rt*e)-1/e);
end
tot_ang = psi_vec + nu_vec + theta;
figure;
hold on;
plot(e_vec,rad2deg(nu_vec))
plot(e_vec,rad2deg(psi_vec))
plot(e_vec,rad2deg(tot_ang))
plot(e_vec,rad2deg(tot_ang-pi))
legend(["$\nu$","$\psi$","$\nu+\psi+\theta$"],"Interpreter","latex","Location","southeast")
xlabel("Eccentricity")
ylabel("Angle [deg]")

f_ang = @(e) acos(-1/e) - acos(a*(1-e^2)/(rt*e)-1/e) + theta - pi;
low = 1;
high = e_max;
while abs(high-low)/2 > 1e-9
    c = (low+high)/2;
    if sign(f_ang(c)) == sign(f_ang(low))
       low = c; 
    else
        high = c;
    end
end
e_des = (low+high)/2;

x = v_inf'/norm(v_inf);
z = cross(x,range_pos_tf_j2000')/norm(cross(x,range_pos_tf_j2000'));
y = cross(z,x);
C_hyp2eci = [x,y,z];

rp_mag = a*(1-e_des);
psi = acos(-1/e);
rp_vec = [rp_mag*cos(pi-psi); rp_mag*sin(pi-psi); 0];
vp_mag = sqrt(2*mu_earth/rp_mag-mu_earth/a);
vp_vec = cross(rp_vec/norm(rp_vec),[0;0;1])*vp_mag;
capsule_pos_tp_j2000 = (C_hyp2eci*rp_vec)';
capsule_vel_tp_j2000 = (C_hyp2eci*vp_vec)';

nuf = -acos(a*(1-e_des^2)/(rt*e_des)-1/e_des);
nui = -acos(a*(1-e_des^2)/(earth_soi*e_des)-1/e_des);
[a,e,i,Om,w,~,~,~,~] = rv2orb(capsule_pos_tp_j2000', capsule_vel_tp_j2000', mu_earth);
[capsule_pos_ti_j2000, capsule_vel_ti_j2000] = orb2rv(a*(1-e^2)/1000,e,i,Om,w,nui);
capsule_pos_ti_j2000 = capsule_pos_ti_j2000'*1000;
capsule_vel_ti_j2000 = capsule_vel_ti_j2000'*1000;

n = sqrt(mu_earth/abs(a)^3);
t_lambert = (nu2M_hyp(nuf,e_des) - nu2M_hyp(nui,e_des))/n;
t_sim = t_lambert + 60;


%% MODEL
% t_lambert = seconds(tf_utc - ti_utc);
% % [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
% % [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_vallado_u(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 's', 0, t_lambert);
% [capsule_pos_ti_j2000,t_lambert,iter] = reentry_shooter_hyperbola(capsule_pos_ti_j2000,t_lambert,v_inf,range_pos_tf_j2000,mu_earth);
% [capsule_vel_ti_j2000, capsule_vel_tf_j2000, error_out] = AA279lambert_curtis(mu_earth, capsule_pos_ti_j2000, range_pos_tf_j2000, 'retro', 0, t_lambert);
% t_sim = t_lambert + 60;
%% FODE
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
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
