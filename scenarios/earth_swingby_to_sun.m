%% EARTH TO SUN SWING-BY
% Given the v_inf at arrival at Earth, find the hyperbola and required delta-V 
% that will put the orbiter on a swing-by path to orbit around the Sun

%% Configurable parameters
% Time of arrival at Sun
t1_mars_dep = tl_min;
t2_earth_dep = ta_min; % assumes that arrival at the planet is approximately the same date as departure
t3_sun_arr = t2_earth_dep + 8; % todo: vary this tof to sweep through possible v_inf_outbounds

%% Compute planetary ephermerides
[r1_mars,v1_mars] = planet_ephemeris_hci(t1_mars_dep,'Mars');
[r2_earth,v2_earth] = planet_ephemeris_hci(t2_earth_dep,'Earth');
%[r3_sun,v3_sun] = planet_ephemeris_hci(t3_sun_arr,'Sun');
%r3_sun = [1,1,1];

% Define some desired orbit around the Sun
a_sun = norm(r2_earth);
e_sun = 0.3;
i_sun = 0;
raan_sun = 0;
aop_sun = 0;
nu_sun = deg2rad(0);

p = a_sun*(1-e_sun^2);

[r3_sun,v3_sun] = orb2rv(p,e_sun,i_sun,raan_sun,aop_sun,nu_sun,0,0,0,mu_Sun);

%% Mars to Earth and Earth to Sun Lambert solvers for v_inf
nrev = 0;
dt12 = seconds(days(t2_earth_dep - t1_mars_dep));
[v1_mars_dep_hci,v2_earth_arr_hci] = AA279lambert_curtis(mu_Sun,r1_mars,r2_earth,'pro',nrev,dt12);

dt23 = seconds(days(t3_sun_arr - t2_earth_dep));
[v2_earth_dep_hci,v3_sun_arr_hci] = AA279lambert_curtis(mu_Sun,r2_earth,r3_sun','pro',1,dt23);
%[v2_earth_dep_hci,v3_sun_arr_hci] = AA279lambert_vallado_u(mu_Sun,r2_earth,r3_sun,'s',0,dt23); 

% Verification
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
n_steps = 100;
tspan12 = linspace(0,seconds(days(t2_earth_dep-t1_mars_dep)),n_steps);

% Capsule heliocentric trajectory (no delta-V applied)
[t1,traj12] = ode113(@fode,tspan12,[r1_mars,v1_mars_dep_hci],options,mu_Sun);
tspan23 = linspace(0,seconds(days(t3_sun_arr-t2_earth_dep)),n_steps);
[t2,traj23] = ode113(@fode,tspan23,[r2_earth,v2_earth_dep_hci],options,mu_Sun);

% Simulate final orbit around Sun
T_sun_orbit = 2*pi*sqrt(a_sun^3 / mu_Sun);
tspan_sun = linspace(0,T_sun_orbit,100);
[t,sun_orbit] = ode113(@fode,tspan_sun,[r3_sun;v3_sun],options,mu_Sun);

% figure(); hold on %visualize_heliocentric_3d() 
sun_yellow = [227,193,2]/255;
% plot3(traj12(:,1),traj12(:,2),traj12(:,3));
% plot3(traj23(:,1),traj23(:,2),traj23(:,3));
% scatter3(r1_mars(1),r1_mars(2),r1_mars(3),36,mars_red,'filled','Marker','o')
% scatter3(r2_earth(1),r2_earth(2),r2_earth(3),36,earth_blue,'filled','Marker','o')
% scatter3(0,0,0,72,sun_yellow,'filled','Marker','o')
% scatter3(r3_sun(1),r3_sun(2),r3_sun(3),20,'g')
% plot3(sun_orbit(:,1),sun_orbit(:,2),sun_orbit(:,3));
% %scatter3(traj23(end,1),traj23(end,2),traj23(end,3))
% legend('Mars to Earth','Earth to Sun','Mars','Earth','Sun','Target heliocentric orbit')
% rotate3d on
% view(3)
% axis equal
% grid on


%% Convert velocity to planetocentric
v_inf_arr_earth = v2_earth_arr_hci - v2_earth; % v-infinity at Earth arrival
v_inf_dep_earth = v2_earth_dep_hci - v2_earth; % v-infinity at Earth departure

%angleBetween(v_inf_arr_earth,v_inf_dep_earth)

% Want desired offset distance at periapsis (closest approach)
rp_des = R_Earth*3;

% Initial guess for the v_inf_arr vector
v_inf_arr = v_inf_arr_earth;

%v_inf_arr = vel_b - v_earth_dv; % take velocity at time of delta-V
%v_inf_dep_earth = norm(v_inf_arr) * v_inf_dep_earth/norm(v_inf_dep_earth);

% Shoot until rp ~ rp_des
tol = 100; % [km]
max_iter = 10;
iter = 0;
dv = 0.0001; % [km/s]
dr = Inf;
while dr > tol && iter < max_iter 
    J = zeros(3,1); % Jacobian matrix
    rp = radiusPeriapsis(v_inf_arr,v_inf_dep_earth,mu_Earth);
    dr = norm(rp_des - rp);
    % Perturb vx
    rp_xhi = radiusPeriapsis(v_inf_arr+[dv,0,0],v_inf_dep_earth,mu_Earth);
    rp_xlo = radiusPeriapsis(v_inf_arr-[dv,0,0],v_inf_dep_earth,mu_Earth);
    J(1) = getColJ(rp_xhi,rp_xlo,dv);
    % Perturb vy
    rp_yhi = radiusPeriapsis(v_inf_arr+[0,dv,0],v_inf_dep_earth,mu_Earth);
    rp_ylo = radiusPeriapsis(v_inf_arr-[0,dv,0],v_inf_dep_earth,mu_Earth);
    J(2) = getColJ(rp_yhi,rp_ylo,dv);
    % Perturb vz
    rp_zhi = radiusPeriapsis(v_inf_arr+[0,0,dv],v_inf_dep_earth,mu_Earth);
    rp_zlo = radiusPeriapsis(v_inf_arr-[0,0,dv],v_inf_dep_earth,mu_Earth);
    J(3) = getColJ(rp_zhi,rp_zlo,dv);

    v_inf_arr = v_inf_arr - pinv(J)*dr;

    iter = iter + 1;
end

rp = radiusPeriapsis(v_inf_arr,v_inf_dep_earth,mu_Earth);
turning_angle = rad2deg(angleBetween(v_inf_arr,v_inf_dep_earth))
%v_inf_flyby_hci = v_inf_arr + v2_earth % at time

% pick a time, find eph of earth, find required hci velocity, subtracxt
% current velocdity

%% Simulate heliocentric trajectory

t2_minus_dv = 10;
t_dv = t2_earth_dep - t2_minus_dv; % fire impulse 7 days before arrival at Earth

tspan_to_dv = seconds(days(t1_mars_dep:.5:t_dv));
tspan_dv_to_end = seconds(days(t_dv:1/24/6:t2_earth_dep));
% todo: could find time to periapsis
%tspan_dv_to_end = seconds(days(t_dv:1:t3_sun_arr));

ic_mars = [r1_mars,v1_mars_dep_hci];
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_traj,traj_to_dv] = ode113(@fode,tspan_to_dv,ic_mars,options,mu_Sun);

% figure(); hold on
% scatter3(r1_mars(1),r1_mars(2),r1_mars(3),36,mars_red,'filled','Marker','o')
% scatter3(r2_earth(1),r2_earth(2),r2_earth(3),36,earth_blue,'filled','Marker','o')
% plot3(traj1(:,1),traj1(:,2),traj1(:,3))

% Position and velocity right before delta-V execution
pos_b = traj_to_dv(end,1:3);
vel_b = traj_to_dv(end,4:6);

% Calculate delta-V required in heliocentric frame
[r_earth_dv,v_earth_dv] = planet_ephemeris_hci(t_dv,'Earth');

% Lambert solver for position at dv to position at earth -- rp?
e = sqrt(1+ (rp_des^2*norm(v_inf_arr)^4)/mu_Earth^2);
psi = acos(-1/e);
rp_vec = [rp_des*cos(pi-psi); rp_des*sin(pi-psi); 0]; % in ECI
rp_vec_hci = eci2hci(rp_vec)'+r2_earth; % position to target should be SOI

E = 1/2*norm(v_inf_arr)^2;
a = -mu_Earth / (2*E);
vp_mag = sqrt(2*mu_Earth/rp_des-mu_Earth/a);
vp_vec = cross([0;0;1],rp_vec/norm(rp_vec))*vp_mag;
vp_vec_hci = (eci2hci(vp_vec)+v2_earth')';

r_soi = 0.929E6;
v_des_hci = v_inf_arr+v2_earth;
v = -v_des_hci / norm(v_des_hci);
b = rp_des*vp_mag/norm(v_inf_arr);
offset_vec = b*cross([0 0 1],v)/norm(cross([0 0 1],v));
r_tgt = r_soi*-v_des_hci / norm(v_des_hci) + r2_earth + offset_vec;
%r_tgt = rp_vec_hci;

% Vary tof
% E = 1/2*norm(v_inf_arr)^2;
% a = -mu_Earth / (2*E);
% nuf = -acos(a*(1-e^2)/(norm(rp_vec_hci)*e)-1/e);
% nui = -acos(a*(1-e^2)/(norm(pos_b)*e)-1/e);
% n = sqrt(mu_Earth/abs(a)^3);
% tof = (nu2M_hyp(nuf,e) - nu2M_hyp(nui,e))/n;
tof = seconds(days(t2_earth_dep - t_dv));
% [v1,v2] = AA279lambert_curtis(mu_Sun,pos_b,rp_vec_hci,'pro',nrev,tof);

% Shoot until v_infs match
% tol = .01; % [km/s]
% max_iter = 100;
% iter = 0;
% dtof = 60; % [s]
% dvinf = Inf;
% while any(abs(dvinf)) > tol && iter < max_iter 
%     [v1,v2] = AA279lambert_curtis(mu_Sun,pos_b,r_tgt,'pro',nrev,tof);
%     vinf = real(v2) - v2_earth;
%     dvinf = (vinf - v_inf_arr)';
%     % Perturb tof
%     [v1hi,v2hi] = AA279lambert_curtis(mu_Sun,pos_b,r_tgt,'pro',nrev,tof+dtof);
%     vinfhi = real(v2hi) - v2_earth;
%     [v1lo,v2lo] = AA279lambert_curtis(mu_Sun,pos_b,r_tgt,'pro',nrev,tof-dtof);
%     vinflo = real(v2lo) - v2_earth;
%     J = getColJ(vinfhi,vinflo,dtof);
% 
%     tof = tof - pinv(J)*dvinf;
% 
%     iter = iter + 1;
% end
% 
% [v1,v2] = AA279lambert_curtis(mu_Sun,pos_b,r_tgt,'pro',nrev,tof);

%tof = seconds(days(28:1/24/60:31));
tof = seconds(days(t2_minus_dv-2:1/24:t2_minus_dv+2));
tol = .5; % [km/s]
deltaVs = NaN(length(tof),1);

for i=1:length(tof)
    [v1,v2] = AA279lambert_curtis(mu_Sun,pos_b,r_tgt,'pro',nrev,tof(i));
    v1 = real(v1);
    v2 = real(v2);
    vinf = v2-v2_earth;
    dvinf = vinf - v_inf_arr;
    % if all(abs(dvinf) < tol)
    %     disp('break');
    %     break;
    % end
    % dvp = v2 - vp_vec_hci;
    % if all(abs(dvp) < tol)
    %     break;
    %     disp('break');
    % end
    if all(abs(dvinf) < tol)
        avoidance_maneuver = v1 - vel_b;
        deltaVs(i) = norm(avoidance_maneuver);
    end
end

v1 = real(v1);
v2 = real(v2)
%vinf = v2-v2_earth
[m,i] = min(deltaVs);
tof = tof(i);
[v1,v2] = AA279lambert_curtis(mu_Sun,pos_b,r_tgt,'pro',nrev,tof);
v1 = real(v1);
v2 = real(v2);
vinf = v2-v2_earth;
dvinf = vinf - v_inf_arr

%tspan_b = linspace(0,tof(i)+3600*24*2,100);

% figure()
% plot(tof,deltaVs)

%%
avoidance_maneuver = v1 - vel_b;
deltaV = norm(avoidance_maneuver)
%[t,traj2] = ode113(@fode,tspan_dv_to_end,[pos_b1,v_hci],options,mu_Sun);
% [t,traj2_burn] = ode113(@fode,tspan_dv_to_end,[pos_b,vel_b+(v_arr_hci - vel_b)],options,mu_Sun);

[t,traj2_burn] = ode113(@fode,tspan_dv_to_end,[pos_b,vel_b+avoidance_maneuver],options,mu_Sun);
[t,traj2_capsule] = ode113(@fode,tspan_dv_to_end,[pos_b,vel_b],options,mu_Sun);

% Simulate departure from Earth fly-by
tspan3 = seconds(days(t2_earth_dep:1:t3_sun_arr));
[t,traj3_flyby] = ode113(@fode,tspan3,[traj2_burn(end,1:3),v2_earth_dep_hci],options,mu_Sun);



% Plot verification of delta-V
figure(); hold on
purple = [0.4940 0.1840 0.5560];
scatter3(r1_mars(1),r1_mars(2),r1_mars(3),36,mars_red,'filled','Marker','o')
scatter3(r2_earth(1),r2_earth(2),r2_earth(3),36,earth_blue,'filled','Marker','o')
plot3(traj_to_dv(:,1),traj_to_dv(:,2),traj_to_dv(:,3),'Color',purple)
scatter3(pos_b(1),pos_b(2),pos_b(3)) % position where the delta-V was fired
plot3(traj2_capsule(:,1),traj2_capsule(:,2),traj2_capsule(:,3),'--','Color',purple)
plot3(traj2_burn(:,1),traj2_burn(:,2),traj2_burn(:,3),'magenta')
scatter3(0,0,0,72,sun_yellow,'filled','Marker','o')
scatter3(r3_sun(1),r3_sun(2),r3_sun(3),'g','Marker','x')
% plot3(traj3_flyby(:,1),traj3_flyby(:,2),traj3_flyby(:,3))
plot3(traj23(:,1),traj23(:,2),traj23(:,3),'Color',purple);
plot3(sun_orbit(:,1),sun_orbit(:,2),sun_orbit(:,3),'Color',[0.3010 0.7450 0.9330]);
scatter3(r_tgt(1),r_tgt(2),r_tgt(3),72,'green','Marker','square')
axis equal
grid on

% Check closest approach
closest_approach = min( vecnorm(traj2_burn(:,1:3) - r2_earth,2,2) )


%% Simulate gravity assist
% Earth-centered coordinates and reference frame
% v_inf_arr_earth = traj2(end,4:6) - v2_earth;
% r_arr_earth = traj2(end,1:3) - r2_earth;
% 
% tspan = seconds(days(0:1/24:24));
% [t,traj_e] = ode113(@fode,tspan,[r_arr_earth,v_inf_arr_earth],options,mu_Earth);
% 
% figure(); hold on
% plot3(traj_e(:,1),traj_e(:,2),traj_e(:,3),'Color','r');
% scatter3(0,0,0,36,earth_blue,'filled','Marker','o');





%% Helpers

function rp = radiusPeriapsis(v_inf_arr,v_inf_dep,mu)
% Computes the radius of closest approach from v_inf of a hyperbolic
% trajectory
    delta = acos( dot(v_inf_arr,v_inf_dep) / (norm(v_inf_arr)*norm(v_inf_dep)) );
    rp = mu/norm(v_inf_arr)^2 * (1/cos((pi-delta)/2) - 1); % from Vallado
end

function col = getColJ(high,low,delta)
col = (high - low)'/(2*delta);
end


% % Energy of the incoming hyperbola at Earth
% v_inf = sqrt(norm(v_inf_arr_earth)*norm(v_inf_dep_earth));
% E = (v_inf^2)/2;
% a = -mu_Earth/(2*E); % [km]

% % Want desired offset distance
% b = R_Earth*5;
% 
% e = sqrt(1+b^2*v_inf^4/mu_Earth^2);
% psi = acos(-1/e) % angle between rp and v after passage
% 
% % What is the periapsis
% rp_mag = (e-1)*mu_Earth/v_inf^2
% rp_vec = [rp_mag*cos(pi-psi); rp_mag*sin(pi-psi); 0];
% vp_mag = sqrt(mu_Earth*(2/rp_mag-1/a))
% vp_vec = -cross(rp_vec/norm(rp_vec),[0;0;1])*vp_mag
% 
% % vp in HCI
% vp_hci = vp_vec + v2_earth;
% 
% % Delta-V needed to be applied at periapsis in order to hit hyperbolic traj
% dV = norm(v_inf_dep_earth - vp_hci)

% Shoot from SOI until we find the tof corresponding to this hyperbola with
% specified rp / apply delta-V to hit this point
% r_earth_soi = earth_soi
% [v_arr_earth,vp_lambert] = AA279lambert_curtis(mu_Earth,r2_earth,r3_sun,'pro',1,dt23);


%%
% With v1_earth_dep, find turning angle between the two vectors
% delta = angleBetween(v_inf_arr_earth,v_inf_dep_earth);
% 
% % Solve for rp by numerical iteration
% rp0 = 10e6; % [km] initial guess
% v_arr = norm(v_inf_arr_earth);
% v_dep = norm(v_inf_dep_earth);
% % Limited to certain angle
% turning_angle = @(rp) delta - ( asin(1/(1+rp*v_arr^2/mu_Earth)) + asin(1/(1+rp*v_dep^2/mu_Earth)) );
% options = optimset('TolX',1e-12);
% [rp_mag,fval,exitflag] = fzero(turning_angle,rp0,options);

% Bisection method
% rp_max = 2*mu_Earth / max(v_arr,v_dep)^2;
% low = 0;
% high = rp_max;
% while abs(high-low)/2 > 1e-9
%     c = (low+high)/2;
%     if sign(turning_angle(c)) == sign(turning_angle(low))
%        low = c; 
%     else
%        high = c;
%     end
% end
% rp = (low+high)/2;
% Check for convergence
% if abs(turning_angle(rp)) > 0.001
%     converge = false;
% else
%     converge = true;
% end


%% Compute delta-V applied at periapsis to put on Sun trajectory
% vB_arr = sqrt(v_arr^2 + 2*mu_earth/rp);
% vB_dep = sqrt(v_dep^2 + 2*mu_earth/rp);
% 
% delta_m = 2*asin(1/(1+rp*v_arr^2/mu_Earth));
% dv_p = abs(vB_dep - vB_arr) % [km/s]
% 
% % Energy of the incoming hyperbola at Earth
% E = (norm(v2_earth_arr_hci)^2)/2;
% a = -mu_Earth/(2*E);
% vp = sqrt(mu_Earth*(2/rp-1/a));

%deltaV = vp-norm(v0);

% u = [0,v_inf_arr_earth*cos(delta_m/2),0];
% u = u / norm(u);
% dv_p_hci = eci2hci([0;dv_p;0])' + v2_earth; % hci

%% Simulation
% ic = [r2_earth, dv_p_hci];
% tspan = linspace(0,seconds(days(t3_sun_arr-t2_earth_dep)),100);
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% [t_traj,traj_e2s] = ode113(@fode,tspan,ic,options,mu_Sun);
% 
% sun_yellow = [227,193,2]/255;
% figure(); hold on
% plot3(traj_e2s(:,1),traj_e2s(:,2),traj_e2s(:,3))
% scatter3(r2_earth(1),r2_earth(2),r2_earth(3),36,earth_blue,'filled','Marker','o')
% scatter3(r3_sun(1),r3_sun(2),r3_sun(3),72,sun_yellow,'filled','Marker','o')
% 
% % steps = 100;
% % for k=1:steps
% %     [t_traj,traj_e2s] = ode113(@fode,0:5:10 + tspan(k),ic,options,mu_Sun);
% %     ic = traj_e2s(end,:)'; % take the last step for next iteration
% %     % Check if at periapsis
% %     r_e = planet_ephemeris_hci(t2_earth_dep+)
% % end
% 
function r_hci = eci2hci(r_eci)
    theta = deg2rad(23.45); % Earth-ecliptic angle

    % Transformation matrix (clockwise rotation about x)
    T_eci2hci = [1 0 0;
                 0 cos(theta) sin(theta);
                 0 -sin(theta) cos(theta)];
    r_hci = (T_eci2hci*r_eci);
end