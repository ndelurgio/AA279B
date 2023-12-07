%% MARS ORBIT TO INJECTION
% Beginning with a nominal low Mars parking orbit, computes the delta-V
% required for injection to a desired v_inf from interplanetary Lambert
% solver
% Requires running mars_to_earth_traj_design 

mars_to_earth_traj_design;

%% Initialize orbit
alt = 325; % [km] source: ESA
ecc = 0; % circular orbit
inc = deg2rad(30);
raan = 0;
aop = 0;
nu = 0; % true anomaly

% Compute position and velocity from orbital elements
sma = alt+R_Mars; % [km]
p = sma*(1-ecc^2);

% Position and velocity in Mars centered inertial frame
[r0,v0] = orb2rv(p,ecc,inc,raan,aop,wrapTo2Pi(nu),0,0,0,mu_Mars);

% Circular orbit: periapsis = a
rp = alt+R_Mars;

%% Loop through all v1s
% Escape velocity
v_esc = sqrt(2*mu_Mars/R_Mars);
deltaVs = NaN*ones(length(tl_range),length(ta_range));

for i=1:length(tl_range)
    for j=1:length(ta_range)
        % 
        % if tl_range(i)==tl_min_dv && ta_range(j)==ta_min_dv
        %     disp('')
        % end

        %v_dep = v1_dep; % from interplanetary analysis
        v_dep = reshape(v1_data(i,j,:),1,[]);

        % Relative velocity of spacecraft with respect to Mars
        v_inf_vec = v_dep - v1_mars_hci;
        v_inf = norm(v_inf_vec); % planetocentric velocity

        % Energy of the hyperbola at Mars
        E = (v_inf^2)/2;
        a = -mu_Mars/(2*E);
        vp = sqrt(mu_Mars*(2/rp-1/a));

        deltaVs(i,j) = vp-norm(v0);

    end
end

%% Plot contours
figure('Name','Delta-V porkchop plot')
[x,y] = meshgrid(tl_range,ta_range);
dv_lvls = [0:1:10]; %[0:.2:4,5:2:40];
[dv_cont,dvh_cont] = contour(x,y,deltaVs',dv_lvls,'ShowText',true,"FaceAlpha",0.3,'Color','b'); hold on;
clabel(dv_cont,dvh_cont,dv_lvls,'Color','b');
c3_lvls = 0:1:20;
teal = [20,180,180]/255;
[c3_cont,c3h_cont] = contour(x,y,c3s',c3_lvls,'ShowText',true,"FaceAlpha",0.3,'Color',teal);
clabel(c3_cont,c3h_cont,c3_lvls,'Color',teal);
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

% Set launch period
zoom2launch = false;

if (zoom2launch)
    period_start = date2JD('2033-01-15 00:00:00');
    period_duration = 20; % [days]
    
    % Zoom in on launch period
    period_end = period_start + period_duration;
    xlim([period_start,period_end])

    % Reset tick marks
    xticks(linspace(period_start,period_end,5))
    xtick_zoomed = linspace(period_start,period_end,5);
    dt_zoomed = datetime(xtick_zoomed,'ConvertFrom','juliandate');
    xticks(xtick_zoomed)
    xticklabels(string(dt_zoomed))
end

% Min. energy arrival
scatter(tl_min,ta_min,75,[0.6350 0.0780 0.1840],'Marker','x','LineWidth',2); % marker at min. energy c3
%yline(ta_min,'r--','LineWidth',1) % line at min. energy arrival time
min_energy_c3 = c3s(i_min,j_min);
min_energy_dv = deltaVs(i_min,j_min);

% Min. dV
min_dv = min(min(deltaVs));
[i_min_dv,j_min_dv] = find(deltaVs==min_dv);
tl_min_dv = tl_range(i_min_dv);
ta_min_dv = ta_range(j_min_dv);
scatter(tl_min_dv,ta_min_dv,75,'magenta','x','LineWidth',2); % marker at min. delta-V
%yline(ta_min_dv,'--','Color','magenta','LineWidth',1) % line at min. delta-V arrival time
% min_dv_c3 = c3s(i_min_dv,j_min_dv);
% min_dv_dv = deltaVs(i_min_dv,j_min_dv);

%legend('Delta-V','C3','Min. energy','','Min. delta-V','','Location','northwest')
legend('Delta-V','C3','Location','northwest')


%% Examine data for fixed arrival date
% date_arr = ta_min; %date2JD('2033-')
%arr_datetext = datestr(datetime(date_arr,'ConvertFrom','juliandate'));
% idx = find(ta_range==date_arr);
% c3_arr = c3s(:,idx);
% dV_arr = deltaVs(:,idx);
% 
% figure()
% yyaxis left
% plot(tl_range,c3_arr);
% yyaxis right
% plot(tl_range,dV_arr);
% legend('c3','delta-V');
% xlim([min(tl_range),max(tl_range)])
% title(['Time of arrival: ',datestr(datetime(date_arr,'ConvertFrom','juliandate'))])

%% Map minimum delta-V trajectory 
[min_dv_v_inf,min_dv_c3,r1_mars_hci,v1_mars_hci,r2_earth_hci,v2_earth_hci] = mars2earth_traj(tl_min_dv,ta_min_dv,mu_Sun);

%% Injection velocity at a specified date

tl = tl_min_dv;
ta = ta_min_dv;
i = find(tl_range==tl); j = find(ta_range==ta);

dV = deltaVs(i,j);
v_inj = dV + norm(v0);
v_inj_pqw = [0;v_inj;0];

R_raan = [cos(raan) -sin(raan) 0;
          sin(raan) cos(raan) 0;
          0 0 1];
R_i = [1 0 0;
       0 cos(inc) -sin(inc);
       0 sin(inc) cos(inc)];
R_aop = [cos(aop) -sin(aop) 0;
         sin(aop) cos(aop) 0;
         0 0 1];

v_inj_vec = (R_raan*R_i*R_aop)*v_inj_pqw;



%% Simulate Mars orbit and visualization

% FODE Mars
% Simulation settings
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% Orbital period
T = 2*pi*sqrt(sma^3/mu_Mars);
n_steps = 50;
tspan = linspace(0,T,n_steps);
[t,traj] = ode113(@fode,tspan,[r0;v0],options,mu_Mars);
r_prk = traj(:,1:3);

ic = [r0;v_inj_vec];
[t,traj_inj] = ode113(@fode,tspan,ic,options,mu_Mars);

% Visualize orbit
figure('Name','Mars injection visualization');
visualize_planet('Mars'); hold on; 
plot3(r_prk(:,1),r_prk(:,2),r_prk(:,3),'Color',mars_red,'LineWidth',1);
plot3(traj_inj(:,1),traj_inj(:,2),traj_inj(:,3),'Color',mars_red,'LineWidth',1);
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
xlim([-5e3 5e3])
ylim([-5e3 5e3])
zlim([-5e3 5e3])
view([434 12])

