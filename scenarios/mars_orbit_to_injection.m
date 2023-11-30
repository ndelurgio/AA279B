%% MARS ORBIT TO INJECTION
% Beginning with a nominal low Mars parking orbit, computes the delta-V
% required for injection to a desired v_inf from interplanetary Lambert
% solver

%setup_constants;

%% Choose inclination to match hyperbolic traj. inclination
% v_p -> rv2orb -> inclin

% Initialize orbit
% alt = 325; % [km] source: ESA
% 
% % a = alt+R_Mars; % [km]
% % p = a*(1-ecc^2);
% 
% % from interplanetary analysis
% v_dep = v1_dep; 
% 
% % Relative velocity of spacecraft with respect to Mars
% v_inf_vec = v_dep - v1_mars_hci;
% v_inf = norm(v_inf_vec);
% 
% rp = -(alt+R_Mars)*v_inf_vec/v_inf; % assumes we launch at periapsis, directly opposite v_inf
% vp = sqrt(mu_Mars/norm(rp))*v_inf_vec/v_inf; % assuming circular Mars orbit
% 
% [a,eMag,inc,raan,aop,nu,truLon,argLat,lonPer,p] = rv2orb(rp,vp,mu_Mars);


%% Initialize orbit
alt = 325; % [km] source: ESA
ecc = 0; % circular orbit
inc = deg2rad(1);
raan = 0;
aop = 0;
nu = 0; % true anomaly

% Compute position and velocity from orbital elements
a = alt+R_Mars; % [km]
p = a*(1-ecc^2);

% Position and velocity in Mars centered inertial frame
[r0,v0] = orb2rv(p,ecc,inc,raan,aop,wrapTo2Pi(nu),0,0,0,mu_Mars);

%%

% v_dep = v1_dep; % from interplanetary analysis
% 
% % Relative velocity of spacecraft with respect to Mars
% v_inf_vec = v_dep - v1_mars_hci;
% v_inf = norm(v_inf_vec);

% Escape velocity
v_esc = sqrt(2*mu_Mars/R_Mars);

v_inj = sqrt(v_inf^2+v_esc^2);

deltaV = v_inj - norm(v0);

%% Loop through all v1s
% Escape velocity
v_esc = sqrt(2*mu_Mars/R_Mars);

%v1s = [2 4 5; 13 5 6; 23 5 3]; % testing
deltaVs = NaN*ones(length(tl_range),length(ta_range));

for i=1:length(tl_range)
    for j=1:length(ta_range)

        %v_dep = v1_dep; % from interplanetary analysis
        v_dep = reshape(v1s(i,j,:),1,[]);

        % Relative velocity of spacecraft with respect to Mars
        v_inf_vec = v_dep - v1_mars_hci;
        v_inf = norm(v_inf_vec);
          
        v_inj = sqrt(v_inf^2+v_esc^2);
        
        deltaVs(i,j) = v_inj - norm(v0);
    end
end

%% Plot contours
figure('Name','Delta-V porkchop plot')
[x,y] = meshgrid(tl_range,ta_range);
dv_lvls = [0:.2:4,5:2:40];
%contour(x,y,deltaVs',dv_lvls)
[c_cont,h_cont] = contour(x,y,deltaVs',dv_lvls,'ShowText',true,"FaceAlpha",0.3,'Color','b'); hold on;
[c_cont,h_cont] = contour(x,y,c3s,c3_lvls,'ShowText',true,"FaceAlpha",0.3,'Color',[20,180,180]/255);
%clabel(c_cont,h_cont,2:2:20);
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




