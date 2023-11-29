function [t_traj,traj] = fast_unperturbed_sim(t_duration,pos_ti,vel_ti,mu_earth)
dt = 6;
% r_thresh = 6378e3;
[a,ecc,incl,RAAN,argp,nu,~,~,~] = rv2orb(pos_ti', vel_ti', mu_earth);
% nu_thresh = acos((a*(1-ecc^2)/r_thresh - 1)/ecc);
M_init = nu2M_hyp(nu,ecc);
% nu_thresh = nu_thresh * sign(M_init);
% M_thresh = nu2M_hyp(nu_thresh,ecc);
n = sqrt(mu_earth/abs(a)^3);
% t_thresh = (M_thresh - M_init)/n;
% Fast kepler sim
% t_traj = 0:dt:t_thresh;
t_traj = 0:dt:(t_duration);
traj = zeros(length(t_traj),6);
for i = 1:length(t_traj)
    M = M_init + n*t_traj(i);
    nu = wrapToPi(M2nu_hyp(M,ecc));
    [r,v] = orb2rv(a*(1-ecc^2)/1000,ecc,incl,RAAN,argp,wrapTo2Pi(nu));
    traj(i,:) = [r;v]'*1000;
end
% Drag-perturbed sim
% tspan = t_thresh:dt:t_duration;
% [t_traj_drag,traj_drag] = ode113(...
%     @fode_drag,...
%     linspace(t_thresh,t_duration,length(tspan)),...
%     traj(end,:),...
%     options,mu_earth,capsule,ti_utc + seconds(t_thresh),space_weather_data);
% t_traj = [t_traj'; t_traj_drag];
% traj = [traj;traj_drag];
end