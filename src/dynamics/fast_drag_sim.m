function [t_traj,traj] = fast_drag_sim(t_duration,pos_ti,vel_ti,mu_earth,capsule,ti_utc,space_weather_data,options)
dt = 6;
r_thresh = 1000e3;
[a,ecc,incl,RAAN,argp,nu,~,~,~] = rv2orb(pos_ti', vel_ti', mu_earth);
nu_thresh = acos((a*(1-ecc^2)/r_thresh - 1)/ecc);
M_thresh = nu2M_hyp(nu_thresh,ecc);
M_init = nu2M_hyp(nu,ecc);
n = sqrt(mu_earth/abs(a)^3);
t_thresh = (M_thresh - M_init)/n;
% Fast kepler sim
t_traj = 0:dt:t_thresh;
traj = zeros(length(t_traj),6);
for i = 1:length(t_traj)
    M = M_init + n*t_traj(i);
    nu = wrapToPi(M2nu_hyp(M,ecc));
    [r,v] = orb2rv(a*(1-ecc^2)/1000,ecc,incl,RAAN,argp,wrapTo2Pi(nu));
    traj(i,:) = [r;v]'*1000;
end
% Drag-perturbed sim
[t_traj_drag,traj_drag] = ode113(...
    @fode_drag,...
    t_thresh:dt:t_duration,...
    [pos_ti,vel_ti],...
    options,mu_earth,capsule,ti_utc + seconds(t_thresh),space_weather_data);
t_traj = [t_traj, t_traj_drag];
traj = [traj;traj_drag];
end