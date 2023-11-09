function [vel_ti,t_duration,iter] = reentry_shooter_dv(vel_ti,t_duration,pos_ti,pos_tf,ti,tf,capsule,mu_earth)
space_weather_data = aeroReadSpaceWeatherData('scenarios/data/SW-Last5Years.csv');
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
time_step = 1;

dv = 1;
dt = 1;
max_iter = 100;

run_ode = @(vi,tf,ti_curr) ode113(...
    @fode_drag,...
    0:time_step:tf,...
    [pos_ti,vi],...
    options,mu_earth,capsule,ti_curr,space_weather_data);

J = zeros(3,4);
for iter = 1:max_iter
    % Get current range inertial position
    % pos_tf = lla2eci(range_LLA,datetime2vec(ti+seconds(t_duration)));
    % Nominal run to get error
    [~,traj_nom]=run_ode(vel_ti,t_duration,ti);
    dr = traj_nom(end,1:3)-pos_tf;
    disp(norm(dr))
    % End if converged
    if norm(dr) < 1000
        break
    end
    % Perturb vx
    [~,traj_high]=run_ode(vel_ti + [dv,0,0],t_duration,ti);
    [~,traj_low]=run_ode(vel_ti - [dv,0,0],t_duration,ti);
    J(:,1) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dv);
    % Perturb vy
    [~,traj_high]=run_ode(vel_ti + [0,dv,0],t_duration,ti);
    [~,traj_low]=run_ode(vel_ti - [0,dv,0],t_duration,ti);
    J(:,2) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dv);
    % Perturb vz
    [~,traj_high]=run_ode(vel_ti + [0,0,dv],t_duration,ti);
    [~,traj_low]=run_ode(vel_ti - [0,0,dv],t_duration,ti);
    J(:,3) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dv);
    % Perturb t
    [~,traj_high]=run_ode(vel_ti,t_duration+dt,ti);
    [~,traj_low]=run_ode(vel_ti,t_duration-dt,ti);
    J(:,4) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dt);
    % Update nominal
    update = [vel_ti';t_duration] - pinv(J)*dr';
    vel_ti = update(1:3)';
    t_duration = update(4);
    ti = tf - seconds(t_duration);
end
end

function col = getColJ(high,low,delta)
col = (high - low)'/(2*delta);
end