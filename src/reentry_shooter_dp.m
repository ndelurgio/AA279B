function [pos_ti,t_duration,iter] = reentry_shooter_dp(pos_ti,t_duration,vel_ti,pos_tf,ti,tf,capsule,mu_earth)
space_weather_data = aeroReadSpaceWeatherData('scenarios/data/SW-Last5Years.csv');
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
time_step = 1;

dp = 1000;
dt = 1;
max_iter = 10;

J = zeros(3,4);
for iter = 1:max_iter
    % Nominal run to get error
    [~,traj_nom]=fast_drag_sim(t_duration,pos_ti,vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    % [~,traj_nom]=run_ode(pos_ti,t_duration,ti);
    dr = traj_nom(end,1:3)-pos_tf;
    disp(norm(dr))
    % End if converged
    if norm(dr) < 1000
        break
    end
    % Perturb vx
    [~,traj_high]=fast_drag_sim(t_duration,pos_ti + [dp,0,0],vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    [~,traj_low]=fast_drag_sim(t_duration,pos_ti - [dp,0,0],vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    J(:,1) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dp);
    % Perturb vy
    [~,traj_high]=fast_drag_sim(t_duration,pos_ti + [0,dp,0],vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    [~,traj_low]=fast_drag_sim(t_duration,pos_ti - [0,dp,0],vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    J(:,2) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dp);
    % Perturb vz
    [~,traj_high]=fast_drag_sim(t_duration,pos_ti + [0,0,dp],vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    [~,traj_low]=fast_drag_sim(t_duration,pos_ti - [0,0,dp],vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    J(:,3) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dp);
    % Perturb t
    [~,traj_high]=fast_drag_sim(t_duration+dt,pos_ti,vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    [~,traj_low]=fast_drag_sim(t_duration-dt,pos_ti,vel_ti,mu_earth,capsule,ti,space_weather_data,options);
    J(:,4) = getColJ(traj_high(end,1:3),traj_low(end,1:3),dt);
    % Update nominal
    update = [pos_ti';t_duration] - pinv(J)*dr';
    pos_ti = update(1:3)';
    t_duration = update(4);
    ti = tf - seconds(t_duration);
end
end

function col = getColJ(high,low,delta)
col = (high - low)'/(2*delta);
end