function [pos_ti,t_duration,iter] = reentry_shooter_hyperbola(pos_ti,t_duration,vel_ti_des,pos_tf,mu_earth)

dp = 100000;
dt = 60;
max_iter = 100;

J = zeros(3,4);
for iter = 1:max_iter
    % Nominal run to get error
    [vi,~,~] = AA279lambert_curtis(mu_earth, pos_ti, pos_tf, 'retro', 0, t_duration);
    dv = vi-vel_ti_des;
    disp(norm(dv))
    % End if converged
    if norm(dv) < 1
        break
    end
    % Perturb px
    [vi_high,~,~] = AA279lambert_curtis(mu_earth, pos_ti + [dp,0,0], pos_tf, 'retro', 0, t_duration);
    [vi_low,~,~] = AA279lambert_curtis(mu_earth, pos_ti - [dp,0,0], pos_tf, 'retro', 0, t_duration);
    J(:,1) = getColJ(vi_high,vi_low,dp);
    % Perturb py
    [vi_high,~,~] = AA279lambert_curtis(mu_earth, pos_ti + [0,dp,0], pos_tf, 'retro', 0, t_duration);
    [vi_low,~,~] = AA279lambert_curtis(mu_earth, pos_ti - [0,dp,0], pos_tf, 'retro', 0, t_duration);
    J(:,2) = getColJ(vi_high,vi_low,dp);
    % Perturb pz
    [vi_high,~,~] = AA279lambert_curtis(mu_earth, pos_ti + [0,0,dp], pos_tf, 'retro', 0, t_duration);
    [vi_low,~,~] = AA279lambert_curtis(mu_earth, pos_ti - [0,0,dp], pos_tf, 'retro', 0, t_duration);
    J(:,3) = getColJ(vi_high,vi_low,dp);
    % Perturb t
    [vi_high,~,~] = AA279lambert_curtis(mu_earth, pos_ti, pos_tf, 'retro', 0, t_duration+dt);
    [vi_low,~,~] = AA279lambert_curtis(mu_earth, pos_ti, pos_tf, 'retro', 0, t_duration-dt);
    J(:,4) = getColJ(vi_high,vi_low,dt);
    % Update nominal
    update = [pos_ti';t_duration] - pinv(J)*dv';
    pos_ti = update(1:3)';
    t_duration = update(4);
    % ti = tf - seconds(t_duration);
end
end

function col = getColJ(high,low,delta)
col = (high - low)'/(2*delta);
end