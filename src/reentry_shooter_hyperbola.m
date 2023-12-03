function [pos_ti,t_duration,iter] = reentry_shooter_hyperbola(pos_ti,t_duration,vel_ti_des,pos_tf,mu_earth)

dp = 100000;
dt = 60;
max_iter = 100;

J = zeros(3,4);
for iter = 1:max_iter
    % Nominal run to get error
    [vi,vf,~] = AA279lambert_curtis(mu_earth, pos_ti, pos_tf, 'retro', 0, t_duration);
    % [~,~,~,~,~,nui,~,~,~] = rv2orb(pos_ti', vi', mu_earth);
    % nui = wrapToPi(nui);
    % [a,ecc,incl,RAAN,argp,nu,~,~,~] = rv2orb(pos_tf', vf', mu_earth);
    % nu = wrapToPi(nu);
    % if nu > 0
    %     nui = -nui;
    %     [p_new,vf]
    % end
    dv = vi-vel_ti_des;
    disp(norm(dv))
    % End if converged
    tol = norm(dv);
    % if wrapToPi(nu) > 0
    %     tol = tol + 
    % end
    if tol < 0.1
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
% [~,~,~,~,~,nuf,~,~,~] = rv2orb(pos_tf', vf', mu_earth);
% nuf = wrapToPi(nuf);
% if nuf > 0
%     [a,ecc,incl,RAAN,argp,nui,~,~,~] = rv2orb(pos_ti', vi', mu_earth);
%     nui = wrapToPi(nui);
%     [pos_ti,~] = orb2rv(a*(1-ecc^2),ecc,incl,RAAN,argp,wrapTo2Pi(-nui));
%     pos_ti = pos_ti';
% end
% [pos_ti,t_duration,iter] = reentry_shooter_hyperbola(pos_ti,t_duration,vel_ti_des,pos_tf,mu_earth);
end

% function [vi,vf,~] = lambert_wrapper(mu_earth, pos_ti, pos_tf, 0, t_duration)
%     [vi,vf,~] = AA279lambert_curtis(mu_earth, pos_ti, pos_tf, 'retro', 0, t_duration);
%     % [~,~,~,~,~,nuf,~,~,~] = rv2orb(pos_tf', vf', mu_earth);
%     % nuf = wrapToPi(nuf);
%     % if nuf > 0
%     %     [a,ecc,incl,RAAN,argp,nui,~,~,~] = rv2orb(pos_ti', vi', mu_earth);
%     %     nui = wrapToPi(nui);
%     %     [pos_ti,~] = orb2rv(a*(1-ecc^2),ecc,incl,RAAN,argp,wrapTo2Pi(-nui));
%     %     pos_ti = pos_ti';
%     % end
% end

function col = getColJ(high,low,delta)
col = (high - low)'/(2*delta);
end