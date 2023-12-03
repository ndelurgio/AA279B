function [v_inf,c3,r1_mars_hci,v1_mars_hci,r2_earth_hci,v2_earth_hci] = mars2earth_traj(tl,ta,mu_Sun)
% Computes the v_inf, c3 of an interplanetary trajectory for given launch
% and arrival dates and produces a heliocentric plot
    dt_sec = (ta-tl)*24*3600;
    [r1_mars_hci,v1_mars_hci] = planet_ephemeris_hci(tl,'Mars');
    [r2_earth_hci,v2_earth_hci] = planet_ephemeris_hci(ta,'Earth'); % arrival
    
    nrev = 0; % number of revolutions / phasing
    [v1_dep,v2_arr] = AA279lambert_curtis(mu_Sun,r1_mars_hci,r2_earth_hci,'pro',nrev,dt_sec);
    
    v_inf = norm(v1_dep - v1_mars_hci);
    c3 = v_inf^2;
    
    % Propagate orbit with Lambert initial conditions
    ic = [r1_mars_hci,v1_dep];
    
    % Simulation settings
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    n_steps = 1000;
    tspan = linspace(0,dt_sec,n_steps);
    [t,traj] = ode113(@fode,tspan,ic,options,mu_Sun);
    
    % Earth and Mars positions
    tspan_d = linspace(tl,ta,n_steps);
    r_mars = planet_ephemeris_hci(tspan_d','Mars');
    r_earth = planet_ephemeris_hci(tspan_d','Earth');
    
    % Trajectory visualization
    earth_blue = [11,98,212]/255;
    mars_red = [252,80,3]/255;
    
    figure('Name','Trajectory'); visualize_heliocentric_3d()
    scatter3(r1_mars_hci(1),r1_mars_hci(2),r1_mars_hci(3),36,mars_red,'filled','Marker','o')
    scatter3(r2_earth_hci(1),r2_earth_hci(2),r2_earth_hci(3),36,earth_blue,'filled','Marker','o')
    plot3(traj(:,1),traj(:,2),traj(:,3))
    plot3(r_mars(:,1),r_mars(:,2),r_mars(:,3),'Color',mars_red)
    plot3(r_earth(:,1),r_earth(:,2),r_earth(:,3),'Color',earth_blue)
    legend('Sun','Mars','Earth','Trajectory')
    axis equal
    
    % Date of min. energy departure and arrival
    dt_dep = datetime(tl,'ConvertFrom','juliandate');
    dt_arr = datetime(ta,'ConvertFrom','juliandate');
    
    disp('Departure date:'); disp(dt_dep)
    disp('Arrival date:'); disp(dt_arr)
    disp('Time of flight (days):'); disp(ta-tl);
    disp('V_inf (km/s):'); disp(v_inf);
    disp('C3 (km^2/s^2:'); disp(c3);
    
    % Display text on plot
    % annstr = sprintf('TOF (days): %0.3g \n $v_{inf}$ (km/s): %0.3g \n C3 (km$^2$/s$^2$): %0.3g', ...
    %     ta_min-tl_min,v_inf,c3'); % annotation text
    % annpos = [0.7 0.5 0.1 0.1]; % annotation position in figure coordinates
    % ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
    % ha.HorizontalAlignment = 'left';
    % ha.BackgroundColor = [1 1 1];


end