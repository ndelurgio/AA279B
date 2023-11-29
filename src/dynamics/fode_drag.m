function xdot = fode_drag(t,x,mu,spacecraft,t_epoch,space_weather_data)
date = t_epoch + seconds(t);
LLA = eci2lla_datetime(x(1:3)',date);
% Density Model
warning('off','aero:atmosnrlmsise00:invalidAltitude');
[f107average,f107daily,magneticIndex] = fluxSolarAndGeomagnetic(...
    date.Year,day(date,'dayofyear'),...
    date.Second,space_weather_data);
[~, rho] = atmosnrlmsise00(LLA(3),LLA(1),LLA(2),...
    date.Year,day(date,'dayofyear'),date.Second,...
    f107average,f107daily,magneticIndex);
rho = rho(6); % This is the total parameter
warning('on','aero:atmosnrlmsise00:invalidAltitude');
% Drag Perturbation
w_earth = [0;0;7.2921150E-5];
v_rel = x(4:6) - cross(w_earth,x(1:3));
a_drag = -1/2*rho*norm(v_rel)*spacecraft.B*v_rel;
% ODE
xdot = [x(4:6); -mu/norm(x(1:3))^3 * x(1:3) + a_drag];
end