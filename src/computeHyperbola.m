function [a,e,i,Om,w] = computeHyperbola(r_target,v_inf,mu)
theta = acos(dot(v_inf,r_target)/(norm(v_inf)*norm(r_target)));
rt = norm(r_target);
energy = norm(v_inf)^2/2;
a = -mu/(2*energy);

e_max = 1 - rt/a;
f_ang = @(e) acos(-1/e) - acos(a*(1-e^2)/(rt*e)-1/e) + theta - pi;
low = 1;
high = e_max;
while abs(high-low)/2 > 1e-9
    c = (low+high)/2;
    if sign(f_ang(c)) == sign(f_ang(low))
       low = c; 
    else
        high = c;
    end
end
e_des = (low+high)/2;

x = v_inf'/norm(v_inf);
z = cross(x,r_target')/norm(cross(x,r_target'));
y = cross(z,x);
C_hyp2eci = [x,y,z];

rp_mag = a*(1-e_des);
psi = acos(-1/e_des);
rp_vec = [rp_mag*cos(pi-psi); rp_mag*sin(pi-psi); 0];
vp_mag = sqrt(2*mu/rp_mag-mu/a);
vp_vec = cross(rp_vec/norm(rp_vec),[0;0;1])*vp_mag;
capsule_pos_tp_j2000 = (C_hyp2eci*rp_vec)';
capsule_vel_tp_j2000 = (C_hyp2eci*vp_vec)';

[a,e,i,Om,w,~,~,~,~] = rv2orb(capsule_pos_tp_j2000', capsule_vel_tp_j2000', mu);
end