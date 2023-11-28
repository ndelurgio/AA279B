function M = nu2M_hyp(nu,e)
H = nu2H(nu,e);
M = H2M(H,e);
% M = e*sqrt(e^2-1)*sin(nu)/(1+e*cos(nu))-...
%     log((sqrt(e+1)+sqrt(e-1)*tan(nu/2))/(sqrt(e+1)-sqrt(e-1)*tan(nu/2)));
end