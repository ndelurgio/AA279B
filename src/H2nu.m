function nu = H2nu(H,e)
% nu = atan2(-sinh(H)*sqrt(e^2-1),(cosh(H)-e));
% nu = acos((cosh(H)-e)/(1-e*cosh(H)));
nu = wrapToPi(2*atan2(sqrt((e+1)/(e-1))*tanh(H/2),1));
end