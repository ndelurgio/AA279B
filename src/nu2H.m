function H = nu2H(nu,e)
% H = atanh(sin(nu)*sqrt(e^2-1)/(e+cos(nu)));
H = asinh(sin(nu)*sqrt(e^2-1)/(1+e*cos(nu)));
end