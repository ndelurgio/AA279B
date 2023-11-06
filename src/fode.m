function xdot = fode(t,x,mu)
xdot = [x(4:6); -mu/norm(x(1:3))^3 * x(1:3)];
end