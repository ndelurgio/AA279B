function H = M2H(M,e)
% Initial Guess
if e < 1.6
    if (M > -pi && M < 0) || (M > pi && M < 2*pi)
        H = M - e;
    else
        H = M + e;
    end
else
    if e < 3.6 && abs(M) > pi
        H = M - sign(M)*e;
    else
        H = M/(e-1);
    end
end
% Iteration
tol = 1e-9;
H_prev = inf;
while abs(H-H_prev) > tol
    H_prev = H;
    H = H + (M-e*sinh(H)+H)/(e*cosh(H)-1);
end
end