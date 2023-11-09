function jd = date2JD(datestring)
% Inputs:
%   datestring - date in UTC
% Outputs:
%   jd - Julian date in days
    D = datetime(datestring);
    jd = juliandate(D);
end