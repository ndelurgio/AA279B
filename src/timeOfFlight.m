function tf = timeOfFlight(jd1,jd2)
% Return time of flight in seconds between 2 Julian dates, converting days to seconds 
    tf = (jd2-jd1)*24*3600;
end
