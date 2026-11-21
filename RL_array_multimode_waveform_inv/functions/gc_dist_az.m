function [d, az] = gc_dist_az(lat1, lon1, lat2, lon2)
% [d, az] = gc_dist_az(lat1, lon1, lat2, lon2)
% Great-circle distance d (degrees) and azimuth az (degrees clockwise from
% north) from point 1 to point 2 on a sphere. Inputs in degrees; arrays
% are expanded against each other.
p1 = lat1(:)*pi/180; p2 = lat2(:)*pi/180; dl = (lon2(:) - lon1(:))*pi/180;
c = sin(p1).*sin(p2) + cos(p1).*cos(p2).*cos(dl);
d = acos(max(-1, min(1, c)))*180/pi;
az = mod(atan2(sin(dl).*cos(p2), cos(p1).*sin(p2) - sin(p1).*cos(p2).*cos(dl))*180/pi, 360);
end
