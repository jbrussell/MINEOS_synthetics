function [arclen, az] = distance(lat1, lon1, lat2, lon2)
% [arclen, az] = distance(lat1, lon1, lat2, lon2)  (fallback when the
% mapping toolbox/package is missing) Great-circle distance (degrees) and
% azimuth from point 1 to point 2 (degrees clockwise from north) on a
% sphere, elementwise with implicit expansion.
p1 = lat1*pi/180; p2 = lat2*pi/180; dl = (lon2 - lon1)*pi/180;
h = sin((p2 - p1)/2).^2 + cos(p1).*cos(p2).*sin(dl/2).^2;
arclen = 2*asin(min(1, sqrt(h))) * 180/pi;
az = mod(atan2(sin(dl).*cos(p2), cos(p1).*sin(p2) - sin(p1).*cos(p2).*cos(dl)) * 180/pi, 360);
end
