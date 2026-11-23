function km = deg2km(deg)
% km = deg2km(deg)  (fallback when the mapping toolbox/package is missing)
% Great-circle degrees to km on a sphere of radius 6371 km.
km = deg * 6371 * pi/180;
end
