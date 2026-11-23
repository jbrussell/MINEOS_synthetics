function deg = km2deg(km)
% deg = km2deg(km)  (fallback when the mapping toolbox/package is missing)
% km to great-circle degrees on a sphere of radius 6371 km.
deg = km / (6371 * pi/180);
end
