function glat = geoc_lat(lat)
% glat = geoc_lat(lat)
% Geographic -> geocentric latitude (degrees), as idagrn6 does before its
% spherical distances and azimuths (azimth.f: flattening 1/298.25, same
% series), so that the distances of the moveout model (event_G_spec.m)
% match those of the idagrn6 reference and synthetic data.
e = 1/298.25;
r = lat*pi/180; w = sin(r);
s = ((2 - e)*w + 4*e*w.^3)*e.*cos(r);
glat = (r - s)*180/pi;
end
