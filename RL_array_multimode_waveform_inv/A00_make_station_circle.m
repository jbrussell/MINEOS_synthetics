% A00_MAKE_STATION_CIRCLE  Data preparation, step 0 (optional, synthetic
% tests): make a station file with opt.circle_N stations evenly spaced in
% azimuth on a circle of diameter opt.circle_aperture_km (km) centred on
% the NoMelt array. Use it for a synthetic test by setting
% opt.stationfile to the new file (with a new opt.aswms_dir and
% opt.prepfile), then A0-A3 and D1-D4 as usual.
%
% Settings (parameters_RL.m, or opt before running):
%   opt.circle_N           number of stations (default 16, as NoMelt)
%   opt.circle_aperture_km diameter of the circle (km, default 1000): the
%                          distance across the array between opposite
%                          stations, along the great circle through the
%                          centre
%   opt.circle_file        station file to write ('' = default
%                          stations_circle_<aperture>km_N<N>.txt here)
%   opt.circle_ref         station file whose centre is used (default
%                          stations_nomelt.txt): the centre is the mean
%                          position of its stations (mean of their unit
%                          vectors); the new stations get the median of its
%                          elevations (km; negative = below sea level)
%
% Stations C01..CNN, clockwise from north (C01 due north of the centre), on
% the sphere (radius 6371 km) at great-circle distance aperture/2 from the
% centre, in the format of stations_nomelt.txt (name, lat, lon, elevation
% in km). A map of the reference and new stations is saved as
% figs_stations/A00_<file name>.png.
%
%   >> A00_make_station_circle                                   % 16 stations, 1000 km
%   >> opt.circle_N = 24; opt.circle_aperture_km = 2000; A00_make_station_circle
%   >> opt.stationfile = 'stations_circle_1000km_N16.txt'; opt.aswms_dir = 'ASWMS_circle/';
%   >> opt.prepfile = 'prep_circle.mat'; A0_synth_make_sac_data   % then A1-A3, D1-D4
%
% Next: A0_synth_make_sac_data.m (synthetic data at these stations).

parameters_RL;
here_ = fileparts(mfilename('fullpath'));
N_ = opt.circle_N; ap_ = opt.circle_aperture_km;
ref_ = opt.circle_ref;
cfile_ = opt.circle_file;
if isempty(cfile_), cfile_ = fullfile(here_, sprintf('stations_circle_%gkm_N%d.txt', ap_, N_)); end
% centre of the reference array: mean of the unit vectors
fid_ = fopen(ref_, 'r');
if fid_ < 0, error('A00_make_station_circle: cannot read %s', ref_); end
S_ = textscan(fid_, '%s %f %f %f'); fclose(fid_);
la_ = S_{2}*pi/180; lo_ = S_{3}*pi/180;
v_ = mean([cos(la_).*cos(lo_), cos(la_).*sin(lo_), sin(la_)], 1);
lat0 = atan2(v_(3), hypot(v_(1), v_(2)))*180/pi;
lon0 = atan2(v_(2), v_(1))*180/pi;
elev_ = median(S_{4});
% stations on the circle: destination points at angular distance a_ (rad)
R_ = 6371;
a_ = (ap_/2)/R_;
az_ = (0:N_-1)'*360/N_*pi/180;
p0_ = lat0*pi/180; l0_ = lon0*pi/180;
lat_ = asin(sin(p0_)*cos(a_) + cos(p0_)*sin(a_)*cos(az_));
lon_ = l0_ + atan2(sin(az_)*sin(a_)*cos(p0_), cos(a_) - sin(p0_)*sin(lat_));
lat_ = lat_*180/pi; lon_ = mod(lon_*180/pi + 180, 360) - 180;
names_ = arrayfun(@(k) sprintf('C%02d', k), (1:N_)', 'uniformoutput', false);
fid_ = fopen(cfile_, 'w');
if fid_ < 0, error('A00_make_station_circle: cannot write %s', cfile_); end
for k_ = 1:N_
    fprintf(fid_, '  %-5s %11.6f %13.6f %7.4f\n', names_{k_}, lat_(k_), lon_(k_), elev_);
end
fclose(fid_);
% check: distance from the centre and largest distance across the array
dc_ = gc_dist_az(lat0*ones(N_, 1), lon0*ones(N_, 1), lat_, lon_)*pi/180*R_;
[I_, J_] = meshgrid(1:N_);
dij_ = gc_dist_az(lat_(I_(:)), lon_(I_(:)), lat_(J_(:)), lon_(J_(:)))*pi/180*R_;
fprintf('%d stations on a circle of diameter %g km about (%.4f, %.4f) (centre of %s):\n', ...
        N_, ap_, lat0, lon0, ref_);
fprintf('  distance from the centre %.1f km, station spacing %.1f km, largest distance %.1f km\n', ...
        mean(dc_), ap_*sin(pi/N_)*(N_ > 1), max(dij_));
fprintf('  written to %s (set opt.stationfile to it)\n', cfile_);
% map
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
fig = numbered_figure(4000, 'position', [50 50 700 650]); hold on; box on;
h1_ = plot(S_{3}, S_{2}, 'k^', 'markerfacecolor', [0.6 0.6 0.6], 'markersize', 7);
h2_ = plot(lon_, lat_, 'ro', 'markerfacecolor', 'r', 'markersize', 6);
plot(lon0, lat0, 'k+', 'markersize', 12, 'linewidth', 1.5);
text(lon_, lat_, strcat({'  '}, names_), 'fontsize', 7);
daspect([1 cosd(lat0) 1]);
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
[~, fn_] = fileparts(cfile_);
legend([h1_ h2_], {sprintf('reference (%d)', numel(S_{2})), sprintf('circle, %g km (%d)', ap_, N_)}, ...
       'location', 'southoutside');
title(fn_, 'interpreter', 'none');
od_ = fullfile(here_, 'figs_stations'); if exist(od_, 'dir') ~= 7, mkdir(od_); end
save_figure(fig, fullfile(od_, [fn_ '.png']), 130);
if fig_close_mode, close(fig); end
clear here_ N_ ap_ ref_ cfile_ fid_ S_ la_ lo_ v_ elev_ R_ a_ az_ p0_ l0_ lat_ lon_ names_ k_ dc_ I_ J_ dij_ h1_ h2_ fn_ od_
