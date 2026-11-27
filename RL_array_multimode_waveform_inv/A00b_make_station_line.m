% A00B_MAKE_STATION_LINE  Data preparation, step 0 (optional, synthetic
% tests): make a station file with stations on a straight line along the
% great circle from one event (opt.line_event) through the centre of the
% NoMelt array, and a catalogue file with that event only - the geometry
% of the single-event line test: every station on the great circle, so
% the waves arrive straight along the line (arrival-angle anomaly 0) and
% the moveout across the line constrains the phase velocities.
%
% Settings (parameters_RL.m, or opt before running):
%   opt.line_event      Global CMT id of the event in opt.eventfile
%                       (default '201210280304A', Haida Gwaii Mw 7.8,
%                       44.7 deg from NoMelt: all branches 0-3 of Love and
%                       Rayleigh radiate toward NoMelt at 0.62-0.91 of
%                       their maximum over azimuth)
%   opt.line_length_km  length of the line (km, default 2000), centred on
%                       the array centre
%   opt.line_spacing_km station spacing along the line (km, default 50;
%                       rounded so that the stations fill the line evenly)
%   opt.line_file       station file to write ('' = default
%                       stations_line_<event>_<length>km.txt here)
%   opt.line_ref        station file whose centre is used (default
%                       stations_nomelt.txt; mean of the unit vectors, as
%                       A00) and whose median elevation the stations get
% The event file events_<event>.cmt (the event's CMTSOLUTION block of
% opt.eventfile) is written next to it.
%
% Station names L<kk>, numbered from the end nearest the event; distances
% along the great circle on the sphere (radius 6371 km). A map is saved as
% figs_stations/A00b_<file name>.png.
%
%   >> A00b_make_station_line
%   >> opt.stationfile = 'stations_line_201210280304A_2000km.txt';
%   >> opt.eventfile = 'events_201210280304A.cmt'; opt.aswms_dir = 'ASWMS_line/';
%   >> opt.prepfile = 'prep_line.mat'; A0_synth_make_sac_data   % then A1-A3, D2-D4
%
% Next: A0_synth_make_sac_data.m (synthetic data at these stations).

parameters_RL;
here_ = fileparts(mfilename('fullpath'));
eid_ = opt.line_event; L_ = opt.line_length_km; ds_ = opt.line_spacing_km;
lfile_ = opt.line_file;
if isempty(lfile_), lfile_ = fullfile(here_, sprintf('stations_line_%s_%gkm.txt', eid_, L_)); end
% centre of the reference array: mean of the unit vectors (as A00)
fid_ = fopen(opt.line_ref, 'r');
if fid_ < 0, error('A00b_make_station_line: cannot read %s', opt.line_ref); end
S_ = textscan(fid_, '%s %f %f %f'); fclose(fid_);
la_ = S_{2}*pi/180; lo_ = S_{3}*pi/180;
v_ = mean([cos(la_).*cos(lo_), cos(la_).*sin(lo_), sin(la_)], 1);
lat0 = atan2(v_(3), hypot(v_(1), v_(2)))*180/pi;
lon0 = atan2(v_(2), v_(1))*180/pi;
elev_ = median(S_{4});
% the event
ev_ = rl_read_cmtsolution(opt.eventfile);
ie_ = find(strcmp({ev_.id}, eid_));
if isempty(ie_), error('A00b_make_station_line: event %s not in %s', eid_, opt.eventfile); end
ev_ = ev_(ie_);
[dc_, baz_] = gc_dist_az(lat0, lon0, ev_.lat, ev_.lon);   % centre -> event
az_ = mod(baz_ + 180, 360);                               % along the path, away from the event
% stations at signed distance s (km) from the centre along the great circle
% (s < 0 toward the event)
R_ = 6371;
ns_ = round(L_/ds_) + 1; s_ = linspace(-L_/2, L_/2, ns_)';
a_ = abs(s_)/R_; th_ = (az_ + 180*(s_ < 0))*pi/180;
p0_ = lat0*pi/180; l0_ = lon0*pi/180;
lat_ = asin(sin(p0_)*cos(a_) + cos(p0_)*sin(a_).*cos(th_));
lon_ = l0_ + atan2(sin(th_).*sin(a_)*cos(p0_), cos(a_) - sin(p0_)*sin(lat_));
lat_ = lat_*180/pi; lon_ = mod(lon_*180/pi + 180, 360) - 180;
fid_ = fopen(lfile_, 'w');
if fid_ < 0, error('A00b_make_station_line: cannot write %s', lfile_); end
for k_ = 1:ns_
    fprintf(fid_, '  %-5s %11.6f %13.6f %7.4f\n', sprintf('L%02d', k_), lat_(k_), lon_(k_), elev_);
end
fclose(fid_);
% the one-event catalogue (the event's block of opt.eventfile, verbatim)
txt_ = fileread(opt.eventfile);
blk_ = regexp(txt_, ['(?m)^[^\n#]*\n(?=event name:\s*' eid_ ')'], 'start');
ts_ = blk_(1); te_ = regexp(txt_(ts_:end), '\n\s*\n', 'once');
if isempty(te_), te_ = numel(txt_) - ts_ + 1; end
efile_ = fullfile(fileparts(lfile_), sprintf('events_%s.cmt', eid_));
fid_ = fopen(efile_, 'w');
[~, en_, ee_] = fileparts(opt.eventfile);
fprintf(fid_, '# %s from %s%s (Global CMT; single-event line test, A00b_make_station_line.m)\n\n', eid_, en_, ee_);
fprintf(fid_, '%s\n', regexprep(txt_(ts_:ts_ + te_ - 1), '\s+$', ''));
fclose(fid_);
% check: distances from the event, offset from the great circle
[de_, aze_] = gc_dist_az(ev_.lat*ones(ns_, 1), ev_.lon*ones(ns_, 1), lat_, lon_);
[~, aze0_] = gc_dist_az(ev_.lat, ev_.lon, lat0, lon0);
fprintf('%d stations every %.1f km on a %g km line through (%.4f, %.4f) (centre of %s)\n', ...
        ns_, L_/(ns_ - 1), L_, lat0, lon0, opt.line_ref);
fprintf('  along the great circle from %s (%.2f, %.2f; %.1f deg from the centre): distances %.2f-%.2f deg,\n', ...
        eid_, ev_.lat, ev_.lon, dc_, min(de_), max(de_));
fprintf('  azimuth from the event %.3f-%.3f deg (to the centre %.3f)\n', min(aze_), max(aze_), aze0_);
fprintf('  written to %s and %s\n', lfile_, efile_);
% map
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
fig = numbered_figure(4001, 'position', [50 50 800 650]); hold on; box on;
h1_ = plot(S_{3}, S_{2}, 'k^', 'markerfacecolor', [0.6 0.6 0.6], 'markersize', 7);
h2_ = plot(lon_, lat_, 'ro', 'markerfacecolor', 'r', 'markersize', 4);
plot(lon0, lat0, 'k+', 'markersize', 12, 'linewidth', 1.5);
text(lon_([1 end]), lat_([1 end]), {'  L01', sprintf('  L%02d', ns_)}, 'fontsize', 8);
daspect([1 cosd(lat0) 1]);
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
[~, fn_] = fileparts(lfile_);
legend([h1_ h2_], {sprintf('reference (%d)', numel(S_{2})), sprintf('line: %g km, %d stations', L_, ns_)}, ...
       'location', 'southoutside');
title(sprintf('%s (toward %s at azimuth %.0f deg)', fn_, eid_, baz_), 'interpreter', 'none');
od_ = fullfile(here_, 'figs_stations'); if exist(od_, 'dir') ~= 7, mkdir(od_); end
save_figure(fig, fullfile(od_, [fn_ '.png']), 130);
if fig_close_mode, close(fig); end
clear here_ eid_ L_ ds_ lfile_ fid_ S_ la_ lo_ v_ elev_ ev_ ie_ dc_ baz_ az_ R_ ns_ s_ a_ th_ p0_ l0_ lat_ lon_ k_ txt_ blk_ ts_ te_ efile_ en_ ee_ de_ aze_ aze0_ h1_ h2_ fn_ od_
