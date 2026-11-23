function parameters = aswms_parameters(opt, component)
% parameters = aswms_parameters(opt, component)
%
% The ASWMS-Q "parameters" structure (setup_parameters.m of
% github.com/jbrussell/ASWMS-Q) for this package, built from opt
% (parameters_RL.m), for one component ('BHT' Love, 'BHZ' Rayleigh, ...).
% The ASWMS-Q functions in functions/aswms call setup_parameters, which
% returns the structure installed with aswms_setup.m.
%
% Defaults are those of ASWMS-Q except where this package fixes them:
%   workingdir      opt.aswms_dir (default './ASWMS/'); SAC data in
%                   <workingdir>sacdata/<event id>/, event list evlist.txt
%   periods         15 log-spaced periods spanning the inversion band
%                   1/opt.fc(4) .. 1/opt.fc(1) (17-112 s by default)
%   lalim, lolim    station bounds of opt.stationfile +- 1 degree
%   refv            4 km/s (ASWMS-Q) for Rayleigh, 4.5 km/s for Love
%                   (component with 'T'): calc_SNR (A1) takes the signal
%                   between dist/(1.1 refv) and dist/(0.5 refv); with 4 km/s
%                   the Love waves (group velocity ~4.4 km/s) fall before it
%                   (synthetic test: Love SNR median 0.3, all stations
%                   rejected; with 4.5 km/s median 9.7, 99% kept)
% Any field can be overridden with opt.aswms (a struct of parameters
% fields).
if nargin < 1 || isempty(opt), opt = struct(); end
if nargin < 2 || isempty(component), component = 'BHZ'; end
og = @(f, d) optget(opt, f, d);
d2k = @(d) d*6371*pi/180;                       % deg2km (6371 km sphere)

parameters.workingdir = og('aswms_dir', ['.' filesep 'ASWMS' filesep]);
if parameters.workingdir(end) ~= '/' && parameters.workingdir(end) ~= filesep
    parameters.workingdir = [parameters.workingdir filesep];
end
parameters.proj_name = og('aswms_net', 'XX');
parameters.component = component;

% stations -> lalim, lolim
fid = fopen(og('stationfile', 'stations_nomelt.txt'));
S = textscan(fid, '%s %f %f %f'); fclose(fid);
parameters.lalim = [min(S{2}) - 1, max(S{2}) + 1];
parameters.lolim = [min(S{3}) - 1, max(S{3}) + 1];
parameters.gridsize = 0.5;

fc = og('fc', [1/112 1/100 1/20 1/17]);
parameters.periods = sort(logspace(log10(1/fc(4)), log10(1/fc(1)), 15));

% data selection (a1_a_sac2eventmat)
parameters.dbpath = [parameters.workingdir 'sacdata/'];
parameters.eventfile = 'evlist.txt';
parameters.minMw = 5.5;
parameters.maxdepth = 50;
parameters.snr_tol = 3;
parameters.resample_delta = og('dt', 4);       % (ASWMS-Q default 1 s; resamples only to a coarser step)

% auto_win_select (a3_a_run_autowinpick)
parameters.largest_epidist_range = 3000;
parameters.cycle_before = 2;
parameters.cycle_after = 5;
parameters.min_dist_tol = d2k(20);
parameters.max_dist_tol = d2k(115);
parameters.min_groupv = 2;
parameters.max_groupv = 5;
parameters.cent_freq = mean(1./parameters.periods);
parameters.min_sta_num = 5;
parameters.is_removebadwindows = 1;

% (used by calc_SNR and the window plots)
parameters.refv = 4;                            % (Rayleigh, ASWMS-Q default)
if any(component == 'T'), parameters.refv = 4.5; end   % Love: calc_SNR's signal window starts at 1.1*refv
parameters.refphv = ones(size(parameters.periods))*parameters.refv;
parameters.min_width = 0.10;
parameters.max_width = parameters.min_width;
parameters.prefilt_frac = 0.25;
parameters.prefilter = [1/(1/min(parameters.periods)*(1+parameters.prefilt_frac)), ...
                        1/(1/max(parameters.periods)*(1-parameters.prefilt_frac))];

% overrides
ov = og('aswms', struct());
fn = fieldnames(ov);
for k = 1:numel(fn), parameters.(fn{k}) = ov.(fn{k}); end
if isfield(ov, 'periods') && ~isfield(ov, 'cent_freq'), parameters.cent_freq = mean(1./parameters.periods); end
end
