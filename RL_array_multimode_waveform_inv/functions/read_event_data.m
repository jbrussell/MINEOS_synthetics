function D = read_event_data(opt, ev, sta_names, comp, tg)
% D = read_event_data(opt, ev, sta_names, comp, tg)
%
% Data of one event and component from the ASWMS-Q event file
% <opt.aswms_dir>eventmat/<ev.id>_<comp>.mat (A1-A3), on the time grid tg
% (s after the CMT centroid time ev.t_cent, i.e. the time axis of the
% idagrn6 reference synthetics), in the station order sta_names.
%
% ev        : one event of read_cmtsolution.m (id, t_cent)
% D         : [] if there is no event file or event.isgood == 0
%             (A3b_run_exam_winpara.m); otherwise
%   U       : [numel(tg) x Nst] traces; NaN where a trace does not reach
%             and NaN columns for stations that are missing or have
%             stadata.isgood < 0 (A3)
%   present : [Nst x 1] logical
%   dist    : [Nst x 1] epicentral distance of the event file (km, vdist)
%   winpara : [v1 t1 v2 t2] of A3 ([] if none)
%   toff    : event.otime - centroid time (s): the windows dist/v1 + t1 ..
%             dist/v2 + t2 are times after event.otime, i.e. tg = t + toff
%   nstd    : [Nst x 1] noise level before the first arrival: std of the
%             samples earlier than dist/13.5 km/s - 60 s after event.otime
%             (dist/13.5 km/s is before the P wave at 20-100 deg; any
%             samples before the origin are included), band-passed like
%             the data (opt.fc; Tukey taper, central 80% used); NaN if fewer
%             than opt.noise_min_n (40) such samples. opt.noise_pre_origin =
%             false uses only the samples after the origin (test option).
% Traces with another sample interval than tg are low-passed (2nd-order
% Butterworth at the new Nyquist frequency, zero phase, as in A1) and
% interpolated (spline).
wd = optget(opt, 'aswms_dir', ['ASWMS' filesep]);
if wd(end) ~= '/' && wd(end) ~= filesep, wd = [wd filesep]; end
f = [wd 'eventmat' filesep ev.id '_' comp '.mat'];
D = [];
if exist(f, 'file') ~= 2, return; end
S = load(f); event = S.event;
if isfield(event, 'isgood') && ~event.isgood, return; end
tg = tg(:); dt = tg(2) - tg(1); Nst = numel(sta_names);
tc = (ev.t_cent - datenum(1970, 1, 1))*86400;    % centroid time (epoch s)
D.U = NaN(numel(tg), Nst); D.present = false(Nst, 1); D.dist = NaN(Nst, 1); D.nstd = NaN(Nst, 1);
fc = optget(opt, 'fc', [1/112 1/100 1/20 1/17]);
D.winpara = []; if isfield(event, 'winpara') && numel(event.winpara) == 4, D.winpara = event.winpara; end
D.toff = event.otime - tc;
for k = 1:numel(event.stadata)
    s = event.stadata(k);
    j = find(strcmp(strtrim(s.stnm), sta_names));
    if isempty(j), continue; end
    if isfield(s, 'isgood') && ~isempty(s.isgood) && s.isgood < 0, continue; end
    x = s.data(:); tk = (s.otime - tc) + (0:numel(x)-1)'*s.delta;
    in = tg >= tk(1) - 1e-6 & tg <= tk(end) + 1e-6;   % part of tg the trace covers (NaN elsewhere)
    if ~any(in), continue; end
    if abs(s.delta - dt) < 1e-9 && abs(rem((tg(1) - tk(1))/dt, 1)) < 1e-6
        ik = round((tg(in) - tk(1))/dt) + 1;      % same grid: no interpolation
        D.U(in, j) = x(ik);
    else
        if s.delta < dt                          % anti-alias (as in A1)
            [b, a] = butter(2, s.delta/dt, 'low');
            x = filtfilt(b, a, x);
        end
        D.U(in, j) = interp1(tk, x, tg(in), 'spline');
    end
    D.present(j) = true; D.dist(j) = s.dist;
    to = tk - D.toff;                            % time after event.otime (s)
    ip = find(to < s.dist/13.5 - 60);            % before the first arrival (P)
    if ~optget(opt, 'noise_pre_origin', true), ip = ip(to(ip) >= 0); end
    if numel(ip) >= optget(opt, 'noise_min_n', 40)
        D.nstd(j) = noise_std_(s.data(ip), s.delta, fc);
    end
end
end

function sd = noise_std_(x, dt, fc)
% std of x band-passed as the data (Tukey taper, one-sided band-pass taper),
% central 80% of the samples
x = x(:) - mean(x); n = numel(x);
f = (0:n-1)'/(n*dt); ib = find(f > fc(1) & f < fc(4));
X = fft(x .* tukeywin_(n, 0.1)); Y = zeros(n, 1); Y(ib) = X(ib) .* bp_taper(f(ib), fc);
u = 2*real(ifft(Y)); k = round(0.1*n) + 1 : round(0.9*n);
sd = std(u(k));
end
