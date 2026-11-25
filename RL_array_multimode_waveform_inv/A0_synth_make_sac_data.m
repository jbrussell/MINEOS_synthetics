% A0_SYNTH_MAKE_SAC_DATA  Synthetic test only (pre-step): make the "data" as
% SAC files, so that they go through the same data preparation as real data
% (A1-A3, ASWMS-Q).
%
% For the true model (CARDS/<opt.trueid>.card, PA5) MINEOS mode tables of
% branches 0..opt.nbranch_data-1 and idagrn6 full mode-sum seismograms
% (transverse and vertical displacement) at the stations of
% opt.stationfile for the events of opt.eventfile (Global CMT), opt.dt
% sampling, opt.synth_len_hr hours from the centroid time, preceded by
% opt.synth_pre_s s of pre-event noise (default 0: the records start at
% the origin; the noise levels are measured before the P wave). Noise is added: white noise
% through the processing band-pass (opt.fc), with a standard deviation of
% opt.noise times the peak band-passed amplitude of the event and wave type,
% times a factor per station 10^(opt.synth_noise_spread*U[-1,1]) (default
% spread 0: the same at all stations).
%
% Output, in opt.aswms_dir (ASWMS-Q layout):
%   sacdata/<event id>/<net>.<station>.<component>.sac
%       components opt.aswms_components = {Love, Rayleigh} ('BHT', 'BHZ');
%       reference time (NZ*) = origin = centroid time, O = 0, B = -pre_s;
%       EVLA, EVLO, EVDP (m),
%       MAG (Mw), STLA, STLO, STEL (m), KSTNM, KCMPNM, KNETWK, DELTA
%   sacdata/evlist.txt  (event ids)
%   <opt.truthfile>     (synthetic test only, never read for real data):
%       noise-free records and the noise, per wave type and event
%       ([Nt x Nsta], station order of opt.stationfile), the true model id,
%       branch count, noise level and times
% Takes about 10 min (MINEOS + idagrn6).
%
% Next: A1_sac2eventmat.m, A2_cleanup_events.m, A3_run_autowinpick.m
% (optional A3b_run_exam_winpara.m), then D2_step2_prep_mineos_idagrn.m.

parameters_RL;
wd = opt.aswms_dir; if wd(end) ~= '/' && wd(end) ~= filesep, wd = [wd filesep]; end
sacdir = [wd 'sacdata' filesep];
if ~exist(sacdir, 'dir'), mkdir(sacdir); end
addpath(fullfile(opt.mdir, 'functions'));

% stations and events
fid = fopen(opt.stationfile); S = textscan(fid, '%s %f %f %f'); fclose(fid);
sta = struct('name', {S{1}}, 'lat', S{2}, 'lon', S{3}, 'elev', S{4});
ev = read_cmtsolution(opt.eventfile);
nev = numel(ev); Nst = numel(sta.name);
% optional arrival-angle anomalies (opt.synth_arrival_shift = [min max], deg;
% [] = none): each source is moved along the small circle about the array
% centre by dth_e = +-U[min, max] (random sign; seed opt.synth_seed + 1), so the
% epicentral distance to the array centre is unchanged and the wave arrives
% with propagation azimuth baz + 180 + dth_e, i.e. an arrival-angle anomaly
% dth_e (the convention of the inversion, event_G_spec.m) relative to the
% catalogue location, which the SAC headers keep. Geometry in geocentric
% latitudes as idagrn6 and the inversion. The true dth_e are saved in the
% truth file (truth.dth_true, deg; truth.ev_shifted).
ev_src = ev; dth_true = zeros(1, nev);
ash = optget(opt, 'synth_arrival_shift', []);
if ~isempty(ash)
    rs_ = rng; rng(opt.synth_seed + 1);
    dth_true = sign(rand(1, nev) - 0.5) .* (ash(1) + (ash(2) - ash(1))*rand(1, nev));
    rng(rs_);
    gl0 = geoc_lat(mean(sta.lat)); lo0 = mean(sta.lon);
    [Dc, bz] = gc_dist_az(gl0, lo0, geoc_lat([ev.lat]), [ev.lon]);
    for e = 1:nev
        a = (bz(e) + dth_true(e))*pi/180; d = Dc(e)*pi/180; p0 = gl0*pi/180;
        p2 = asin(sin(p0)*cos(d) + cos(p0)*sin(d)*cos(a));
        l2 = lo0 + atan2(sin(a)*sin(d)*cos(p0), cos(d) - sin(p0)*sin(p2))*180/pi;
        % geocentric -> geographic latitude: tan(lat_g) = tan(lat_c)/(1-f)^2,
        % f = 1/298.25 as geoc_lat.m (idagrn6's series; they agree to 2e-6 deg)
        lg = atan(tan(p2)/(1 - 1/298.25)^2)*180/pi;
        ev_src(e).lat = lg; ev_src(e).lon = mod(l2 + 180, 360) - 180;
    end
    fprintf('arrival-angle anomalies (deg, by event): %s\n', mat2str(round(dth_true*100)/100));
end
Tdisp = optget(opt, 'Tdisp', logspace(log10(0.95*17), log10(1.02*112), 10));

% MINEOS + idagrn6 of the true model (all data branches)
U = cell(2, nev); TYs = 'TS';
tic;
for wv = 1:2
    TY = TYs(wv);
    mineos_lock(opt.mdir, 'acquire');
    try
        mineos_modes(opt.mdir, opt.trueid, TY, opt.nbranch_data, opt.maxF, Tdisp, []);
        D = idagrn_synth(opt.mdir, opt.trueid, TY, opt.maxF, opt.nbranch_data, ev_src, sta, ...
                         opt.synth_len_hr, opt.dt, []);
    catch err
        mineos_lock(opt.mdir, 'release'); rethrow(err);
    end
    mineos_lock(opt.mdir, 'release');
    U(wv, :) = D.u(:, 1)';
    t = D.t;
    fprintf('idagrn6 %c (true model %s, branches 0-%d): %.0f s\n', TY, opt.trueid, opt.nbranch_data - 1, toc);
end

% noise: white noise through the processing band-pass, std = opt.noise x
% the peak band-passed amplitude of the event and wave type, times a factor
% per station (and event, wave type) 10^(opt.synth_noise_spread*U[-1,1])
% (default spread 0: the same level at all stations); the records start
% opt.synth_pre_s s before the origin with noise only (pre-event noise, for
% the per-record noise levels of the inversion)
rng(opt.synth_seed);
fc = optget(opt, 'fc', [1/112 1/100 1/20 1/17]);
Nt = numel(t); f = (0:Nt-1)'/(Nt*opt.dt);
bpf = zeros(Nt, 1); ib = find(f > fc(1) & f < fc(4)); bpf(ib) = bp_taper(f(ib), fc);
bp = @(X) 2*real(ifft(fft(X) .* bpf));             % one-sided band-pass (real output)
Npre = round(optget(opt, 'synth_pre_s', 0)/opt.dt); Ne = Npre + Nt;
fe = (0:Ne-1)'/(Ne*opt.dt);
bpe = zeros(Ne, 1); ibe = find(fe > fc(1) & fe < fc(4)); bpe(ibe) = bp_taper(fe(ibe), fc);
spread = optget(opt, 'synth_noise_spread', 0);
noise = cell(2, nev); noise_pre = noise; nfac = noise;
for wv = 1:2
    for e = 1:nev
        A = max(max(abs(bp(U{wv, e}))));
        nz = 2*real(ifft(fft(randn(Ne, Nst)) .* bpe)); nz = nz / std(nz(:));
        nfac{wv, e} = 10.^(spread*(2*rand(1, Nst) - 1));
        nz = opt.noise * A * nz .* nfac{wv, e};
        noise_pre{wv, e} = nz(1:Npre, :); noise{wv, e} = nz(Npre+1:end, :);
    end
end

% SAC files and event list
comps = opt.aswms_components;
fid = fopen([sacdir 'evlist.txt'], 'w');
for e = 1:nev
    fprintf(fid, '%s\n', ev(e).id);
    edir = [sacdir ev(e).id filesep];
    if ~exist(edir, 'dir'), mkdir(edir); end
    t0 = ev(e).t_cent;                                % SAC reference time = origin (centroid)
    for wv = 1:2
        for j = 1:Nst
            H = struct('DELTA', opt.dt, 'KSTNM', sta.name{j}, 'KCMPNM', comps{wv}, ...
                       'KNETWK', opt.aswms_net, 'STLA', sta.lat(j), 'STLO', sta.lon(j), ...
                       'STEL', 1000*sta.elev(j), 'EVLA', ev(e).lat, 'EVLO', ev(e).lon, ...
                       'EVDP', 1000*ev(e).dep, 'MAG', ev(e).Mw);
            fn = sprintf('%s%s.%s.%s.sac', edir, opt.aswms_net, sta.name{j}, comps{wv});
            x = [noise_pre{wv, e}(:, j); U{wv, e}(:, j) + noise{wv, e}(:, j)];
            mksac(fn, x, t0, H);
            sac_set_times(fn, t(1) - Npre*opt.dt, opt.dt, numel(x));   % B < 0: pre-event noise
        end
    end
end
fclose(fid);

% truth (synthetic test only)
truth = struct('trueid', opt.trueid, 'nbranch_data', opt.nbranch_data, 'noise', opt.noise, ...
               'seed', opt.synth_seed, 'dt', opt.dt, 't', t, 'comps', {comps}, ...
               'ev_id', {{ev.id}}, 't_cent', [ev.t_cent], 'sta', {sta.name}, ...
               'clean', {U}, 'noise_tr', {noise}, 'pre_s', Npre*opt.dt, ...
               'noise_spread', spread, 'noise_fac', {nfac});
truth.dth_true = dth_true; truth.ev_shifted = struct('lat', [ev_src.lat], 'lon', [ev_src.lon]);
save([wd opt.truthfile], 'truth', '-v7');
fprintf('%d events x %d stations x 2 components written to %s; truth in %s%s\n', ...
    nev, Nst, sacdir, wd, opt.truthfile);
