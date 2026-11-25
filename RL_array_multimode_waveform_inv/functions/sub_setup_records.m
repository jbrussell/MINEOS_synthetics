% SUB_SETUP_RECORDS  (drivee script, called by D3_step2_waveform_inv.m and
% D4_step2_plot_results.m) Configuration, NoMelt geometry, real events,
% MINEOS dispersion and idagrn6 data for the step-2 multimode array
% waveform inversion. Needs opt from parameters_RL.m and the MINEOS/idagrn6
% results cached by D2_step2_prep_mineos_idagrn.m (opt.prepfile).
%
% Model: Vsv(z) and xi(z) = (Vsh/Vsv)^2 are independent transdimensional
% cubic B-spline profiles between the Moho (7 km below the seafloor) and
% 350 km (no discontinuities), with interior control points kx <= kv. The
% bottom end values (350 km) are tied to the reference model.
% Vsh = Vsv*sqrt(xi). Vpv, Vph and density follow Vsv and Vsh with the
% Vpv/Vsv, Vph/Vsh and rho/Vsv ratios of the reference model; eta, the
% water, the crust and the structure below 350 km are held at the
% reference model.
%
% Reference model: the isotropic step-1 model (D1), which has a lid and a
% low-velocity zone. True model: PA5 (Gaherty et al., 1996).
%
% Dispersion: MINEOS. Phase velocities of branches n = 0..opt.nm-1 are the
% MINEOS reference values plus the MINEOS kernels times the model
% perturbation (mineos_dispersion.m, mineos_lin_setup.m); Rayleigh waves
% are radially anisotropic (they see Vsh through the MINEOS kernels).
%
% Data (opt.data_source = 'aswms', default): the ASWMS-Q event files
% <opt.aswms_dir>eventmat/<event id>_<component>.mat made by A1-A3 from
% SAC files (read_event_data.m), as for real data; for the synthetic test
% the SAC files are idagrn6 seismograms of the true model with branches
% 0..opt.nbranch_data-1 plus band-limited noise (A0_synth_make_sac_data.m),
% for 15 real 2012 Pacific-rim earthquakes (events_nomelt_2012.cmt, Global
% CMT) at 16 NoMelt OBSs, transverse (Love) and vertical (Rayleigh)
% displacement. Events removed by A1-A3 and stations without data are
% left out. Each record is band-passed (fc) and normalized by its peak in
% the windows. ('internal': the true-model synthetics of the prep file,
% with band-limited noise added here.) Reference spectra Uref: idagrn6_mask
% single-branch seismograms (n = 0..opt.nm-1) of the reference model,
% processed the same way.
%
% Windows: from the group velocities of the reference model (MINEOS):
% overtone window from Delta/(max U of branches 1..nm-1 + 0.6) to
% Delta/(max U0 + 0.2), fundamental-mode window from there to
% Delta/(min U0 - 0.35) + 60 s, for each station's epicentral distance.
%
% Nuisance model for each record and mode (event_G_spec): a smooth complex
% excitation/path correction B_n(w) on opt.spec_K hat functions in log w,
% a frequency-dependent time shift tau_n(w) on opt.spec_Kt Legendre terms
% and a frequency-dependent arrival angle dth_n(w) on opt.spec_Ka Legendre
% terms in normalized log w. By default (opt.nuisance = 'marginal') all
% three are integrated out (B exactly, tau and dth by linearization about
% points that are updated every opt.relin iterations; nuis_relin.m).
%
% Noise: in the data ('aswms'; the synthetic level and correlation are
% measured against the noise-free records of A0) or band-limited white
% noise added here ('internal'),
% modelled with an exponential (AR(1)) covariance with a level sigma and
% correlation r for each of the 4 window types (Love/Rayleigh x
% overtone/fundamental), all sampled hierarchically. Fixed seed.

if ~exist('opt', 'var'), opt = struct(); end
og = @(f, d) optget(opt, f, d);
rng(3);
outdir = og('outdir', 'figs');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ---- MCMC settings ----
mc.nit   = og('nit', 8000);  mc.burn = og('burn', 3000);  mc.thin = 5;
mc.kvmin = 0;  mc.kvmax = og('kvmax', 15);   % interior Vsv control points
mc.kxmin = 0;  mc.kxmax = og('kxmax', 15);   % interior xi control points
mc.vsmin = 4.0; mc.vsmax = 5.0;              % U prior on Vsv (km/s)
mc.ximin = og('ximin', 0.85); mc.ximax = og('ximax', 1.15);   % U prior on xi (control points; default symmetric about 1)
mc.p_xi  = 0.4;                              % fraction of moves on xi
mc.smooth_sig = og('smooth_sig', [0 0]);     % curvature prior on the smooth_param profiles (0 = off; lr_common)
mc.smooth_len = og('smooth_len', 10);        % ... one N(0, sig^2) curvature term per smooth_len km
mc.smooth_param = og('smooth_param', 'vsv_vsh'); % ... profiles smoothed: 'vsv_vsh' or 'vsv_xi'
mc.vsv_ref_sig = og('vsv_ref_sig', 0);       % prior toward the reference Vsv (km/s, or 'step1'; 0 = off; lr_common)
mc.vsv_ref_len = og('vsv_ref_len', 10);      % ... one N(Vsv_ref, sig^2) term per vsv_ref_len km
mc.xi_ref_sig = og('xi_ref_sig', 0);         % Gaussian prior on xi toward xi_ref_mean (0 = off; lr_common)
mc.xi_ref_mean = og('xi_ref_mean', 1);       % ... its mean (1 = isotropic)
mc.xi_ref_len = og('xi_ref_len', 10);        % ... one N(xi_ref_mean, sig^2) term per xi_ref_len km
mc.xi_ref_hier = og('xi_ref_hier', false);   % ... strength sampled (hierarchical) instead of xi_ref_sig
mc.xi_ref_sig_range = og('xi_ref_sig_range', [0.005 0.2]);   % ... log-uniform hyperprior on it
mc.smooth_hier = og('smooth_hier', false);   % curvature prior strengths sampled (hierarchical) instead of smooth_sig
mc.smooth_sig_range = og('smooth_sig_range', [5e-5 5e-3]);   % ... log-uniform hyperprior on them
% Vsv-only inversion (opt.xi_fixed, step 1 waveform version): xi held at
% the reference model's (cfg.xi_fixed below; no xi moves, no xi control points)
if og('xi_fixed', false), mc.p_xi = 0; mc.kxmax = 0; end
% eta (opt.eta_free): held at the reference model by default; if free, a
% third transdimensional B-spline profile (ke <= kv, bottom tied to the
% reference, U prior), updated in a fraction p_eta of the moves
mc.eta_free = og('eta_free', false);
mc.kemin = 0; mc.kemax = og('kemax', 15);
mc.etamin = og('etamin', 0.80); mc.etamax = og('etamax', 1.05);
mc.p_eta = og('p_eta', 0.2); mc.step_eta = og('step_eta', 0.01);
mc.step_vs = 0.03; mc.step_xi = 0.01; mc.step_z = 10;
mc.step_tau = 0.4; mc.step_dth = 0.25; mc.step_lsig = 0.05; mc.step_rho = 0.02;
mc.rhomax = 0.98;                            % U[0, rhomax] prior on correlation
mc.tau_prior = og('tau_prior', 15); mc.dth_prior = 5; mc.sb = og('sb', 0.4);
% spectral nuisance model (event_G_spec)
mc.spec_K  = og('spec_K', 3);                % hat functions for B_n(w)
mc.spec_Ka = og('spec_Ka', 2);               % Legendre terms for dth_n(w)
mc.spec_Kt = og('spec_Kt', 3);               % Legendre terms for tau_n(w)
mc.dthp = [mc.dth_prior, og('dth1_prior', 5)*ones(1, mc.spec_Ka - 1)];   % prior std (deg)
mc.taup = [mc.tau_prior, og('dtau_prior', 10)*ones(1, mc.spec_Kt - 1)];  % prior std (s)
mc.relin = og('relin', 1);                   % relinearize tau, dth, B every relin iterations
mc.nuis_start = og('nuis_start', 'best');   % start of the nuisance fit: 'best' per record of 'direct', 'continuation_global',
                                             %   'continuation_banded', or one of them (spec_refine)
mc.fcont_stages = og('fcont_stages', 4);     % stages of the frequency continuation
mc.jump_every = og('jump_every', 100);       % cycle-skip check of tau_n(w) every jump_every burn-in iterations (0 = off; cycle_jump)
mc.jump_tol = og('jump_tol', 0);             % log-likelihood gain required to keep a cycle change
mc.max_dtau = 2; mc.max_dth = 3;             % max step of a relinearization (s, deg)
mc.nchains = og('nchains', 6);

%% ---- MINEOS and idagrn6 (cached) ----
P = prep_mineos_idagrn(opt);
nm = P.nm;
Tdisp = P.Tdisp(:)';
fc = og('fc', [1/112 1/100 1/20 1/17]);       % processing band (Hz)
cfg.fc = fc;

%% ---- model space ----
cfg.zw = max(P.card_true.z(P.card_true.vsv == 0 & P.card_true.z < 20));   % water depth (km)
cfg.moho = 7; cfg.zmin = 7; cfg.zmax = 350;
cfg.lin = mineos_lin_setup(P.card_ref, P.K_ref_T, P.K_ref_S, P.phv_ref_T, P.phv_ref_S, ...
                           cfg.zw, cfg.zmin, cfg.zmax);
cfg.need = isfinite(P.phv_ref_T) & isfinite(P.phv_ref_S);   % branch exists at the node
cfg.Tdisp = Tdisp;
% bottom end values tied to the reference model at zmax
[~, ib_] = max(cfg.lin.z);
cfg.vbot = cfg.lin.vsv0(ib_); cfg.xbot = (cfg.lin.vsh0(ib_)/cfg.lin.vsv0(ib_))^2;
cfg.xi_fixed = og('xi_fixed', false);       % xi = the reference model's (Vsv-only; lr_common.m, mineos_dispersion.m)
mc.vsv_ref_sigz = vsv_ref_sigz(opt, cfg);   % depth-dependent Vsv prior width (opt.vsv_ref_sig = 'step1'), else []
cfg.ebot = cfg.lin.eta0(ib_);
CrefL = P.phv_ref_T; CrefR = P.phv_ref_S;     % reference (step-1 model) dispersion
CtrueL = P.phv_true_T; CtrueR = P.phv_true_S; % true (PA5) dispersion, MINEOS
% true and reference profiles (below the seafloor) for plotting
zplot = (0:2:400)';
card_prof = @(c, f) card_at_depth(c, f, zplot + cfg.zw)/1000;
tru.vsv = card_prof(P.card_true, 'vsv'); tru.vsh = card_prof(P.card_true, 'vsh');
tru.xi = (tru.vsh ./ tru.vsv).^2;
refp.vsv = card_prof(P.card_ref, 'vsv'); refp.vsh = card_prof(P.card_ref, 'vsh');
refp.xi = (refp.vsh ./ refp.vsv).^2;
tru.eta = card_at_depth(P.card_true, 'eta', zplot + cfg.zw);      % dimensionless
refp.eta = card_at_depth(P.card_ref, 'eta', zplot + cfg.zw);
cfg.refp = struct('z', zplot, 'vsv', refp.vsv, 'vsh', refp.vsh, 'eta', refp.eta);
% linearized MINEOS dispersion of the true model's Vsv, Vsh between zmin
% and zmax (the rest of the model at the reference): how well the
% parameterization and the linearization can reproduce the true dispersion
mt.vsv = card_at_depth(P.card_true, 'vsv', cfg.lin.z + cfg.zw)/1000;   % mantle side at the Moho
mt.vsh = card_at_depth(P.card_true, 'vsh', cfg.lin.z + cfg.zw)/1000;
[ClinL, ClinR] = mineos_lin_disp(cfg.lin, mt.vsv, mt.vsh, [], 'both');   % eta at the reference
fprintf('true - reference phase velocity (rms, km/s) Love %s Rayleigh %s\n', ...
    mat2str(round(sqrt(mean((CtrueL - CrefL).^2))*1000)/1000), mat2str(round(sqrt(mean((CtrueR - CrefR).^2))*1000)/1000));
fprintf('linearized PA5 - true PA5 (rms, km/s)      Love %s Rayleigh %s\n', ...
    mat2str(round(sqrt(mean((ClinL - CtrueL).^2))*1000)/1000), mat2str(round(sqrt(mean((ClinR - CtrueR).^2))*1000)/1000));

%% ---- stations ----
sta = P.sta.name; lat = P.sta.lat; lon = P.sta.lon;
lat0 = mean(lat); lon0 = mean(lon);
cfg.xy = og('aperture', 1) * [(lon - lon0)*111.195*cosd(lat0), (lat - lat0)*111.195];
cfg.exact_geom = og('exact_geom', true);     % moveout from great-circle distances (event_G_spec.m)
Nst = numel(sta);

%% ---- time / frequency ----
cfg.dt = P.dt; cfg.N = og('N', 800);
f = (0:cfg.N-1)'/(cfg.N*cfg.dt);
cfg.ib = find(f > fc(1) & f < fc(4));
fb = f(cfg.ib); cfg.w = 2*pi*fb;
taper = bp_taper(fb, fc);
disp_band = @(C) disp_to_band(C, Tdisp, fb);
Cref_f = {disp_band(CrefL), disp_band(CrefR)};
Ctrue_f = {disp_band(CtrueL), disp_band(CtrueR)};

%% ---- events and records ----
ev = P.ev; nev = numel(ev);
% distances and azimuths on the sphere from GEOCENTRIC latitudes (geoc_lat.m),
% as in idagrn6 (reference spectra, synthetic data)
[Del_c, baz] = gc_dist_az(geoc_lat(lat0), lon0, geoc_lat([ev.lat]), [ev.lon]);   % array center -> event
Del0 = Del_c * 6371*pi/180;                                   % km
wname = {'Love', 'Rayleigh'}; uname = {'T', 'S'};
rec = struct([]); r = 0;
% windows (opt.windows): 'single' (default) = one surface-wave window per
% wave type from Delta/U_start to Delta/U_f + 60 s (overtones and the
% fundamental mode; one noise level and correlation per wave type, labelled
% 2*wv-1, the slots 2*wv-1 and 2*wv tied in the samplers); 'split' = an
% overtone window (to Delta/U_b) and a fundamental-mode window, each with
% its own noise parameters; 'fund' = the fundamental-mode window only
% (from Delta/U_w, U_w = max fundamental-mode group velocity of the
% reference + opt.fund_margin(wave) km/s, default [0.4 0.2] for [Love
% Rayleigh], to Delta/U_f + 60 s): cuts the early, fast higher modes (also
% the branches >= nm that are not modelled). opt.U_start (km/s; scalar or [Love
% Rayleigh]): start of the window ('single'); default (empty) max group
% velocity of branches 1..nm-1 of the reference + 0.6 km/s. 'auto' = one
% window per event picked from the data by ASWMS-Q (A3_run_autowinpick.m,
% auto_win_select.m: narrow-band envelope peaks tracked across the band,
% group delays fitted against distance, window from 2 cycles before to 5
% cycles after the arrivals): dist/v1 + t1 .. dist/v2 + t2 after the event
% origin of the event file (rec.winpara = event.winpara = [v1 t1 v2 t2]).
wmode_ = og('windows', 'single');
% data source (opt.data_source): 'aswms' = the ASWMS-Q event files of A1-A3
% (<opt.aswms_dir>eventmat/; read_event_data.m), with the noise already in
% the data; for the synthetic test, the noise-free records of A0
% (opt.truthfile) give the noise level; 'internal' = true-model synthetics
% of the prep file with noise added here
aswms_ = strcmp(og('data_source', 'internal'), 'aswms');
nlev_mode_ = og('noise_levels', 'pre'); nlev_n_ = 0;   % per-record noise levels ('aswms' data)
comps_ = og('aswms_components', {'BHT', 'BHZ'});
synth_ = false;
if aswms_
    wd_ = og('aswms_dir', ['ASWMS' filesep]); if wd_(end) ~= '/' && wd_(end) ~= filesep, wd_ = [wd_ filesep]; end
    tf_ = [wd_ og('truthfile', 'synth_truth.mat')];
    synth_ = exist(tf_, 'file') == 2;
    if synth_
        truth = load(tf_); truth = truth.truth;
        it_tr = round((P.t(1) - truth.t(1))/cfg.dt);      % A0 time axis -> idagrn time axis
    end
end
single_ = any(strcmp(wmode_, {'single', 'fund', 'auto'}));
U_start_ = og('U_start', []);
mc.single_window = single_;
segtap = tukeywin_(cfg.N, 0.1);          % taper of the time segments (data and Uref)
waves_ = og('waves', [1 2]);             % wave types used: 1 = Love, 2 = Rayleigh ([2]: Rayleigh only, step 1)
for wv = waves_(:)'
    % reference group velocities -> window boundaries for this wave type
    U = P.(['grv_ref_' uname{wv}]);
    U0 = U(:,1); Uov = U(:, 2:end);
    if isempty(Uov), Uov = U0; end                   % (fundamental mode only: no overtones)
    win(wv).Ub = max(U0) + 0.2; win(wv).Uov1 = max(Uov(:)) + 0.6; win(wv).Uf2 = min(U0) - 0.35;
    if ~isempty(U_start_), win(wv).Uov1 = U_start_(min(wv, numel(U_start_))); end   % opt.U_start
    if strcmp(wmode_, 'fund')                                   % no overtone window
        fm_ = og('fund_margin', [0.4 0.2]);
        win(wv).Ub = max(U0) + fm_(min(wv, numel(fm_))); win(wv).Uov1 = win(wv).Ub;
    end
    fprintf('%s ref. group vel: fund %.2f-%.2f, overtones %.2f-%.2f km/s\n', ...
        wname{wv}, min(U0), max(U0), min(Uov(:)), max(Uov(:)));
    for e = 1:nev
        Dj = gc_dist_az(geoc_lat(lat), lon, geoc_lat(ev(e).lat)*ones(Nst,1), ev(e).lon*ones(Nst,1)) * 6371*pi/180;
        if aswms_                                   % data of the ASWMS-Q event file
            Dd = read_event_data(opt, ev(e), sta, comps_{wv}, P.t);
            if isempty(Dd)
                fprintf('%s %s: no event file (or event.isgood = 0): not used\n', wname{wv}, ev(e).id); continue
            end
            if strcmp(wmode_, 'auto') && isempty(Dd.winpara)
                fprintf('%s %s: no window from A3 (winpara): not used\n', wname{wv}, ev(e).id); continue
            end
        end
        r = r + 1;
        rec(r).wave = wv; rec(r).event = e; rec(r).id = ev(e).id;
        rec(r).az_prop = mod(baz(e) + 180, 360);
        rec(r).Del0 = Del0(e);
        % window start / overtone-fundamental boundary / end per station (s)
        tw1 = Dj/win(wv).Uov1; twb = Dj/win(wv).Ub; tw2 = Dj/win(wv).Uf2 + 60;
        if strcmp(wmode_, 'auto')                       % picked from the data (A3, ASWMS-Q)
            if aswms_
                wp = Dd.winpara; Dw = Dd.dist; toff = Dd.toff;  % times after event.otime -> after the centroid
            else
                error('opt.windows = ''auto'' needs opt.data_source = ''aswms'' (windows picked by A3_run_autowinpick.m)');
            end
            rec(r).winpara = wp;
            tw1 = Dw/wp(1) + wp(2) + toff; twb = tw1;
            tw2 = Dw/wp(3) + wp(4) + toff;
        end
        if strcmp(wmode_, 'auto'), t0_ = min(tw1(isfinite(tw1))); else, t0_ = Del0(e)/win(wv).Uov1; end
        rec(r).tred = cfg.dt*round((t0_ - 150)/cfg.dt);
        % (a late window start, opt.U_start, must not push the segment past
        % the end of the idagrn records)
        rec(r).tred = min(rec(r).tred, P.t(1) + cfg.dt*floor((P.t(end) - P.t(1))/cfg.dt - (cfg.N - 1)));
        t  = (0:cfg.N-1)'*cfg.dt + rec(r).tred;
        it = round((t - P.t(1))/cfg.dt) + 1;          % samples of the idagrn traces
        if it(end) > numel(P.t) || it(1) < 1, error('idagrn records too short: increase opt.len_hr'); end
        if max(tw2) > t(end) || min(tw1) < t(1)
            warning('record %d (%s): windows cut by the %d-sample segment (opt.N, opt.len_hr)', r, ev(e).id, cfg.N);
        end
        proc = @(U) proc_segment(U, segtap, taper, cfg);   % segment taper + band-pass (as Uref)
        if aswms_
            Ud = Dd.U(it, :);
            pres = Dd.present(:)' & all(isfinite(Ud), 1);   % stations with data over the whole segment
            Ud(:, ~pres) = 0;
            dat = proc(Ud);
            if synth_                                   % noise-free records of A0 (diagnostics, noise level)
                ke = find(strcmp(truth.ev_id, ev(e).id));
                Uc = truth.clean{wv, ke}(it + it_tr, :); Uc(:, ~pres) = 0;
                clean = proc(Uc);
            else
                clean = NaN(size(dat));
            end
        else
            % data: full idagrn6 synthetics of the true model, same processing
            clean = proc(P.(['u_true_' uname{wv}]){e}(it, :));
            pres = true(1, Nst); rec(r).nacf_ = []; rec(r).npsd_ = []; rec(r).nstd_ = [];
        end
        idx = []; wt = [];
        for j = find(pres)
            io = find(t >= tw1(j) & t < twb(j));
            iF = find(t >= twb(j) & t < tw2(j));
            idx = [idx; io + (j-1)*cfg.N; iF + (j-1)*cfg.N];      %#ok<AGROW>
            if single_
                wt = [wt; (2*wv-1)*ones(numel(io) + numel(iF), 1)]; %#ok<AGROW>   % one surface-wave window
            else
                wt = [wt; (2*wv-1)*ones(size(io)); 2*wv*ones(size(iF))]; %#ok<AGROW>
            end
        end
        rec(r).idx = idx; rec(r).wt = wt; rec(r).present = pres(:);
        rec(r).first = [true; diff(idx) ~= 1 | diff(wt) ~= 0];   % segment starts
        % normalize each record (data and, in apply_reference.m, the
        % reference spectra) by its peak amplitude in the windows, as is
        % usual with real data
        if aswms_
            rec(r).scale = max(abs(dat(idx)));
            rec(r).d = dat(idx) / rec(r).scale;
            rec(r).clean = clean(idx) / rec(r).scale;
            % per-record, per-station noise levels from the noise before the
            % first arrival (opt.noise_levels = 'pre'; read_event_data.m): noise std of a sample = sigma x the
            % measured level of its station / the record's peak
            rec(r).nlev = []; rec(r).nacf_ = Dd.nacf;    % (pre-P noise autocorrelation per station)
            rec(r).npsd_ = Dd.npsd; rec(r).nstd_ = Dd.nstd(:)';   % (pre-P noise spectra and levels)
            if strcmp(nlev_mode_, 'pre')
                lev = Dd.nstd(:)' / rec(r).scale;          % 1 x Nst
                okl = pres & isfinite(lev) & lev > 0;
                if any(okl)
                    lev(~okl) = median(lev(okl));          % (stations without a measurement)
                    jst = ceil(idx / cfg.N);               % station of each sample
                    rec(r).nlev = lev(jst)';
                    nlev_n_ = nlev_n_ + 1;
                else
                    warning('record %d (%s): no pre-event noise measured; no per-station noise levels', r, ev(e).id);
                end
            end
        else
            rec(r).scale = max(abs(clean(idx)));
            rec(r).clean = clean(idx) / rec(r).scale;
        end
        rec(r).Dj = Dj;
        % for the waveform figures: the step-1 starting model (all branches,
        % no nuisance correction)
        if isfield(P, ['u_start_' uname{wv}])
            u = proc(P.(['u_start_' uname{wv}]){e}(it, :));
            rec(r).start = u(idx) / rec(r).scale;
        else
            rec(r).start = NaN(size(idx));
        end
    end
end
% reference spectra Uref (single-branch idagrn6 synthetics of the
% reference model, processed like the data), reference phase velocities and
% the MINEOS linearization (apply_reference.m)
[rec, cfg] = apply_reference(rec, cfg, P, Tdisp, fb);
nrec = numel(rec);
wave = [rec.wave]';
[cfg.Phi, cfg.Psi, cfg.Psit] = spec_basis(cfg.w, mc.spec_K, mc.spec_Ka, mc.spec_Kt);

%% ---- noise correlation model (opt.noise_model) ----
% 'ar1' (default): AR(1) per window type with the correlation r sampled;
% 'ar2': AR(2) per wave type fitted to the pooled pre-P noise
% autocorrelation (read_event_data.m), fixed; 'arp': the same with order
% opt.noise_arp; 'emp': each station of each
% record its own AR(p) (opt.noise_arp, default 8) from its own pre-P noise
% (the pooled fit where a station's noise is too short), fixed. The fits
% (levinson_ar.m, white nugget opt.noise_nugget = 0.05) give each record a
% sparse whitening operator rec.W with log det rec.wld (whiten_operator.m);
% only sigma is then sampled.
nmodel_ = og('noise_model', 'ar1');
if any(strcmp(nmodel_, {'ar2', 'arp', 'emp'}))
    if ~aswms_, error('opt.noise_model = ''%s'' needs opt.data_source = ''aswms'' (pre-P noise)', nmodel_); end
    if strcmp(nmodel_, 'ar2'), p_ = 2; else, p_ = og('noise_arp', 8); end
    nug_ = og('noise_nugget', 0.05);
    pooled_ = NaN(2, size(rec(1).nacf_, 2));
    for wv = 1:2
        A_ = vertcat(rec([rec.wave] == wv).nacf_); A_ = A_(isfinite(A_(:, 1)), :);
        for L_ = 1:size(A_, 2), a_ = A_(isfinite(A_(:, L_)), L_); if ~isempty(a_), pooled_(wv, L_) = mean(a_); end, end
        [ph_, v_] = levinson_ar(pooled_(wv, :), p_, nug_);
        Mp_(wv) = struct('phi', {ph_}, 'v', v_);                    %#ok<SAGROW>
        fprintf('%s noise, pooled pre-P autocorrelation (lags 0-%d) %s; AR(%d) coefficients %s, innovation var %.3f\n', ...
            wname{wv}, p_, mat2str(round(pooled_(wv, 1:p_+1)*100)/100), numel(v_) - 1, mat2str(round(ph_{end}*1000)/1000), v_(end));
    end
    nown_ = 0; ntot_ = 0;
    for r = 1:nrec
        jst_ = ceil(rec(r).idx / cfg.N);
        if any(strcmp(nmodel_, {'ar2', 'arp'}))     % one model per wave type (pooled)
            [rec(r).W, rec(r).wld] = whiten_operator(rec(r).first, ones(size(jst_)), Mp_(rec(r).wave));
        else
            M_ = repmat(Mp_(rec(r).wave), 1, Nst);
            for j = unique(jst_)'
                a_ = rec(r).nacf_(j, :); ntot_ = ntot_ + 1;
                if all(isfinite(a_(1:p_+1)))
                    [ph_, v_] = levinson_ar(a_, p_, nug_); M_(j) = struct('phi', {ph_}, 'v', v_); nown_ = nown_ + 1;
                end
            end
            [rec(r).W, rec(r).wld] = whiten_operator(rec(r).first, jst_, M_);
        end
    end
    if strcmp(nmodel_, 'emp')
        fprintf('empirical AR(%d) noise models: %d of %d record-stations from their own pre-P noise (others pooled)\n', p_, nown_, ntot_);
    end
    mc.fixed_corr = true; mc.noise_model = nmodel_;
    rho_emp_model_ = pooled_(:, 2)';
end

% 'fdom': frequency-domain likelihood (fdom_operator.m): each station's
% window tapered (Tukey 0.2) and transformed to its in-band DFT
% coefficients, scaled by the noise power spectrum of the station from its
% pre-P noise (read_event_data.m): opt.noise_psd = 'shape' (default; the
% pooled spectral shape of the wave type x the station's own level, robust
% with short pre-P noise) or 'station' (each station's own smoothed
% spectrum; for real data with station-specific spectra). sigma is a scale
% factor (~1); the per-station levels are in the spectra (rec.nlev unused).
if strcmp(nmodel_, 'fdom')
    if ~aswms_, error('opt.noise_model = ''fdom'' needs opt.data_source = ''aswms'' (pre-P noise)'); end
    fg_ = linspace(fc(1), fc(4), 60); psdm_ = og('noise_psd', 'shape');
    shape_ = NaN(2, numel(fg_));
    for wv = 1:2
        A_ = [];
        for r = find([rec.wave] == wv)
            A_ = [A_; rec(r).npsd_ ./ (rec(r).nstd_(:).^2)];          %#ok<AGROW>
        end
        for k_ = 1:numel(fg_), a_ = A_(isfinite(A_(:, k_)), k_); if ~isempty(a_), shape_(wv, k_) = median(a_); end, end
    end
    nown_ = 0; ntot_ = 0;
    for r = 1:nrec
        jst_ = ceil(rec(r).idx / cfg.N);
        lv_ = rec(r).nstd_(:).^2; lv_(~isfinite(lv_)) = median(lv_(isfinite(lv_)));
        Srow_ = lv_ * shape_(rec(r).wave, :);                             % [Nst x 60]
        if strcmp(psdm_, 'station')
            for j = unique(jst_)'
                ntot_ = ntot_ + 1;
                if all(isfinite(rec(r).npsd_(j, :))), Srow_(j, :) = rec(r).npsd_(j, :); nown_ = nown_ + 1; end
            end
        end
        Srow_ = Srow_ / rec(r).scale^2;                                   % (records are normalized by their peak)
        Srow_ = max(Srow_, 1e-3*max(Srow_, [], 2));                       % (floor: band edges)
        Sfun_ = @(j, f) interp1(fg_, Srow_(j, :), f, 'linear', 'extrap');
        rec(r).W = fdom_operator(rec(r).first, rec(r).wt, jst_, cfg.dt, [fc(1) fc(4)], Sfun_);
        rec(r).wld = rec(r).W.ld; rec(r).nlev = [];
    end
    fprintf('frequency-domain noise model (%s spectra): %d coefficients per record on average (%d samples)\n', psdm_, ...
        round(mean(arrayfun(@(x) numel(x.W.wt), rec))), round(mean(arrayfun(@(x) numel(x.d), rec))));
    if strcmp(psdm_, 'station'), fprintf('   %d of %d record-stations with their own spectrum\n', nown_, ntot_); end
    mc.fixed_corr = true; mc.noise_model = 'fdom';
end

%% ---- noise ----
% band-limited noise: white noise through the processing band-pass with a
% standard deviation of opt.noise (default 3%) of the peak amplitude of
% each (normalized) record
% ('aswms': the noise is in the data; its level and correlation are
% measured from data - noise-free records for the synthetic test, and set
% to opt.noise_guess and 0.5 for real data, as starting values only: both
% are sampled)
sig_true = og('noise', 0.03)*ones(1, 4);
bpf = zeros(cfg.N, 1); bpf(cfg.ib) = taper;
cc = zeros(4, 2);                         % lag-1 autocorrelation accumulators
ss = zeros(4, 2);                         % noise power accumulators ('aswms')
ssw_ = zeros(4, 2);                       % (whitened, fixed noise models)
for r = 1:nrec
    if aswms_
        nr = rec(r).d - rec(r).clean;     % (NaN for real data)
        if isfield(rec, 'nlev') && ~isempty(rec(r).nlev), nr = nr ./ rec(r).nlev; end   % relative to the measured levels
        if isfield(rec, 'W') && ~isempty(rec(r).W)    % fixed noise model: whitened units (sigma ~ 1)
            Wr_ = rec(r).W;
            if isstruct(Wr_), yw_ = Wr_.M*nr; wtw_ = Wr_.wt; else, yw_ = Wr_*nr; wtw_ = rec(r).wt; end
            for w = unique(wtw_)', ssw_(w, :) = ssw_(w, :) + [sum(yw_(wtw_ == w).^2), sum(wtw_ == w)]; end
        end
    else
        X = fft(randn(cfg.N, Nst));
        nz = 2*real(ifft(X .* bpf));      % one-sided filter -> real band-limited noise
        nz = nz / std(nz(:));
        nr = sig_true(2*rec(r).wave) * nz(rec(r).idx);
        rec(r).d = rec(r).clean + nr;
    end
    nf = find(~rec(r).first);
    for w = unique(rec(r).wt)'
        i = nf(rec(r).wt(nf) == w);
        cc(w, :) = cc(w, :) + [sum(nr(i).*nr(i-1)), sum(nr(i-1).^2)];
        ss(w, :) = ss(w, :) + [sum(nr(rec(r).wt == w).^2), sum(rec(r).wt == w)];
    end
end
rho_emp = (cc(:,1) ./ cc(:,2))';
if aswms_
    if synth_
        sig_true = sqrt(ss(:,1) ./ ss(:,2))';
        if any(ssw_(:, 2) > 0), sig_true = sqrt(ssw_(:,1) ./ max(ssw_(:,2), 1))'; end   % (whitened units)
        if single_, sig_true([2 4]) = sig_true([1 3]); end
        fprintf('noise of the A0 data (std; with per-station levels: relative to them) per window %s\n', mat2str(round(sig_true*1e4)/1e4));
    elseif nlev_n_ > 0
        sig_true = ones(1, 4); rho_emp = 0.5*ones(1, 4);      % (scale on the measured levels)
    else
        sig_true = og('noise_guess', 0.03)*ones(1, 4); rho_emp = 0.5*ones(1, 4);
    end
    if nlev_n_ > 0
        fprintf('per-station noise levels from the pre-event noise in %d of %d records (sigma = scale on them)\n', nlev_n_, nrec);
    end
end
if single_, rho_emp([2 4]) = rho_emp([1 3]); end   % tied slots (single window)
fprintf('band-limited noise: empirical lag-1 correlation per window %s\n', mat2str(round(rho_emp*100)/100));
fprintf('%d stations, %d events, %d records, %d data samples\n', Nst, numel(unique([rec.event])), nrec, ...
    sum(arrayfun(@(x) numel(x.d), rec)));

%% ---- path delay of the data relative to the reference (diagnostic) ----
% MINEOS true - reference phase velocities over the distance to the array
% center: the part of tau_n(w) the nuisance model has to absorb
pd.true_delay = cell(1, nrec);
for r = 1:nrec
    pd.true_delay{r} = rec(r).Del0 * (1./Ctrue_f{rec(r).wave} - 1./Cref_f{rec(r).wave});   % [Nf x nm] (s)
end
% There is no "true" nuisance model for idagrn data, so the log-likelihood
% of the true structure is evaluated with nuisance parameters fitted to it
% (fit_nuisance_LR.m, called by D3); NaN until then.
Ltrue = NaN;
% synthetic test only (opt.fix_tau): time shifts held at the true path
% delays, exactly at every frequency (apply_fixed_tau.m)
if og('fix_tau', false), [rec, mc] = apply_fixed_tau(rec, cfg, pd, mc); end
% synthetic test only (opt.fix_nuisance): B_n(w), tau_n(w) and dth_n(w) held
% at B = 1, tau = 0, dth = 0 (prior std 1e-3) - their true values when the
% reference model of the prep file is the TRUE model and the data have the
% modelled branches only (oracle nuisance test; set opt.jump_every = 0)
if og('fix_nuisance', false)
    mc.sb = 1e-3; mc.taup = 1e-3*ones(size(mc.taup)); mc.dthp = 1e-3*ones(size(mc.dthp));
    mc.fix_nuisance = true;
    fprintf('nuisance parameters held at B = 1, tau = 0, dth = 0 (opt.fix_nuisance)\n');
end
