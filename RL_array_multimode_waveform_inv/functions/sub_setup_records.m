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
% Data: idagrn6 seismograms of the true model with branches
% 0..opt.nbranch_data-1, 15 real 2012 Pacific-rim earthquakes
% (events_nomelt_2012.cmt, Global CMT), 16 NoMelt OBSs, transverse (Love)
% and vertical (Rayleigh) displacement, band-passed (fc), normalized by
% each record's peak, plus band-limited noise. Reference spectra Uref:
% idagrn6_mask single-branch seismograms (n = 0..opt.nm-1) of the
% reference model, processed the same way.
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
% Noise: band-limited (white noise through the processing band-pass),
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
mc.ximin = 0.90; mc.ximax = 1.20;            % U prior on xi
mc.p_xi  = 0.4;                              % fraction of moves on xi
mc.smooth_sig = og('smooth_sig', [0 0]);     % curvature prior on [Vsv Vsh] (km/s/km^2; 0 = off; lr_common)
mc.smooth_len = og('smooth_len', 10);        % ... one N(0, sig^2) curvature term per smooth_len km
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
mc.nuis_start = og('nuis_start', 'best');   % start of the nuisance fit: 'best' of old/full/band per record, or one of them (spec_refine)
mc.fcont_stages = og('fcont_stages', 4);     % stages of the frequency continuation ('full', 'band' and 'best' starts)
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
[Del_c, baz] = gc_dist_az(lat0, lon0, [ev.lat], [ev.lon]);   % array center -> event
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
% window per event picked from the data (auto_win_select.m, the ASWMS-Q
% picker: narrow-band envelope peaks tracked across the band, group delays
% fitted against distance, window from 2 cycles before to 5 cycles after
% the arrivals; options opt.autowin, a struct of auto_win_select options),
% Delta/v1 + t1 .. Delta/v2 + t2 (rec.winpara = [v1 t1 v2 t2]).
wmode_ = og('windows', 'single');
single_ = any(strcmp(wmode_, {'single', 'fund', 'auto'}));
autowin_ = og('autowin', struct());
U_start_ = og('U_start', []);
mc.single_window = single_;
segtap = tukeywin_(cfg.N, 0.1);          % taper of the time segments (data and Uref)
for wv = 1:2
    % reference group velocities -> window boundaries for this wave type
    U = P.(['grv_ref_' uname{wv}]);
    U0 = U(:,1); Uov = U(:, 2:end);
    win(wv).Ub = max(U0) + 0.2; win(wv).Uov1 = max(Uov(:)) + 0.6; win(wv).Uf2 = min(U0) - 0.35;
    if ~isempty(U_start_), win(wv).Uov1 = U_start_(min(wv, numel(U_start_))); end   % opt.U_start
    if strcmp(wmode_, 'fund')                                   % no overtone window
        fm_ = og('fund_margin', [0.4 0.2]);
        win(wv).Ub = max(U0) + fm_(min(wv, numel(fm_))); win(wv).Uov1 = win(wv).Ub;
    end
    fprintf('%s ref. group vel: fund %.2f-%.2f, overtones %.2f-%.2f km/s\n', ...
        wname{wv}, min(U0), max(U0), min(Uov(:)), max(Uov(:)));
    for e = 1:nev
        r = r + 1;
        rec(r).wave = wv; rec(r).event = e; rec(r).id = ev(e).id;
        rec(r).az_prop = mod(baz(e) + 180, 360);
        rec(r).Del0 = Del0(e);
        Dj = gc_dist_az(lat, lon, ev(e).lat*ones(Nst,1), ev(e).lon*ones(Nst,1)) * 6371*pi/180;
        % window start / overtone-fundamental boundary / end per station (s)
        tw1 = Dj/win(wv).Uov1; twb = Dj/win(wv).Ub; tw2 = Dj/win(wv).Uf2 + 60;
        if strcmp(wmode_, 'auto')                       % picked from the data
            rec(r).winpara = auto_win_select(P.(['u_true_' uname{wv}]){e}, P.t, Dj, 1./fc([4 1]), autowin_);
            tw1 = Dj/rec(r).winpara(1) + rec(r).winpara(2); twb = tw1;
            tw2 = Dj/rec(r).winpara(3) + rec(r).winpara(4);
        end
        if strcmp(wmode_, 'auto'), t0_ = min(tw1); else, t0_ = Del0(e)/win(wv).Uov1; end
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
        % data: full idagrn6 synthetics of the true model, same processing
        X = fft(P.(['u_true_' uname{wv}]){e}(it, :) .* segtap);
        Y = zeros(cfg.N, Nst); Y(cfg.ib, :) = X(cfg.ib, :) .* taper;
        clean = 2*real(ifft(Y));
        idx = []; wt = [];
        for j = 1:Nst
            io = find(t >= tw1(j) & t < twb(j));
            iF = find(t >= twb(j) & t < tw2(j));
            idx = [idx; io + (j-1)*cfg.N; iF + (j-1)*cfg.N];      %#ok<AGROW>
            if single_
                wt = [wt; (2*wv-1)*ones(numel(io) + numel(iF), 1)]; %#ok<AGROW>   % one surface-wave window
            else
                wt = [wt; (2*wv-1)*ones(size(io)); 2*wv*ones(size(iF))]; %#ok<AGROW>
            end
        end
        rec(r).idx = idx; rec(r).wt = wt;
        rec(r).first = [true; diff(idx) ~= 1 | diff(wt) ~= 0];   % segment starts
        % normalize each record (data and, in apply_reference.m, the
        % reference spectra) by its peak amplitude in the windows, as is
        % usual with real data
        rec(r).scale = max(abs(clean(idx)));
        rec(r).clean = clean(idx) / rec(r).scale;
        rec(r).Dj = Dj;
        % for the waveform figures: the step-1 starting model (all branches,
        % no nuisance correction)
        if isfield(P, ['u_start_' uname{wv}])
            X = fft(P.(['u_start_' uname{wv}]){e}(it, :) .* segtap);
            Y = zeros(cfg.N, Nst); Y(cfg.ib, :) = X(cfg.ib, :) .* taper;
            u = 2*real(ifft(Y)); rec(r).start = u(idx) / rec(r).scale;
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

%% ---- noise ----
% band-limited noise: white noise through the processing band-pass with a
% standard deviation of opt.noise (default 3%) of the peak amplitude of
% each (normalized) record
sig_true = og('noise', 0.03)*ones(1, 4);
bpf = zeros(cfg.N, 1); bpf(cfg.ib) = taper;
cc = zeros(4, 2);                         % lag-1 autocorrelation accumulators
for r = 1:nrec
    X = fft(randn(cfg.N, Nst));
    nz = 2*real(ifft(X .* bpf));          % one-sided filter -> real band-limited noise
    nz = nz / std(nz(:));
    nr = sig_true(2*rec(r).wave) * nz(rec(r).idx);
    rec(r).d = rec(r).clean + nr;
    nf = find(~rec(r).first);
    for w = unique(rec(r).wt)'
        i = nf(rec(r).wt(nf) == w);
        cc(w, :) = cc(w, :) + [sum(nr(i).*nr(i-1)), sum(nr(i-1).^2)];
    end
end
rho_emp = (cc(:,1) ./ cc(:,2))';
if single_, rho_emp([2 4]) = rho_emp([1 3]); end   % tied slots (single window)
fprintf('band-limited noise: empirical lag-1 correlation per window %s\n', mat2str(round(rho_emp*100)/100));
fprintf('%d stations, %d events, %d records, %d data samples\n', Nst, nev, nrec, ...
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
