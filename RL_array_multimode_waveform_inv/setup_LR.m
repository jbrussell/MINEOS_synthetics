% SETUP_LR  Configuration, NoMelt geometry, real events, MINEOS dispersion
% and idagrn6 synthetic data for the joint radially anisotropic Love +
% Rayleigh multimode array waveform inversion (recovery test of PA5 with a
% PREM reference model).
%
% Model: Vsv(z) and xi(z) = (Vsh/Vsv)^2 are independent transdimensional
% cubic B-spline profiles between the Moho (7 km below the seafloor) and
% 350 km (no discontinuities), with interior control points kx <= kv. The
% bottom end values (350 km) are tied to the reference model.
% Vsh = Vsv*sqrt(xi). Vp, density, eta, the water, the crust and the
% structure below 350 km are held at the reference model.
%
% Reference ("input") model: PREM below the Moho, with the NoMelt (PA5)
% water, sediments and crust (mineos_tools/make_ref_card.m).
% True model: PA5 (Gaherty et al., 1996), pa5_5km.card.
%
% Dispersion: MINEOS. Phase velocities of branches n = 0..2 are the MINEOS
% reference values plus the MINEOS Vsv and Vsh kernels times the model
% perturbation (mineos_dispersion.m); Rayleigh waves are radially
% anisotropic (they see Vsh through the MINEOS kernels).
%
% Data: idagrn6 synthetics of the true model (all mode branches in the
% table, opt.nbranch_data) for 15 real 2012 Pacific-rim earthquakes
% (events_nomelt_2012.cmt, Global CMT) at the 16 NoMelt OBSs, transverse
% (Love) and vertical (Rayleigh) displacement, band-passed (fc), plus
% band-limited noise. Reference spectra Uref: idagrn6_mask single-branch
% synthetics (n = 0..2) of the reference model, processed the same way.
% Everything from MINEOS and idagrn6 is computed once by
% prep_mineos_idagrn.m and cached in opt.prepfile.
%
% Path correction: the source-to-array path is PA5 in the data but PREM in
% Uref, which delays the fundamental Rayleigh wave by up to ~100 s at 16 s
% (many cycles). Before the inversion, the array-average path phase and
% amplitude of each record (fundamental and overtone windows) are measured
% from the data and folded into Uref (path_precorrect.m), so the nuisance
% model only has to absorb the residual. The in-array moveout, which
% carries the information on the local structure, is not changed.
%
% Nuisance model for each record and mode (event_G_spec): a smooth complex
% excitation/path correction B_n(w) on opt.spec_K hat functions in log w
% (marginalized), a time shift tau_n and a frequency-dependent arrival
% angle dth_n(w) on opt.spec_Ka Legendre terms in normalized log w.
%
% Noise: band-limited (white noise through the processing band-pass),
% modelled with an exponential (AR(1)) covariance with a level sigma and
% correlation r for each of the 4 window types (Love/Rayleigh x
% overtone/fundamental), all sampled hierarchically. Fixed seed.
%
% Optional overrides: define a struct OPT before calling.

if ~exist('opt', 'var'), opt = struct(); end
og = @(f, d) optget(opt, f, d);
rng(3);
outdir = og('outdir', 'figs');
if ~exist(outdir, 'dir'), mkdir(outdir); end
addpath(fullfile(fileparts(mfilename('fullpath')), 'mineos_tools'));

%% ---- MCMC settings ----
mc.nit   = og('nit', 8000);  mc.burn = og('burn', 3000);  mc.thin = 5;
mc.kvmin = 0;  mc.kvmax = og('kvmax', 15);   % interior Vsv control points
mc.kxmin = 0;  mc.kxmax = og('kxmax', 15);   % interior xi control points
mc.vsmin = 4.0; mc.vsmax = 5.0;              % U prior on Vsv (km/s)
mc.ximin = 0.90; mc.ximax = 1.20;            % U prior on xi
mc.p_xi  = 0.4;                              % fraction of moves on xi
mc.step_vs = 0.03; mc.step_xi = 0.01; mc.step_z = 10;
mc.step_tau = 0.4; mc.step_dth = 0.25; mc.step_lsig = 0.05; mc.step_rho = 0.02;
mc.rhomax = 0.98;                            % U[0, rhomax] prior on correlation
mc.tau_prior = og('tau_prior', 15); mc.dth_prior = 5; mc.sb = og('sb', 0.4);
% spectral nuisance model (event_G_spec)
mc.spec_K  = og('spec_K', 3);                % hat functions for B_n(w)
mc.spec_Ka = og('spec_Ka', 2);               % Legendre terms for dth_n(w)
mc.dthp = [mc.dth_prior, og('dth1_prior', 5)*ones(1, mc.spec_Ka - 1)];  % prior std (deg)
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
CrefL = P.phv_ref_T; CrefR = P.phv_ref_S;     % reference (PREM) dispersion
CtrueL = P.phv_true_T; CtrueR = P.phv_true_S; % true (PA5) dispersion, MINEOS
% true and reference profiles (below the seafloor) for plotting
zplot = (0:2:400)';
card_prof = @(c, f) card_at_depth(c, f, zplot + cfg.zw)/1000;
tru.vsv = card_prof(P.card_true, 'vsv'); tru.vsh = card_prof(P.card_true, 'vsh');
tru.xi = (tru.vsh ./ tru.vsv).^2;
refp.vsv = card_prof(P.card_ref, 'vsv'); refp.vsh = card_prof(P.card_ref, 'vsh');
refp.xi = (refp.vsh ./ refp.vsv).^2;
cfg.refp = struct('z', zplot, 'vsv', refp.vsv, 'vsh', refp.vsh);
% linearized MINEOS dispersion of the true model's Vsv, Vsh between zmin
% and zmax (the rest of the model at the reference): how well the
% parameterization and the linearization can reproduce the true dispersion
mt.vsv = card_at_depth(P.card_true, 'vsv', cfg.lin.z + cfg.zw)/1000;   % mantle side at the Moho
mt.vsh = card_at_depth(P.card_true, 'vsh', cfg.lin.z + cfg.zw)/1000;
dm_ = [mt.vsv - cfg.lin.vsv0; mt.vsh - cfg.lin.vsh0];
ClinL = CrefL; ClinR = CrefR;
for n = 1:nm, ClinL(:,n) = ClinL(:,n) + cfg.lin.GL{n}*dm_; ClinR(:,n) = ClinR(:,n) + cfg.lin.GR{n}*dm_; end
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
cfg.dt = P.dt; cfg.N = 512;
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
segtap = tukeywin_(cfg.N, 0.1);          % taper of the time segments (data and Uref)
for wv = 1:2
    % reference group velocities -> window boundaries for this wave type
    U = P.(['grv_ref_' uname{wv}]);
    U0 = U(:,1); Uov = U(:, 2:end);
    win(wv).Ub = max(U0) + 0.2; win(wv).Uov1 = max(Uov(:)) + 0.6; win(wv).Uf2 = min(U0) - 0.35;
    fprintf('%s ref. group vel: fund %.2f-%.2f, overtones %.2f-%.2f km/s\n', ...
        wname{wv}, min(U0), max(U0), min(Uov(:)), max(Uov(:)));
    for e = 1:nev
        r = r + 1;
        rec(r).wave = wv; rec(r).event = e; rec(r).id = ev(e).id;
        rec(r).az_prop = mod(baz(e) + 180, 360);
        rec(r).Del0 = Del0(e);
        rec(r).tred = cfg.dt*round((Del0(e)/win(wv).Uov1 - 150)/cfg.dt);
        Dj = gc_dist_az(lat, lon, ev(e).lat*ones(Nst,1), ev(e).lon*ones(Nst,1)) * 6371*pi/180;
        t  = (0:cfg.N-1)'*cfg.dt + rec(r).tred;
        it = round((t - P.t(1))/cfg.dt) + 1;          % samples of the idagrn traces
        if it(end) > numel(P.t), error('idagrn records too short: increase opt.len_hr'); end
        % reference spectra (single-branch idagrn6 synthetics of the reference model)
        Uref = cell(1, nm);
        for n = 1:nm
            X = fft(P.(['u_ref_' uname{wv}]){e, n}(it, :) .* segtap);
            Uref{n} = X(cfg.ib, :) .* taper;
        end
        rec(r).Uref = Uref; rec(r).cref = Cref_f{wv};
        % data: full idagrn6 synthetics of the true model, same processing
        X = fft(P.(['u_true_' uname{wv}]){e}(it, :) .* segtap);
        Y = zeros(cfg.N, Nst); Y(cfg.ib, :) = X(cfg.ib, :) .* taper;
        clean = 2*real(ifft(Y));
        idx = []; wt = [];
        for j = 1:Nst
            io = find(t >= Dj(j)/win(wv).Uov1 & t < Dj(j)/win(wv).Ub);
            iF = find(t >= Dj(j)/win(wv).Ub   & t < Dj(j)/win(wv).Uf2 + 60);
            idx = [idx; io + (j-1)*cfg.N; iF + (j-1)*cfg.N];      %#ok<AGROW>
            wt  = [wt; (2*wv-1)*ones(size(io)); 2*wv*ones(size(iF))]; %#ok<AGROW>
        end
        rec(r).idx = idx; rec(r).wt = wt;
        rec(r).first = [true; diff(idx) ~= 1 | diff(wt) ~= 0];   % segment starts
        % normalize each record (data and reference spectra) by its peak
        % amplitude in the windows, as is usual with real data
        rec(r).scale = max(abs(clean(idx)));
        rec(r).clean = clean(idx) / rec(r).scale;
        for n = 1:nm, rec(r).Uref{n} = rec(r).Uref{n} / rec(r).scale; end
        rec(r).Dj = Dj;
    end
end
nrec = numel(rec);
wave = [rec.wave]';
[cfg.Phi, cfg.Psi] = spec_basis(cfg.w, mc.spec_K, mc.spec_Ka);

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
fprintf('band-limited noise: empirical lag-1 correlation per window %s\n', mat2str(round(rho_emp*100)/100));
fprintf('%d stations, %d events, %d records, %d data samples\n', Nst, nev, nrec, ...
    sum(arrayfun(@(x) numel(x.d), rec)));

%% ---- path correction of the reference spectra (from the data) ----

if og('path_correct', true)
    [rec, pc] = path_precorrect(rec, cfg, nm, og('pc', struct()));
    % compare with the true path delay (array center, MINEOS PA5 - reference)
    for r = 1:nrec
        Ct = Ctrue_f{rec(r).wave}; Cr = Cref_f{rec(r).wave};
        pc.true_delay{r} = rec(r).Del0 * (1./Ct - 1./Cr);     % [Nf x nm] (s)
        pc.err(r) = sqrt(mean((pc.tau{r}(taper > 0.5, 1) - pc.true_delay{r}(taper > 0.5, 1)).^2));
    end
    pc.dropped = {};
    for r = find(~pc.ok)
        pc.dropped{end+1} = sprintf('%s %s', wname{rec(r).wave}, rec(r).id);
    end
    fprintf('path correction (Rayleigh mode 0): rms error vs true path delay %s s\n', ...
        mat2str(round(pc.err(wave == 2 & pc.ok')*10)/10));
    if ~isempty(pc.dropped)
        fprintf('dropped (phase unwrapping failed, low coherence): %s\n', strjoin(pc.dropped, ', '));
    end
    keep = find(pc.ok);
    for f = {'tau', 'tau_meas', 'coh', 'amp', 'true_delay'}, pc.(f{1}) = pc.(f{1})(keep); end
    pc.err = pc.err(keep);
    rec = rec(keep); nrec = numel(rec); wave = [rec.wave]';
else
    pc = [];
end
% There is no "true" nuisance model for idagrn data, so the log-likelihood
% of the true structure is evaluated with nuisance parameters fitted to it
% (fit_true_LR.m, called by the drivers); NaN until then.
Ltrue = NaN;
tru.tau = []; tru.dth = []; tru.dth1 = [];
