% PARAMETERS_RL  Settings for the two-step inversion. Called by every
% driver (D1_..., D2_..., D3_..., D4_...). Edit the values here, or define
% a struct OPT with any of these fields before running a driver.
%
% Step 1 (D1): fundamental-mode Rayleigh phase velocities (MINEOS, 20-100 s)
%   -> isotropic Vsv model (xi = 1), linearized inversion as in
%   ../run_MINEOS/c2_RLphv_inv_vsvh_vpvh_rho_eta.m. True model PA5, starting
%   model an isotropic PA5 without its lid/LVZ (Vs linear in depth from the
%   Moho to 250 km) with the Vp/Vs, rho/Vs, eta and Q of PA5
%   (opt.s1_start_model = 'pa5_iso'), or isotropic PREM ('prem').
% Step 2 (D2-D4): multimode Love + Rayleigh array waveform inversion
%   (rjMcMC) for Vsv, Vsh and xi, with the step-1 model as the reference.
%   Data: idagrn6 seismograms of PA5 with branches 0..nbranch_data-1;
%   modeled: branches 0..nm-1.
% In both steps Vpv, Vph and density follow Vsv and Vsh with the
% Vpv/Vsv, Vph/Vsh and rho/Vsv ratios of the starting/reference model.

if ~exist('opt', 'var'), opt = struct(); end
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, 'functions'));
addpath(fullfile(here_, 'functions', 'aswms'));                 % ASWMS-Q functions (A1-A3)
addpath(fullfile(here_, 'functions', 'octave_shims'), '-end');  % only where a toolbox lacks them
mdir_ = fullfile(here_, '..', 'run_MINEOS');
addpath(fullfile(mdir_, 'functions'));
opt = opt_defaults(opt, struct( ...
    ... % ---- project: name and archive of a run ----
    'project_tag', '', ...               % NAME OF THE PROJECT (run): when set, the run tag is '_<project_tag>'
    ...                                  %   (setup_<project_tag>.mat, chain_<project_tag>_XX.mat,
    ...                                  %   figs_<project_tag>/) and D4 archives the run in
    ...                                  %   <results_root>/<project_tag>/ (save_project_results.m: figures,
    ...                                  %   summaries, mat files, waveform data, this file and the exact
    ...                                  %   options; rerun with opt = load_project('<project_tag>')).
    ...                                  %   '' = from the run tag (opt.tag, default '_<sampler>') without the '_'
    'results_root', 'results', ...       % folder of the project archives
    'save_results', true, ...            % archive the run at the end of D4
    'results_data', true, ...            % also archive the waveform data (opt.aswms_dir), the reference prep
    ...                                  %   file and the model cards (~50-150 MB per project)
    ... % ---- models and files ----
    'mdir',      mdir_, ...              % run_MINEOS (MINEOS cards, drivers, ../FORTRAN/bin)
    'trueid',    'pa5_5km', ...          % true model (CARDS/<trueid>.card)
    'premid',    'prem_35', ...          % PREM card
    'premref',   'prem_nomelt', ...      % PREM below the true Moho, true water/crust above (make_ref_card)
    'zmoho',     12, ...                 % Moho of the true model (km below sea level)
    'eventfile', fullfile(here_, 'events_nomelt_2012.cmt'), ...
    'stationfile', fullfile(here_, 'stations_nomelt.txt'), ...
    ... % ---- step 1: isotropic Rayleigh inversion (D1) ----
    's1_periods',   logspace(log10(20), log10(150), 24), ...   % (s) step-1 dispersion periods (24, log-spaced, 20-150 s)
    's1_cstd',      0.0025, ...          % data uncertainty (fraction of c); also the std of the synthetic noise
    's1b_noise',    true, ...            % synthetic step-1 dispersion data: add Gaussian noise of std s1_cstd*c
    's1b_noise_seed', 1, ...             %   (D1 linearized and D1_step1_Bayesian_disp.m; seed of the noise)
    's1_start_model', 'pa5_iso', ...     % starting model of step 1: 'pa5_iso' = isotropic PA5 without lid/LVZ
    ...                                  %   (Vs linear in depth from the Moho to s1_zlin km, PA5 Vp/Vs, rho/Vs,
    ...                                  %   eta and Q; functions/pa5_iso_ref_card.m); 'prem' = isotropic PREM
    's1_zlin',      250, ...             % (km below sea level) 'pa5_iso': PA5 below this depth
    's1_start',     '', ...              % starting card (written to CARDS; '' = RLinv_pa5iso or RLinv_prem_iso)
    's1_work',      'RLinv_s1_work', ... % card used for the MINEOS recalculations
    's1_final',     'RLinv_step1', ...   % step-1 model (written to CARDS): the step-2 reference
    's1_file',      'step1_result.mat', ...
    's1_eps',       struct('data', 1, 'H', 0.02, 'J', 0.001, 'F', 2, 'vpvs', 10, ...
                           'rhovs', 10, 'xi', 100, 'eta', 10), ...
    's1_z_dampbot', 300, 's1_zmax', 400, ...   % (km below sea level)
    's1_nit',       4, ...               % iterations (MINEOS recalculated after each)
    ... % ---- data: ASWMS-Q event files (A0 synthetic only; A1-A3 for synthetic and real data) ----
    'data_source',  'aswms', ...         % 'aswms': the data are read from the ASWMS-Q event files
    ...                                  %   <aswms_dir>eventmat/<event id>_<component>.mat (A1-A3; for the
    ...                                  %   synthetic test written from the SAC files of A0); 'internal': the
    ...                                  %   previous behaviour (true-model synthetics computed in step 2, noise
    ...                                  %   added to the processed records)
    'aswms_dir',    ['ASWMS' filesep], ... % ASWMS-Q working directory (sacdata/, eventmat/, winpara/, branchmat/)
    'aswms_components', {{'BHT', 'BHZ'}}, ... % channel names of the [Love Rayleigh] records (SAC KCMPNM)
    'aswms_net',    'XX', ...            % network code of the synthetic SAC files
    'aswms',        struct(), ...        % overrides of the ASWMS-Q parameters (aswms_parameters.m), e.g.
    ...                                  %   struct('snr_tol', 3, 'maxdepth', 50, 'min_sta_num', 5)
    'synth_len_hr', 2, ...               % A0: length of the synthetic SAC records (h; ASWMS-Q needs dist/2 km/s)
    'synth_seed',   3, ...               % A0: random seed of the noise
    'truthfile',    'synth_truth.mat', ... % A0: noise-free records, noise and true model (in aswms_dir)
    'noise_guess',  0.03, ...            % real data without noise levels: starting noise level (fraction of peak)
    'noise_levels', 'pre', ...           % 'pre': noise std of each record and station measured before the first
    ...                                  %   arrival (dist/13.5 km/s - 60 s after the origin, plus any samples
    ...                                  %   before the origin); the sampled sigma is a scale factor on it;
    ...                                  %   'none': one sigma per wave type
    ... % -- noise model: how the correlation between neighbouring noise samples enters the likelihood --
    ...   %   noise_model picks ONE model; the three settings below it apply only to some models
    ...   %   (marked [used by: ...]) and are IGNORED otherwise (with the default 'ar1': all three).
    'noise_model',  'ar1', ...           % 'ar1'  (DEFAULT) AR(1): sample i correlates with sample i-1 through r;
    ...                                  %        r and sigma are SAMPLED in the inversion (the "r = [...]" in the logs)
    ...                                  % 'fdom' frequency-domain likelihood: the in-band spectrum of each station's
    ...                                  %        window weighted by its noise spectrum from the pre-P noise
    ...                                  %        (fdom_operator.m); meant for real data. FIXED (r not sampled)
    ...                                  % EXPERIMENTAL, AR models FITTED ONCE to the pre-P noise and then FIXED:
    ...                                  % 'ar2'  AR(2), one per wave type, pooled over stations
    ...                                  % 'arp'  AR(p), p = noise_arp, one per wave type, pooled over stations
    ...                                  % 'emp'  AR(p), p = noise_arp, one per station of each record
    ...                                  %   (test with band-passed noise, 4 s sampling: none of the fixed AR models
    ...                                  %   whitens it better than AR(1) - residual lag-1 correlation 0.6-0.65 vs
    ...                                  %   0.56 - because the band is ~1/3 of the Nyquist range, so the noise
    ...                                  %   covariance is nearly singular; high orders are unstable)
    'noise_psd',    'shape', ...         % [used by: 'fdom'] noise spectra: 'shape' = one spectral shape pooled over
    ...                                  %   stations x each station's level; 'station' = each station's own spectrum
    'noise_arp',    8, ...               % [used by: 'arp', 'emp'] AR order p: each noise sample predicted from the
    ...                                  %   p samples before it ('ar2' always uses p = 2)
    'noise_nugget', 0.05, ...            % [used by: 'ar2', 'arp', 'emp'] white-noise fraction added to the pre-P
    ...                                  %   noise autocorrelation before the AR fit (keeps the fit stable)
    'synth_pre_s',  0, ...               % A0: noise before the origin (s; 0: records start at the origin)
    'synth_noise_spread', 0, ...         % A0: noise level per station x 10^(spread*U[-1,1]) (0: all equal)
    'synth_arrival_shift', [], ...       % A0: arrival-angle anomalies: sources moved along the small circle about the
    ...                                  %   array centre by +-U[min max] deg (e.g. [2 5]); [] = none (truth.dth_true)
    ... % ---- step 2: MINEOS + idagrn6 (D2) ----
    'prepfile',     'prep_step2.mat', ...
    'refid',        'RLinv_step1', ...   % reference card of step 2 (= s1_final)
    'ref_qmod',     '', ...              % Q (.qmod) of the step-1 and step-2 cards ('' = trueid for
    ...                                  %   s1_start_model = 'pa5_iso', premid for 'prem')
    'nbranch_data', 4, ...               % branches 0..3 in the PA5 data (synthetic test)
    'nm',           4, ...               % branches 0..3 modeled
    'maxF',         70, ...              % mode tables up to (mHz)
    'dt',           4, 'len_hr', 1.25, ...
    'start_synth',  1, ...               % idagrn6 seismograms of the step-1 starting model (all branches)
    ...                                  %   for the waveform figures and record sections
    ... % ---- step 2: inversion (D3) ----
    'sampler',  'pt', ...                % 'pt' (parallel tempering; default) or 'plain' (independent rjMcMC chains)
    'nuisance', 'marginal', ...          % tau_n(w), dth_n(w), B_n(w): 'marginal' (integrated out,
    ...                                  %   linearized) or 'sample' (tau, dth sampled; needs sampler 'plain')
    'relin',    1, ...                   % relinearize tau, dth and B every relin iterations (1: per proposal)
    'nuis_start', 'best', ...            % start of the nuisance fit (spec_refine), against cycle skipping of
    ...                                  %   tau_n(w): 'best' = for every record the best of 'direct' (whole band
    ...                                  %   at once: constant tau, then all terms), 'continuation_global' and
    ...                                  %   'continuation_banded' (frequency continuation, long to short periods,
    ...                                  %   with a tau basis over the whole band or over each stage's band)
    'fcont_stages', 4, ...               % stages of the frequency continuation
    'fix_tau', false, ...                % synthetic test only: tau_n(w) held at the true path delays, exactly at
    ...                                  %   every frequency (apply_fixed_tau.m; update_ref 'none' or 'lin')
    'fix_nuisance', false, ...           % synthetic test only: B, tau, dth held at 1, 0, 0 (true values when the
    ...                                  %   reference is the true model and data = modelled branches)
    'ckpt_every', 100, ...               % checkpoint every this many iterations (resume after an interruption)
    'jump_every', 100, ...               % every this many burn-in iterations, try moving each tau_n(w) by one
    ...                                  %   cycle at the short periods and keep it if log L improves (0 = off)
    'eta_free', false, ...               % eta: false = held at the reference model; true = a free
    ...                                  %   B-spline profile (prior U[etamin, etamax], default 0.80-1.05)
    'update_ref', 'burnin', ...          % 'burnin': recompute the reference (MINEOS + idagrn6) once at the
    ...                                  %   end of burn-in from the mean model of all chains; 'lin': only the
    ...                                  %   MINEOS linearization of the dispersion (waveform reference of
    ...                                  %   step 1 kept); 'none': keep step 1
    'lin_updates', 1, ...                % 'lin': number of linearization updates, evenly spaced over the burn-in
    'upd_frac', 0.2, ...                 % fraction of the burn-in (its end) averaged for the update
    'upd_ngn',  10, ...                  % Gauss-Newton steps per state to move tau, dth, B to the new reference
    'ref_prep', '', ...                  % optional prep file of another starting reference (D6)
    'spec_Kt',  3, 'tau_prior', 15, 'dtau_prior', 10, ...   % tau_n(w): Legendre terms, prior std (s)
    'windows',  'auto', ...              % 'auto': the window picked from the data per event by ASWMS-Q
    ...                                  %   (A3_run_autowinpick.m: event.winpara of the event files; picker
    ...                                  %   options through 'aswms', e.g. struct('cycle_before', 2, 'cycle_after', 5));
    ...                                  %   'fund': fundamental-mode window only, Delta/U_w to Delta/U_f + 60 s
    ...                                  %   (cuts the fast higher modes, incl. the unmodelled branches >= nm);
    ...                                  %   U_w = max fundamental-mode U of the reference + fund_margin
    ...                                  %   'single': one window per wave type from Delta/U_start (overtones +
    ...                                  %   fundamental mode, one noise level); 'split': overtone and fundamental-mode windows
    'smooth_sig', [2e-4 2e-4], ...       % smoothness prior [sig_1 sig_2] on the curvature of the two
    ...                                  %   smooth_param profiles: log p = -1/2 sum_c int (X_c''/sig_c)^2 dz/smooth_len.
    ...                                  %   sig is the STD of the allowed curvature, not a weight:
    ...                                  %     LARGER sig = WEAKER prior = LESS smoothing (rougher profiles allowed);
    ...                                  %     SMALLER sig = STRONGER prior = MORE smoothing (penalty scales as 1/sig^2,
    ...                                  %     so halving sig quadruples it); 0 = prior OFF for that profile (not
    ...                                  %     'infinitely smooth'). E.g. 4e-4 km/s/km^2: a curvature of 4e-4 over
    ...                                  %     smooth_len costs 1/2 in log-prior. See lr_common.m.
    ...                                  %   Default 2e-4 (4-branch PA5 tests: better than 4e-4 and than off)
    'smooth_param', 'vsv_vsh', ...       % profiles smoothed by smooth_sig:
    ...                                  %   'vsv_vsh' (default): [Vsv Vsh], both km/s/km^2; xi smoothed only through Vsh
    ...                                  %   'vsv_xi': [Vsv (km/s/km^2), xi (1/km^2)], independent; sig_xi ~ 2*sig_sv/Vsv
    ...                                  %     (about 2e-4 for 4e-4) is comparable to the 'vsv_vsh' prior
    'smooth_len', 10, ...                % depth interval (km) of one Gaussian curvature term (larger = weaker prior)
    'smooth_hier', false, ...            % true: HIERARCHICAL smoothing - the two curvature strengths are sampled (log-
    ...                                  %   uniform hyperprior on smooth_sig_range) and smooth_sig is not used; false:
    ...                                  %   smooth_sig as set (a user-chosen fixed strength, or 0 = off)
    'smooth_sig_range', [5e-5 5e-3], ... % range of the sampled curvature strengths (smooth_hier = true)
    'vsv_ref_sig', 'step1', ...          % prior toward the reference Vsv (km/s; 0 = off): N(Vsv_ref(z), vsv_ref_sig^2)
    ...                                  %   every vsv_ref_len km between the Moho and zmax, Vsv_ref = the reference
    ...                                  %   model of the prep file (step-1 model; true model in an oracle test).
    ...                                  %   LARGER vsv_ref_sig = WEAKER prior. xi is not affected (lr_common.m)
    'vsv_ref_len', 10, ...               % depth interval (km) of one Gaussian term of that prior (larger = weaker prior)
    ...                                  %   vsv_ref_sig = 'step1': width = the step-1 posterior std (s1_file s1_unc, from a
    ...                                  %   Bayesian step 1) x vsv_ref_scale, floored at vsv_ref_sig_min (vsv_ref_sigz.m);
    ...                                  %   DEFAULT 'step1': the step-1 (D1) posterior carried forward into step 2
    'vsv_ref_scale', 1, 'vsv_ref_sig_min', 0.01, ...
    'ximin', 0.85, 'ximax', 1.15, ...    % uniform prior on the xi control points (default symmetric about 1, prior
    ...                                  %   mean 1 = isotropic; runs before 11/2026 used 0.90-1.20, prior mean 1.05)
    'xi_ref_sig', 0, ...                 % Gaussian prior on xi toward xi_ref_mean (0 = OFF): N(xi_ref_mean, xi_ref_sig^2)
    ...                                  %   every xi_ref_len km between the Moho and zmax, on top of the uniform
    ...                                  %   ximin-ximax. LARGER xi_ref_sig = WEAKER prior (e.g. 0.03-0.05; lr_common.m)
    'xi_ref_mean', 1, ...                % its mean (1 = isotropic; a fixed value, not the reference model's xi)
    'xi_ref_len', 10, ...                % depth interval (km) of one Gaussian term of that prior (larger = weaker prior)
    'xi_ref_hier', true, ...             % true (default): HIERARCHICAL xi prior - its strength sigma_xi is sampled (log-uniform
    ...                                  %   hyperprior on xi_ref_sig_range), so the data set how strongly xi is pulled
    ...                                  %   toward xi_ref_mean; xi_ref_sig is then not used. false: xi_ref_sig as set
    ...                                  %   (a user-chosen fixed strength, or 0 = off)
    'xi_ref_sig_range', [0.005 0.2], ... % range of the sampled sigma_xi (xi_ref_hier = true)
    'waves', [1 2], ...                  % wave types inverted: 1 = Love, 2 = Rayleigh ([2] = Rayleigh only)
    'xi_fixed', false, ...               % true: xi held at the reference model, Vsv only (step-1 waveform version)
    'fund_margin', [0.4 0.2], ...        % 'fund' window start margin (km/s) [Love Rayleigh]
    'U_start',  [], ...                  % window start group velocity (km/s; scalar or [Love Rayleigh]);
    ...                                  %   [] = max U of branches 1..nm-1 of the reference + 0.6 km/s
    'exact_geom', true, ...              % moveout across the array from great-circle distances (false: local
    ...                                  %   plane wave in x-y, as before; event_G_spec.m)
    'N',        800, ...                 % samples per record (N*dt = 3200 s: earliest overtone to slowest fundamental)
    'noise',    0.03, ...                % synthetic noise std, fraction of the peak band-passed amplitude of each
    ...                                  %   event and wave type (A0; 'internal': of each record's peak)
    'nchains',  4, 'nit', 3000, 'burn', 1500, ...           % plain sampler
    'nens',     4, 'K', 6, 'Tmax', 20, 'pt_nit', 1500, 'pt_burn', 750, ... % parallel tempering (4 ensembles: R-hat between them)
    'parallel', [], ...                  % D3 chains/ensembles over a MATLAB pool: [] = MATLAB's rule (open or
    ...                                  %   auto-created pool), true = start one (par_workers.m), false = serial
    'par_workers', [], ...               % maximum pool workers for D3 ([] = one per chain/ensemble)
    'close_figs', [], ...                % close figures after saving them (and 'close all' at the start):
    ...                                  %   [] = auto (true headless: Octave or MATLAB -batch/-nodisplay;
    ...                                  %   false in interactive MATLAB), or true / false (fig_close_mode.m)
    'plot_fontsize', 15 ...              % font size (pt) of all plot axes (tick labels, axis labels, titles);
    ...                                  %   set as the default for new axes (plot_font_defaults.m)
    ));
plot_font_defaults(opt.plot_fontsize);
setappdata(0, 'RL_close_figs', opt.close_figs);      % read by fig_close_mode.m
pa5iso_ = strcmp(opt.s1_start_model, 'pa5_iso');
if ~pa5iso_ && ~strcmp(opt.s1_start_model, 'prem'), error('opt.s1_start_model must be ''pa5_iso'' or ''prem'''); end
if isempty(opt.s1_start), if pa5iso_, opt.s1_start = 'RLinv_pa5iso'; else, opt.s1_start = 'RLinv_prem_iso'; end, end
if isempty(opt.ref_qmod), if pa5iso_, opt.ref_qmod = opt.trueid; else, opt.ref_qmod = opt.premid; end, end
clear pa5iso_
% run tag, project tag and figure folder: set automatically unless given.
% An automatic value is remembered (opt.tag_auto_, opt.outdir_auto_) and
% recomputed on the next call, so an opt left in the (MATLAB) workspace by an
% earlier run does not keep a stale tag after opt.sampler changes; a tag or
% outdir set by hand is kept.
% With opt.project_tag set (by hand), the automatic tag is '_<project_tag>';
% otherwise '_<sampler>', and the project tag is the run tag without '_'.
pt_user_ = ~isempty(opt.project_tag) && ~(isfield(opt, 'project_tag_auto_') && strcmp(opt.project_tag, opt.project_tag_auto_));
if ~isfield(opt, 'tag') || (isfield(opt, 'tag_auto_') && strcmp(opt.tag, opt.tag_auto_))
    if pt_user_, opt.tag = ['_' opt.project_tag]; else, opt.tag = ['_' opt.sampler]; end
    opt.tag_auto_ = opt.tag;
end
if ~pt_user_
    opt.project_tag = regexprep(opt.tag, '^_', ''); opt.project_tag_auto_ = opt.project_tag;
end
clear pt_user_
if ~isfield(opt, 'outdir') || (isfield(opt, 'outdir_auto_') && strcmp(opt.outdir, opt.outdir_auto_))
    opt.outdir = ['figs' opt.tag]; opt.outdir_auto_ = opt.outdir;
end
clear here_ mdir_
