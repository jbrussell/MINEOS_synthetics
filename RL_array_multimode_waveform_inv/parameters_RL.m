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
    ... % ---- models and files ----
    'mdir',      mdir_, ...              % run_MINEOS (MINEOS cards, drivers, ../FORTRAN/bin)
    'trueid',    'pa5_5km', ...          % true model (CARDS/<trueid>.card)
    'premid',    'prem_35', ...          % PREM card
    'premref',   'prem_nomelt', ...      % PREM below the true Moho, true water/crust above (make_ref_card)
    'zmoho',     12, ...                 % Moho of the true model (km below sea level)
    'eventfile', fullfile(here_, 'events_nomelt_2012.cmt'), ...
    'stationfile', fullfile(here_, 'stations_nomelt.txt'), ...
    ... % ---- step 1: isotropic Rayleigh inversion (D1) ----
    's1_periods',   logspace(log10(20), log10(100), 12), ...   % (s)
    's1_cstd',      0.005, ...           % data uncertainty (fraction of c)
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
    ... % ---- data: ASWMS-Q event files (D0 synthetic only; A1-A3 for synthetic and real data) ----
    'data_source',  'aswms', ...         % 'aswms': the data are read from the ASWMS-Q event files
    ...                                  %   <aswms_dir>eventmat/<event id>_<component>.mat (A1-A3; for the
    ...                                  %   synthetic test written from the SAC files of D0); 'internal': the
    ...                                  %   previous behaviour (true-model synthetics computed in step 2, noise
    ...                                  %   added to the processed records)
    'aswms_dir',    ['ASWMS' filesep], ... % ASWMS-Q working directory (sacdata/, eventmat/, winpara/, branchmat/)
    'aswms_components', {{'BHT', 'BHZ'}}, ... % channel names of the [Love Rayleigh] records (SAC KCMPNM)
    'aswms_net',    'XX', ...            % network code of the synthetic SAC files
    'aswms',        struct(), ...        % overrides of the ASWMS-Q parameters (aswms_parameters.m), e.g.
    ...                                  %   struct('snr_tol', 3, 'maxdepth', 50, 'min_sta_num', 5)
    'synth_len_hr', 2, ...               % D0: length of the synthetic SAC records (h; ASWMS-Q needs dist/2 km/s)
    'synth_seed',   3, ...               % D0: random seed of the noise
    'truthfile',    'synth_truth.mat', ... % D0: noise-free records, noise and true model (in aswms_dir)
    'noise_guess',  0.03, ...            % real data without noise levels: starting noise level (fraction of peak)
    'noise_levels', 'pre', ...           % 'pre': noise std of each record and station measured before the first
    ...                                  %   arrival (dist/13.5 km/s - 60 s after the origin, plus any samples
    ...                                  %   before the origin); the sampled sigma is a scale factor on it;
    ...                                  %   'none': one sigma per wave type
    'noise_model',  'ar1', ...           % noise correlation: 'ar1' (AR(1), correlation sampled); EXPERIMENTAL, fixed
    ...                                  %   from the pre-P noise: 'ar2' (AR(2) per wave type, pooled), 'arp' (AR(p),
    ...                                  %   pooled), 'emp' (AR(p) per station and record). Test (band-passed noise,
    ...                                  %   4 s sampling): none whitens it well (residual lag-1 correlation 0.6-0.65
    ...                                  %   vs 0.56 for AR(1)); the noise covariance is nearly singular (the band is
    ...                                  %   ~1/3 of the Nyquist range), high orders are unstable.
    ...                                  %   'fdom': frequency-domain likelihood (tapered in-band DFT coefficients of
    ...                                  %   each station's window, scaled by its noise spectrum from the pre-P
    ...                                  %   noise; fdom_operator.m); for real data
    'noise_psd',    'shape', ...         % 'fdom' spectra: 'shape' (pooled shape x station level) or 'station'
    'noise_arp',    8, ...               % 'arp', 'emp': AR order
    'noise_nugget', 0.05, ...            % white-noise fraction added before the AR fits (stability)
    'synth_pre_s',  0, ...               % D0: noise before the origin (s; 0: records start at the origin)
    'synth_noise_spread', 0, ...         % D0: noise level per station x 10^(spread*U[-1,1]) (0: all equal)
    ... % ---- step 2: MINEOS + idagrn6 (D2) ----
    'prepfile',     'prep_step2.mat', ...
    'refid',        'RLinv_step1', ...   % reference card of step 2 (= s1_final)
    'ref_qmod',     '', ...              % Q (.qmod) of the step-1 and step-2 cards ('' = trueid for
    ...                                  %   s1_start_model = 'pa5_iso', premid for 'prem')
    'nbranch_data', 16, ...              % branches 0..15 in the PA5 data
    'nm',           8, ...               % branches 0..7 modeled
    'maxF',         70, ...              % mode tables up to (mHz)
    'dt',           4, 'len_hr', 1.25, ...
    'start_synth',  1, ...               % idagrn6 seismograms of the step-1 starting model (all branches)
    ...                                  %   for the waveform figures and record sections
    ... % ---- step 2: inversion (D3) ----
    'sampler',  'plain', ...             % 'plain' (independent rjMcMC chains) or 'pt' (parallel tempering)
    'nuisance', 'marginal', ...          % tau_n(w), dth_n(w), B_n(w): 'marginal' (integrated out,
    ...                                  %   linearized) or 'sample' (tau, dth sampled; plain only)
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
    'smooth_sig', [0 0], ...             % smoothness prior [sig_1 sig_2] on the curvature of the two
    ...                                  %   smooth_param profiles: log p = -1/2 sum_c int (X_c''/sig_c)^2 dz/smooth_len.
    ...                                  %   sig is the STD of the allowed curvature, not a weight:
    ...                                  %     LARGER sig = WEAKER prior = LESS smoothing (rougher profiles allowed);
    ...                                  %     SMALLER sig = STRONGER prior = MORE smoothing (penalty scales as 1/sig^2,
    ...                                  %     so halving sig quadruples it); 0 = prior OFF for that profile (not
    ...                                  %     'infinitely smooth'). E.g. 4e-4 km/s/km^2: a curvature of 4e-4 over
    ...                                  %     smooth_len costs 1/2 in log-prior. See lr_common.m.
    'smooth_param', 'vsv_vsh', ...       % profiles smoothed by smooth_sig:
    ...                                  %   'vsv_vsh' (default): [Vsv Vsh], both km/s/km^2; xi smoothed only through Vsh
    ...                                  %   'vsv_xi': [Vsv (km/s/km^2), xi (1/km^2)], independent; sig_xi ~ 2*sig_sv/Vsv
    ...                                  %     (about 2e-4 for 4e-4) is comparable to the 'vsv_vsh' prior
    'smooth_len', 10, ...                % depth interval (km) of one Gaussian curvature term (larger = weaker prior)
    'vsv_ref_sig', 0, ...                % prior toward the reference Vsv (km/s; 0 = off): N(Vsv_ref(z), vsv_ref_sig^2)
    ...                                  %   every vsv_ref_len km between the Moho and zmax, Vsv_ref = the reference
    ...                                  %   model of the prep file (step-1 model; true model in an oracle test).
    ...                                  %   LARGER vsv_ref_sig = WEAKER prior. xi is not affected (lr_common.m)
    'vsv_ref_len', 10, ...               % depth interval (km) of one Gaussian term of that prior (larger = weaker prior)
    'fund_margin', [0.4 0.2], ...        % 'fund' window start margin (km/s) [Love Rayleigh]
    'U_start',  [], ...                  % window start group velocity (km/s; scalar or [Love Rayleigh]);
    ...                                  %   [] = max U of branches 1..nm-1 of the reference + 0.6 km/s
    'exact_geom', true, ...              % moveout across the array from great-circle distances (false: local
    ...                                  %   plane wave in x-y, as before; event_G_spec.m)
    'N',        800, ...                 % samples per record (N*dt = 3200 s: earliest overtone to slowest fundamental)
    'noise',    0.03, ...                % synthetic noise std, fraction of the peak band-passed amplitude of each
    ...                                  %   event and wave type (D0; 'internal': of each record's peak)
    'nchains',  4, 'nit', 3000, 'burn', 1500, ...           % plain sampler
    'nens',     2, 'K', 6, 'Tmax', 20, 'pt_nit', 1500, 'pt_burn', 750, ... % parallel tempering
    'plot_fontsize', 15 ...              % font size (pt) of all plot axes (tick labels, axis labels, titles);
    ...                                  %   set as the default for new axes (plot_font_defaults.m)
    ));
plot_font_defaults(opt.plot_fontsize);
pa5iso_ = strcmp(opt.s1_start_model, 'pa5_iso');
if ~pa5iso_ && ~strcmp(opt.s1_start_model, 'prem'), error('opt.s1_start_model must be ''pa5_iso'' or ''prem'''); end
if isempty(opt.s1_start), if pa5iso_, opt.s1_start = 'RLinv_pa5iso'; else, opt.s1_start = 'RLinv_prem_iso'; end, end
if isempty(opt.ref_qmod), if pa5iso_, opt.ref_qmod = opt.trueid; else, opt.ref_qmod = opt.premid; end, end
clear pa5iso_
if ~isfield(opt, 'tag'), opt.tag = ['_' opt.sampler]; end
if ~isfield(opt, 'outdir'), opt.outdir = ['figs' opt.tag]; end
clear here_ mdir_
