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
    'fcont_stages', 4, ...               % initial nuisance fit from long to short periods in this many band
    ...                                  %   stages, against cycle skipping of tau_n(w) (1 = off; spec_refine)
    'jump_every', 100, ...               % every this many burn-in iterations, try moving each tau_n(w) by one
    ...                                  %   cycle at the short periods and keep it if log L improves (0 = off)
    'eta_free', false, ...               % eta: false = held at the reference model; true = a free
    ...                                  %   B-spline profile (prior U[etamin, etamax], default 0.80-1.05)
    'update_ref', 'burnin', ...          % 'burnin': recompute the reference (MINEOS + idagrn6) once at the
    ...                                  %   end of burn-in from the mean model of all chains; 'lin': only the
    ...                                  %   MINEOS linearization of the dispersion (waveform reference of
    ...                                  %   step 1 kept); 'none': keep step 1
    'upd_frac', 0.2, ...                 % fraction of the burn-in (its end) averaged for the update
    'upd_ngn',  10, ...                  % Gauss-Newton steps per state to move tau, dth, B to the new reference
    'ref_prep', '', ...                  % optional prep file of another starting reference (D6)
    'spec_Kt',  3, 'tau_prior', 15, 'dtau_prior', 10, ...   % tau_n(w): Legendre terms, prior std (s)
    'windows',  'auto', ...              % 'auto': window picked from the data per event (auto_win_select.m, the
    ...                                  %   ASWMS-Q picker; options in 'autowin', e.g. struct('cycle_before', 2,
    ...                                  %   'cycle_after', 5, 'min_groupv', 2, 'max_groupv', 5));
    ...                                  %   'fund': fundamental-mode window only, Delta/U_w to Delta/U_f + 60 s
    ...                                  %   (cuts the fast higher modes, incl. the unmodelled branches >= nm);
    ...                                  %   U_w = max fundamental-mode U of the reference + fund_margin
    ...                                  %   'single': one window per wave type from Delta/U_start (overtones +
    ...                                  %   fundamental mode, one noise level); 'split': overtone and fundamental-mode windows
    'smooth_sig', [0 0], ...             % smoothness prior: std of the curvature of [Vsv Vsh] (km/s/km^2),
    ...                                  %   one Gaussian term per smooth_len km (0 = off; see lr_common.m)
    'smooth_len', 10, ...
    'autowin',  struct(), ...            % auto_win_select options ('auto' windows); empty = ASWMS-Q defaults
    'fund_margin', [0.4 0.2], ...        % 'fund' window start margin (km/s) [Love Rayleigh]
    'U_start',  [], ...                  % window start group velocity (km/s; scalar or [Love Rayleigh]);
    ...                                  %   [] = max U of branches 1..nm-1 of the reference + 0.6 km/s
    'N',        800, ...                 % samples per record (N*dt = 3200 s: earliest overtone to slowest fundamental)
    'noise',    0.03, ...                % noise std, fraction of each record's peak
    'nchains',  4, 'nit', 3000, 'burn', 1500, ...           % plain sampler
    'nens',     2, 'K', 6, 'Tmax', 20, 'pt_nit', 1500, 'pt_burn', 750 ... % parallel tempering
    ));
pa5iso_ = strcmp(opt.s1_start_model, 'pa5_iso');
if ~pa5iso_ && ~strcmp(opt.s1_start_model, 'prem'), error('opt.s1_start_model must be ''pa5_iso'' or ''prem'''); end
if isempty(opt.s1_start), if pa5iso_, opt.s1_start = 'RLinv_pa5iso'; else, opt.s1_start = 'RLinv_prem_iso'; end, end
if isempty(opt.ref_qmod), if pa5iso_, opt.ref_qmod = opt.trueid; else, opt.ref_qmod = opt.premid; end, end
clear pa5iso_
if ~isfield(opt, 'tag'), opt.tag = ['_' opt.sampler]; end
if ~isfield(opt, 'outdir'), opt.outdir = ['figs' opt.tag]; end
clear here_ mdir_
