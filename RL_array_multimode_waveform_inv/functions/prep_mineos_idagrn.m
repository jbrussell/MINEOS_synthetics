function P = prep_mineos_idagrn(opt)
% P = prep_mineos_idagrn(opt)
%
% Run MINEOS and idagrn6 for the recovery test and cache the results in
% opt.prepfile (default 'prep_mineos.mat'):
%
%  1. the reference ("input") card opt.ref_card (a read_model_card
%     structure: the step-1 model in D2, or an updated reference from a
%     previous inversion stage), written to CARDS/<opt.refid>.card with the
%     Q of CARDS/<opt.ref_qmod>.qmod;
%  2. MINEOS mode tables of the true model (PA5; opt.nbranch_data branches)
%     and of the reference model (branches 0..nm-1), Love (T) and Rayleigh
%     (S), up to opt.maxF mHz;
%  3. reference phase/group velocities and MINEOS phase-velocity kernels of
%     branches 0..nm-1 at the dispersion nodes Tdisp; true (PA5) phase
%     velocities for comparison;
%  4. idagrn6 synthetics at the 16 NoMelt stations for the 15 events of
%     events_nomelt_2012.cmt: the full PA5 synthetics (all branches in the
%     table: the "data") and single-branch synthetics of the reference model
%     for branches 0..nm-1 (idagrn6_mask: the reference spectra Uref),
%     transverse (T) and vertical (Z) components, displacement;
%  5. (opt.start_synth = 1) full synthetics of the step-1 starting card
%     (opt.s1_start, written by D1, same branches as the data), for the
%     record sections (D2) and the waveform figures (fig4).
%
% Steps 2-4 for the true model are skipped when opt.P_true (a previous
% output) is given; only the reference part is then recomputed.
%
% opt.data_source = 'aswms' (default of parameters_RL.m): the data are the
% ASWMS-Q event files (A1-A3; for the synthetic test from the SAC files of
% A0_synth_make_sac_data.m), so the full true-model synthetics of step 4
% are NOT computed here (P.u_true_* absent). The true-model MINEOS phase
% velocities (diagnostics of the synthetic test) are computed when
% opt.trueid is not empty (set it to '' for real data). The single-branch
% reference synthetics are also written in the ASWMS-Q event-file format
% (write_branchmat.m): <opt.aswms_dir>branchmat/<refid>/<event id>_<component>.mat.
% opt.lin_only = true: MINEOS (dispersion, kernels) only, no idagrn6
% synthetics (linearization-only reference updates, D3).
%
% Needs the MINEOS and idagrn6 binaries in ../FORTRAN/bin and SAC (sacio)
% for idagrn6; see ../FORTRAN/README. Run from this directory.
%
% opt fields: ref_card (required), refid ('RLinv_step1'), ref_qmod
% (trueid); optional: prepfile, mdir (path to run_MINEOS), trueid
% ('pa5_5km'), zmoho (12 km below sea level), maxF (70 mHz), nbranch_data
% (8), nm (3), Tdisp, dt (4 s), len_hr (1 h), eventfile, stationfile, P_true,
% start_synth (0), s1_start ('RLinv_pa5iso')
%
% jbrussell 11/2026

if nargin < 1, opt = struct(); end
og = @(f, d) optget(opt, f, d);
prepfile = og('prepfile', 'prep_mineos.mat');
here = fileparts(fileparts(mfilename('fullpath')));   % package folder
if exist(prepfile, 'file') == 2
    P = load(prepfile); P = P.P;
    fprintf('loaded %s\n', prepfile);
    % prep files made before the event reader had the centroid time t_cent
    % (needed by read_event_data.m): re-read the events (cheap, same order)
    if ~isfield(P.ev, 't_cent')
        ev_ = rl_read_cmtsolution(og('eventfile', fullfile(here, 'events_nomelt_2012.cmt')));
        if numel(ev_) ~= numel(P.ev) || ~all(strcmp({ev_.id}, {P.ev.id}))
            error('prep_mineos_idagrn: %s is from another event list; delete it and rerun D2', prepfile);
        end
        P.ev = ev_;
        fprintf('  (event list re-read: %s predates t_cent)\n', prepfile);
    end
    return
end
P.mdir = og('mdir', fullfile(here, '..', 'run_MINEOS'));
addpath(fullfile(P.mdir, 'functions'));
P.trueid = og('trueid', 'pa5_5km');
P.ref_qmod = og('ref_qmod', P.trueid);       % Q of the reference card
P.refid  = og('refid', 'RLinv_step1');
P.zmoho  = og('zmoho', 12);
P.maxF   = og('maxF', 70);
P.nbranch_data = og('nbranch_data', 8);
P.nm     = og('nm', 3);
fc_ = wave_band(og('fc', []), 0);              % waveform band(s) (Hz): the dispersion nodes span both
P.Tdisp  = og('Tdisp', logspace(log10(0.95/fc_(4)), log10(1.02/fc_(1)), 10));
P.dt     = og('dt', 4);
P.len_hr = og('len_hr', 1);

% stations and events
fid = fopen(og('stationfile', fullfile(here, 'stations_nomelt.txt')));
S = textscan(fid, '%s %f %f %f'); fclose(fid);
P.sta = struct('name', {S{1}}, 'lat', S{2}, 'lon', S{3}, 'elev', S{4});
P.ev = rl_read_cmtsolution(og('eventfile', fullfile(here, 'events_nomelt_2012.cmt')));

% 1. reference card
rc = og('ref_card', []);
if isempty(rc), error('prep_mineos_idagrn: opt.ref_card (the reference model, e.g. the step-1 card) is required'); end
ws_ = warning;                                  % write_MINEOS_mod turns all warnings on
write_MINEOS_mod(rc, fullfile(P.mdir, 'CARDS', [P.refid, '.card']));
warning(ws_);
copyfile(fullfile(P.mdir, 'CARDS', [P.ref_qmod, '.qmod']), fullfile(P.mdir, 'CARDS', [P.refid, '.qmod']));
P.card_ref = read_model_card(fullfile(P.mdir, 'CARDS', [P.refid, '.card']));
if ~isempty(P.trueid)
    P.card_true = read_model_card(fullfile(P.mdir, 'CARDS', [P.trueid, '.card']));
else
    P.card_true = P.card_ref;                     % (real data: no true model; placeholder for the plots)
end

% MINEOS/idagrn6 runs share scratch files in P.mdir: serialized between
% processes (mineos_lock.m) from here until the synthetics are done
mineos_lock(P.mdir, 'acquire');
try
% 2-3. mode tables, dispersion and kernels
Pt = og('P_true', []);
aswms_ = strcmp(og('data_source', 'internal'), 'aswms');
true_ = ~isempty(P.trueid);
tic;
for TY = 'TS'
    if ~true_
        P.(['phv_true_', TY]) = NaN(numel(P.Tdisp), P.nm);           % (real data: no true model)
    elseif isempty(Pt) || size(Pt.(['phv_true_', TY]), 2) < P.nm
        o = mineos_modes(P.mdir, P.trueid, TY, P.nbranch_data, P.maxF, P.Tdisp, []);
        P.(['phv_true_', TY]) = o.phv(:, 1:P.nm);
    else
        P.(['phv_true_', TY]) = Pt.(['phv_true_', TY]);
    end
    o = mineos_modes(P.mdir, P.refid, TY, P.nm, P.maxF, P.Tdisp, 0:P.nm-1);
    P.(['phv_ref_', TY]) = o.phv(:, 1:P.nm);
    P.(['grv_ref_', TY]) = o.grv(:, 1:P.nm);
    P.(['K_ref_', TY]) = o.K;
    fprintf('MINEOS %c done (%.0f s)\n', TY, toc);
end

% 4. synthetics (none with opt.lin_only: a linearization-only reference
% update needs the MINEOS dispersion and kernels only)
lin_only_ = og('lin_only', false);
for TY = 'TS'
    if lin_only_, break; end
    if aswms_
        % the data are the ASWMS-Q event files (no true-model synthetics here)
    elseif isempty(Pt)
        D = idagrn_synth(P.mdir, P.trueid, TY, P.maxF, P.nbranch_data, P.ev, P.sta, P.len_hr, P.dt, []);
        P.t = D.t;
        P.(['u_true_', TY]) = D.u;                   % {nev x 1} [Nt x Nsta]
    else
        P.t = Pt.t; P.(['u_true_', TY]) = Pt.(['u_true_', TY]);
    end
    D = idagrn_synth(P.mdir, P.refid, TY, P.maxF, P.nm, P.ev, P.sta, P.len_hr, P.dt, 0:P.nm-1);
    P.(['u_ref_', TY]) = D.u;                        % {nev x nm}
    P.t = D.t;
    fprintf('idagrn %c done (%.0f s)\n', TY, toc);
end

% 5. synthetics of the step-1 STARTING model (CARDS/<s1_start>.card, with
% its .qmod, written by D1), all branches (for the waveform figures)
if og('start_synth', 0) && ~lin_only_
    sid = og('s1_start', 'RLinv_pa5iso');
    if exist(fullfile(P.mdir, 'CARDS', [sid, '.card']), 'file') ~= 2
        error('prep_mineos_idagrn: %s.card (the step-1 starting model) not found; run D1 first', sid);
    end
    for TY = 'TS'
        mineos_modes(P.mdir, sid, TY, P.nbranch_data, P.maxF, P.Tdisp, []);
        D = idagrn_synth(P.mdir, sid, TY, P.maxF, P.nbranch_data, P.ev, P.sta, P.len_hr, P.dt, []);
        P.(['u_start_', TY]) = D.u;
        fprintf('step-1 start idagrn %c done (%.0f s)\n', TY, toc);
    end
end

if (aswms_ || ~isempty(og('aswms_dir', ''))) && ~lin_only_   % reference branch synthetics, ASWMS-Q event-file format
    write_branchmat(P, opt);
end
catch err
    mineos_lock(P.mdir, 'release'); rethrow(err);
end
mineos_lock(P.mdir, 'release');
tmpf = tmp_name('prep');
save(tmpf, 'P'); movefile(tmpf, prepfile);
end
