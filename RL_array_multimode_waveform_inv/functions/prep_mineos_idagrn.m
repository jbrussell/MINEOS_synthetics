function P = prep_mineos_idagrn(opt)
% P = prep_mineos_idagrn(opt)
%
% Run MINEOS and idagrn6 for the recovery test and cache the results in
% opt.prepfile (default 'prep_mineos.mat'):
%
%  1. the reference ("input") card: by default PREM below the PA5 Moho with
%     the PA5 water, sediments and crust above it (mineos_tools/
%     make_ref_card.m), or the card opt.ref_card (a read_model_card
%     structure, e.g. an updated reference from a previous inversion stage);
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
%  5. (opt.start_synth = 1) full synthetics of the step-1 starting card (opt.s1_start,
%     written by D1, same branches as the data), for comparison
%     in the waveform figures.
%
% Steps 2-4 for the true model are skipped when opt.P_true (a previous
% output) is given; only the reference part is then recomputed.
%
% Needs the MINEOS and idagrn6 binaries in ../FORTRAN/bin and SAC (sacio)
% for idagrn6; see ../FORTRAN/README. Run from this directory.
%
% opt fields (all optional): prepfile, mdir (path to run_MINEOS), trueid
% ('pa5_5km'), premid ('prem_35'), refid ('prem_nomelt'), ref_card, zmoho
% (12 km below sea level), maxF (70 mHz), nbranch_data (8), nm (3), Tdisp,
% dt (4 s), len_hr (1 h), eventfile, stationfile, P_true, start_synth (0),
% premref ('prem_nomelt')
%
% jbrussell 11/2026

if nargin < 1, opt = struct(); end
og = @(f, d) optget(opt, f, d);
prepfile = og('prepfile', 'prep_mineos.mat');
if exist(prepfile, 'file') == 2
    P = load(prepfile); P = P.P;
    fprintf('loaded %s\n', prepfile);
    return
end
here = fileparts(fileparts(mfilename('fullpath')));   % package folder
P.mdir = og('mdir', fullfile(here, '..', 'run_MINEOS'));
addpath(fullfile(P.mdir, 'functions'));
P.trueid = og('trueid', 'pa5_5km');
P.premid = og('premid', 'prem_35');
P.ref_qmod = og('ref_qmod', P.premid);       % Q of the reference card
P.refid  = og('refid', 'prem_nomelt');
P.zmoho  = og('zmoho', 12);
P.maxF   = og('maxF', 70);
P.nbranch_data = og('nbranch_data', 8);
P.nm     = og('nm', 3);
P.Tdisp  = og('Tdisp', logspace(log10(0.95*17), log10(1.02*112), 10));
P.dt     = og('dt', 4);
P.len_hr = og('len_hr', 1);

% stations and events
fid = fopen(og('stationfile', fullfile(here, 'stations_nomelt.txt')));
S = textscan(fid, '%s %f %f %f'); fclose(fid);
P.sta = struct('name', {S{1}}, 'lat', S{2}, 'lon', S{3}, 'elev', S{4});
P.ev = read_cmtsolution(og('eventfile', fullfile(here, 'events_nomelt_2012.cmt')));

% 1. reference card
rc = og('ref_card', []);
if isempty(rc)
    P.card_ref = make_ref_card(P.mdir, P.refid, P.premid, P.trueid, P.zmoho);
else
    write_MINEOS_mod(rc, fullfile(P.mdir, 'CARDS', [P.refid, '.card']));
    copyfile(fullfile(P.mdir, 'CARDS', [P.ref_qmod, '.qmod']), fullfile(P.mdir, 'CARDS', [P.refid, '.qmod']));
    P.card_ref = read_model_card(fullfile(P.mdir, 'CARDS', [P.refid, '.card']));
end
P.card_true = read_model_card(fullfile(P.mdir, 'CARDS', [P.trueid, '.card']));

% 2-3. mode tables, dispersion and kernels
Pt = og('P_true', []);
tic;
for TY = 'TS'
    if isempty(Pt) || size(Pt.(['phv_true_', TY]), 2) < P.nm
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

% 4. synthetics
for TY = 'TS'
    if isempty(Pt)
        D = idagrn_synth(P.mdir, P.trueid, TY, P.maxF, P.nbranch_data, P.ev, P.sta, P.len_hr, P.dt, []);
        P.t = D.t;
        P.(['u_true_', TY]) = D.u;                   % {nev x 1} [Nt x Nsta]
    else
        P.t = Pt.t; P.(['u_true_', TY]) = Pt.(['u_true_', TY]);
    end
    D = idagrn_synth(P.mdir, P.refid, TY, P.maxF, P.nm, P.ev, P.sta, P.len_hr, P.dt, 0:P.nm-1);
    P.(['u_ref_', TY]) = D.u;                        % {nev x nm}
    fprintf('idagrn %c done (%.0f s)\n', TY, toc);
end

% 5. synthetics of the step-1 STARTING model (CARDS/<s1_start>.card, with
% its .qmod, written by D1), all branches (for the waveform figures)
if og('start_synth', 0)
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

tmpf = tmp_name('prep');
save(tmpf, 'P'); movefile(tmpf, prepfile);
end
