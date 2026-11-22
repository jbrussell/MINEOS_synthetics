% D2_STEP2_PREP_MINEOS_IDAGRN  Step 2, part 1 (driver): MINEOS and idagrn6.
%
% With the step-1 model (CARDS/<opt.s1_final>.card, from D1) or, with
% opt.ref_model = 'pa5_iso', an isotropic PA5 without its lid/LVZ
% (functions/pa5_iso_ref_card.m) as the reference model, runs once and
% caches in opt.prepfile:
%   - MINEOS mode tables of the true model (PA5, branches
%     0..opt.nbranch_data-1) and of the reference (branches 0..opt.nm-1)
%   - reference phase/group velocities and Vsv, Vsh, Vpv, Vph, rho kernels
%     of branches 0..opt.nm-1 at the dispersion nodes
%   - idagrn6 seismograms at the 16 NoMelt OBSs for the 15 events:
%       data      : PA5, all opt.nbranch_data branches (T and Z)
%       Uref      : reference model, one branch at a time (idagrn6_mask)
%       PREM      : PREM, all branches (for the waveform figures)
% (functions/prep_mineos_idagrn.m). Takes a few minutes.
% Then record sections of the full mode-sum synthetics (PA5 data and PREM,
% T and Z, inversion band, with the window bounds) in figs_step2/
% (functions/plot_record_sections.m).
%
% Needs the MINEOS and idagrn6 binaries (../FORTRAN/bin) and SAC (sacio).
% Previous: D1_step1_rayleigh_isotropic_inv.m. Next: D3_step2_waveform_inv.m

parameters_RL;
o = opt;
if strcmp(opt.ref_model, 'pa5_iso')         % isotropic PA5 without lid/LVZ, PA5 ratios and eta
    ct = read_model_card(fullfile(opt.mdir, 'CARDS', [opt.trueid '.card']));
    o.ref_card = pa5_iso_ref_card(ct, opt.zmoho, opt.ref_zlin);
else
    S1 = load(opt.s1_file, 'card_final');
    o.ref_card = S1.card_final;             % the step-1 isotropic model
end
o.refid = opt.refid;
P = prep_mineos_idagrn(o);
fprintf('cached in %s: data branches 0-%d, modeled branches 0-%d\n', opt.prepfile, opt.nbranch_data-1, opt.nm-1);
plot_record_sections(P, 'figs_step2', optget(opt, 'fc', []), optget(opt, 'U_start', []));
fprintf('record sections in figs_step2/\n');
