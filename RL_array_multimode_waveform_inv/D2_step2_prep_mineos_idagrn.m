% D2_STEP2_PREP_MINEOS_IDAGRN  Step 2, part 1 (driver): MINEOS and idagrn6.
%
% With the step-1 model (CARDS/<opt.s1_final>.card, from D1) as the
% reference model, runs once and caches in opt.prepfile:
%   - MINEOS mode tables of the true model (PA5, branches
%     0..opt.nbranch_data-1) and of the reference (branches 0..opt.nm-1)
%   - reference phase/group velocities and Vsv, Vsh, Vpv, Vph, rho kernels
%     of branches 0..opt.nm-1 at the dispersion nodes
%   - idagrn6 seismograms at the 16 NoMelt OBSs for the 15 events:
%       data      : PA5, all opt.nbranch_data branches (T and Z)
%       Uref      : reference model, one branch at a time (idagrn6_mask)
%       PREM      : PREM, all branches (for the waveform figures)
% (functions/prep_mineos_idagrn.m). Takes a few minutes.
%
% Needs the MINEOS and idagrn6 binaries (../FORTRAN/bin) and SAC (sacio).
% Previous: D1_step1_rayleigh_isotropic_inv.m. Next: D3_step2_waveform_inv.m

parameters_RL;
S1 = load(opt.s1_file, 'card_final');
o = opt;
o.ref_card = S1.card_final;                 % the step-1 isotropic model
o.refid = opt.s1_final;
P = prep_mineos_idagrn(o);
fprintf('cached in %s: data branches 0-%d, modeled branches 0-%d\n', opt.prepfile, opt.nbranch_data-1, opt.nm-1);
