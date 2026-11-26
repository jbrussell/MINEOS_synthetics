function P = update_reference(vsv, vsh, zp, cfg, opt, refid_new, prepfile_new, eta)
% P = update_reference(vsv, vsh, zp, cfg, opt, refid_new, prepfile_new, eta)
%
% New reference model for the step-2 waveform inversion: the current
% reference card (cfg.card_ref) with Vsv and Vsh between cfg.zmin and
% cfg.zmax replaced by the profiles vsv(zp), vsh(zp) (km/s, zp in km below
% the seafloor; e.g. a posterior mean), and Vpv, Vph and density scaled
% with the reference ratios (ref_card_update.m); eta(zp) replaces eta if
% given (when eta is a model parameter). MINEOS (mode tables,
% phase/group velocities and kernels of branches 0..nm-1) and idagrn6_mask
% (single-branch reference synthetics) are rerun for it with
% prep_mineos_idagrn.m; the true-model data and the PREM synthetics of
% opt.prepfile are reused, not recomputed. The result is cached in
% prepfile_new (CARDS/<refid_new>.card); install it with apply_reference.m.
%
% The MINEOS and idagrn6 runs share scratch files in opt.mdir, so they are
% serialized between processes (mineos_lock.m, in prep_mineos_idagrn.m).
%
% Used at the end of burn-in by D3_step2_waveform_inv.m
% (opt.update_ref = 'burnin') and by D6_step2_update_reference_and_rerun.m.
if exist(prepfile_new, 'file') == 2
    P = load(prepfile_new); P = P.P; fprintf('loaded %s\n', prepfile_new); return
end
Pold = prep_mineos_idagrn(opt);                 % cached: true data, step-1 start synthetics
if nargin < 8, eta = []; end
c = ref_card_update(cfg.card_ref, cfg.zw, cfg.zmin, cfg.zmax, zp, vsv, vsh, eta);
o = opt; o.prepfile = prepfile_new; o.refid = refid_new; o.ref_card = c;
o.P_true = Pold; o.start_synth = 0;
P = prep_mineos_idagrn(o);                       % (MINEOS/idagrn6 serialized inside: mineos_lock.m)
for TY = 'TS'                                    % step-1 start synthetics for the figures
    f = ['u_start_' TY];
    if isfield(Pold, f), P.(f) = Pold.(f); end
end
save(prepfile_new, 'P');
end

