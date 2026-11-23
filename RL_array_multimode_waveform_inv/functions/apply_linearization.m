function cfg = apply_linearization(cfg, P)
% cfg = apply_linearization(cfg, P)
%
% Install ONLY the MINEOS linearization of the reference card of the prep
% structure P (prep_mineos_idagrn.m / update_reference.m): cfg.lin (phase
% velocities and kernels of P.card_ref, mineos_lin_setup.m), used by
% mineos_dispersion.m to predict the dispersion C(m) of a model.
%
% Unlike apply_reference.m, the waveform reference is NOT changed: the
% records keep their reference spectra rec.Uref and phase velocities
% rec.cref, and cfg.need, cfg.vbot, cfg.xbot, cfg.ebot and cfg.card_ref stay
% those of the setup. The nuisance parameters tau_n(w), dth_n(w), B_n(w)
% are relative to that waveform reference, so they keep their meaning and
% the likelihood of a given dispersion is unchanged; only the mapping from
% the model to its dispersion becomes more accurate near P.card_ref.
%
% Used for opt.update_ref = 'lin' (D3_step2_waveform_inv.m) and its figures.
lin = mineos_lin_setup(P.card_ref, P.K_ref_T, P.K_ref_S, P.phv_ref_T, P.phv_ref_S, ...
                       cfg.zw, cfg.zmin, cfg.zmax);
need = isfinite(P.phv_ref_T) & isfinite(P.phv_ref_S);
if any(cfg.need(:) & ~need(:))
    warning('apply_linearization:need', ...
        'the new linearization lacks %d (period, branch) nodes of the reference; predictions there are NaN', ...
        sum(cfg.need(:) & ~need(:)));
end
cfg.lin = lin;
cfg.lin_card = P.card_ref;                  % card of the linearization (for reference)
end
