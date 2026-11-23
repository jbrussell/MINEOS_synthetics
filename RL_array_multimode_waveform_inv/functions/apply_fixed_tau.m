function [rec, mc] = apply_fixed_tau(rec, cfg, pd, mc)
% [rec, mc] = apply_fixed_tau(rec, cfg, pd, mc)
%
% Synthetic test only (opt.fix_tau = true): hold the nuisance time shifts
% tau_n(w) at the TRUE path delays, exactly at every frequency of the band
% (not through the tau basis cfg.Psit), to test how much of the bias and
% of the convergence problems come from the time shifts (cycle skipping,
% the limited 3-term basis):
%   - each reference spectrum is multiplied by exp(-i w tau_true_n(w)),
%     tau_true = pd.true_delay{r} = Del0 (1/C_true - 1/C_ref) (MINEOS, at
%     the array centre), i.e. the reference waveforms are moved to the
%     true path delay;
%   - the sampled tau_n(w) (now relative to that) are held at zero by a
%     prior std of 1e-3 s on all its terms (mc.taup).
% B_n(w) and the arrival angles stay free. With the true dispersion the
% prediction then has the true phase at every station (plane-wave moveout
% across the array). Works with opt.update_ref = 'none' or 'lin' (a full
% reference update would replace the reference spectra).
w = cfg.w(:);
for r = 1:numel(rec)
    td = pd.true_delay{r};                       % [Nf x nm] (s)
    for n = 1:numel(rec(r).Uref)
        tn = td(:, n); tn(~isfinite(tn)) = 0;
        rec(r).Uref{n} = rec(r).Uref{n} .* exp(-1i * w .* tn);
    end
end
mc.taup = 1e-3 * ones(size(mc.taup));
mc.fix_tau = true;
end
