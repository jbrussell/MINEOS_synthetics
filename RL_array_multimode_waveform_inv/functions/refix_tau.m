function rec = refix_tau(rec, cfg, Ctrue_f, Cref_f)
% rec = refix_tau(rec, cfg, Ctrue_f, Cref_f)
%
% Synthetic test only (opt.fix_tau = true) with a full reference update
% (opt.update_ref = 'burnin', D3): apply_reference.m installs the new
% reference spectra without the true path delays of apply_fixed_tau.m, so
% they are applied again, now relative to the new reference:
%   rec(r).Uref{n} <- rec(r).Uref{n} .* exp(-i w tau_n(w)),
%   tau = Del0 (1/C_true - 1/C_ref_new)   (Ctrue_f, Cref_f: {Love, Rayleigh}
%                                          phase velocities at the band frequencies)
% The sampled tau_n(w) stay held at zero (mc.taup of apply_fixed_tau.m).
w = cfg.w(:);
for r = 1:numel(rec)
    wv = rec(r).wave;
    td = rec(r).Del0 * (1./Ctrue_f{wv} - 1./Cref_f{wv});   % [Nf x nm] (s)
    for n = 1:numel(rec(r).Uref)
        tn = td(:, n); tn(~isfinite(tn)) = 0;
        rec(r).Uref{n} = rec(r).Uref{n} .* exp(-1i * w .* tn);
    end
end
end
