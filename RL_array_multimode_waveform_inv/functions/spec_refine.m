function lin = spec_refine(lin, rec, Cf, sig, rho, cfg, mc, K, Ka, dth_prior, ngn)
% SPEC_REFINE  Extend a pilot linearization point (pt_pilot_LR: tau, dth
% coefficients, scalar b) to the event_G_spec parameterization (K spectral nodes,
% Ka arrival-angle coefficients) by damped Gauss-Newton at the pilot
% structure Cf = {Love, Rayleigh} band dispersion. Records need .Uref, .cref.
%
% Output: lin.tau0 nrec x nm, lin.dth0 nrec x nm x Ka, lin.b0 nrec x nm x K,
% plus lin.K, lin.Ka, lin.dth_prior (run_pt_LR switches to event_G_spec
% when lin.K exists).
if nargin < 11, ngn = 12; end
nrec = numel(rec); nm = size(lin.tau0, 2);
[cfg.Phi, cfg.Psi] = spec_basis(cfg.w, K, Ka);
T0 = lin.tau0; D0 = zeros(nrec, nm, Ka); B0 = zeros(nrec, nm, K);
for r = 1:nrec
    ev = rec(r); sg = sig(2*ev.wave-1:2*ev.wave);
    L.tau0 = lin.tau0(r, :);
    if isfield(lin, 'dthh')                       % pilot already has dth(w) terms
        L.dth0 = [lin.dth0(r, :)', reshape(lin.dthh(r, :, 1:Ka-1), nm, Ka-1)];
    else
        L.dth0 = [lin.dth0(r, :)', zeros(nm, Ka-1)];
    end
    L.b0   = repmat(lin.b0(r, :).', 1, K);
    for itr = 1:ngn
        G = event_G_spec(ev, Cf{ev.wave}, L, cfg, sg, mc.sb);
        [mu, sb2, it, ith] = spec_layout(L, nm, K, Ka, mc.sb, mc.tau_prior, dth_prior);
        [~, xh] = event_loglike_ar(ev.d, G, ev.wt, ev.first, sig, rho, mu, sb2);
        nb = 2*nm*K;
        L.b0   = reshape(xh(1:2:nb) + 1i*xh(2:2:nb), K, nm).';
        L.tau0 = L.tau0 + max(min(xh(it)', 2), -2);
        L.dth0 = L.dth0 + max(min(xh(ith), 3), -3);
    end
    T0(r, :) = L.tau0; D0(r, :, :) = reshape(L.dth0, 1, nm, Ka); B0(r, :, :) = reshape(L.b0, 1, nm, K);
end
lin.tau0 = T0; lin.dth0 = D0; lin.b0 = B0;
lin.K = K; lin.Ka = Ka; lin.dth_prior = dth_prior;
end
