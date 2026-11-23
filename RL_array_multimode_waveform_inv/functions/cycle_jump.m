function [lr, L, G, xh, njump] = cycle_jump(ev, Cf, lr, L, G, xh, cfg, sig, rho, mc, nm)
% [lr, L, G, xh, njump] = cycle_jump(ev, Cf, lr, L, G, xh, cfg, sig, rho, mc, nm)
%
% Cycle-skip check of one record's nuisance linearization point (option 3
% of the cycle-skipping guards; run_pt_LR.m every mc.jump_every
% iterations during burn-in).
%
% A time shift tau_n(w) that is off by the same whole number of cycles at
% all frequencies (k/f) is harmless (the spectra are unchanged). The
% harmful kind is off by a cycle at the short periods but not at the long
% ones: Gauss-Newton started from the wrong side of half a cycle at the
% short periods locks onto it. For every branch n and both signs, this
% tries
%       tau_n(w) -> tau_n(w) +- T(w) r(w),
% r(w) = 0 at the longest and 1 at the shortest period of the band (linear
% in log f), i.e. one cycle at the short end, projected onto the tau basis
% (cfg.Psit); one Gauss-Newton step (nuis_relin) lets B_n, the other
% branches and dth adjust, and the change is kept if the record's
% log-likelihood increases by more than mc.jump_tol. Accepted changes are
% applied one after the other (greedy).
%
% lr, L, G, xh : current linearization point, log-likelihood, operator and
%                posterior mean of the record (as from nuis_relin)
% njump        : number of accepted changes
f = cfg.w(:)/(2*pi); T = 1./f;
r = (log(f) - log(min(f))) / (log(max(f)) - log(min(f)));
Kt = size(lr.tau0, 2);
Ps = cfg.Psit(:, 1:Kt);
dc = (Ps \ (T .* r)).';                         % 1 x Kt coefficients of one short-period cycle
tol = 0; if isfield(mc, 'jump_tol'), tol = mc.jump_tol; end
njump = 0;
for n = 1:nm
    if all(isnan(ev.cref(:, n))), continue; end
    for sgn = [1 -1]
        lt = lr; lt.tau0(n, :) = lt.tau0(n, :) + sgn*dc;
        lt = nuis_relin(ev, Cf, lt, cfg, sig, rho, mc, nm, true);
        [~, Lt, Gt, xt] = nuis_relin(ev, Cf, lt, cfg, sig, rho, mc, nm, false);
        if Lt > L + tol
            lr = lt; L = Lt; G = Gt; xh = xt; njump = njump + 1;
            break                                   % (the other sign undoes it)
        end
    end
end
end
