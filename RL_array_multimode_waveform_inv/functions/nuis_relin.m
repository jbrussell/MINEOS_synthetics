function [lr, L, G, xh] = nuis_relin(ev, Cf, lr, cfg, sig, rho, mc, nm, update)
% [lr, L, G, xh] = nuis_relin(ev, Cf, lr, cfg, sig, rho, mc, nm, update)
%
% Linearized (Laplace) marginal log-likelihood of one record for a fixed
% structure (band dispersion Cf, Nf x nm), with the nuisance parameters
% B_n(w) (exactly linear), tau_n(w) and dth_n(w) (linearized about lr)
% integrated out with Gaussian priors (mc.sb, mc.taup, mc.dthp), and
% optionally one Gauss-Newton update of the linearization point.
%
% lr : linearization point: tau0 (nm x Kt, s), dth0 (nm x Ka, deg),
%      b0 (nm x K complex)
% update = true: move lr to the posterior mean of tau, dth and B (steps
%      limited to mc.max_dtau s and mc.max_dth deg per update); L and G
%      are those at the OLD point, i.e. the likelihood that was evaluated.
% Returns the log-likelihood L, the operator G at lr (before the update)
% and the posterior mean xh of [beta; dtau; dth].
if nargin < 9, update = false; end
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); nb = 2*nm*K;
G = event_G_spec(ev, Cf, lr, cfg, sig(2*ev.wave-1:2*ev.wave), mc.sb);
[mu, sb2, it, ith] = spec_layout(lr, nm, K, Ka, mc.sb, mc.taup, mc.dthp);
[L, xh] = event_loglike_ar(ev.d, G, ev.wt, ev.first, sig, rho, mu, sb2);
if update
    lr.b0   = reshape(xh(1:2:nb) + 1i*xh(2:2:nb), K, nm).';
    lr.tau0 = lr.tau0 + max(min(reshape(xh(it), size(it)), mc.max_dtau), -mc.max_dtau);
    lr.dth0 = lr.dth0 + max(min(reshape(xh(ith), size(ith)), mc.max_dth), -mc.max_dth);
end
end
