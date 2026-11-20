function [Ltot, chi, lin] = fit_nuisance_LR(rec, Cf, sig, rho, cfg, mc, nm, ngn)
% [Ltot, chi, lin] = fit_nuisance_LR(rec, Cf, sig, rho, cfg, mc, nm, ngn)
%
% Fit the nuisance parameters of every record (time shifts, frequency-
% dependent arrival angles and the smooth spectral amplitudes B_n(w)) for
% a FIXED structure, given as band dispersion Cf = {Love, Rayleigh}
% [Nf x nm], by damped Gauss-Newton (spec_refine.m), and return the
% linearized-marginal log-likelihood (as in run_pt_LR) and the reduced
% chi-square of each record (residual of the best fit over sigma).
%
% Used for the log-likelihood of the true structure (idagrn data have no
% "true" nuisance parameters) and to check how well a dispersion model can
% fit the data at all.
if nargin < 8, ngn = 12; end
nrec = numel(rec);
lin.tau0 = zeros(nrec, nm); lin.dth0 = zeros(nrec, nm); lin.b0 = ones(nrec, nm);
lin = spec_refine(lin, rec, Cf, sig, rho, cfg, mc, mc.spec_K, mc.spec_Ka, mc.dthp, ngn);
cfgs = cfg; [cfgs.Phi, cfgs.Psi] = spec_basis(cfg.w, lin.K, lin.Ka);
Ltot = 0; chi = zeros(nrec, 1);
for r = 1:nrec
    sgr = sig(2*rec(r).wave-1:2*rec(r).wave);
    le = struct('tau0', lin.tau0(r,:), 'dth0', reshape(lin.dth0(r,:,:), nm, lin.Ka), ...
                'b0', reshape(lin.b0(r,:,:), nm, lin.K));
    G = event_G_spec(rec(r), Cf{rec(r).wave}, le, cfgs, sgr, mc.sb);
    [mul, sbl] = spec_layout(le, nm, lin.K, lin.Ka, mc.sb, mc.tau_prior, lin.dth_prior);
    [L, xh] = event_loglike_ar(rec(r).d, G, rec(r).wt, rec(r).first, sig, rho, mul, sbl);
    Ltot = Ltot + L;
    res = rec(r).d - G*xh;
    chi(r) = mean((res ./ sig(rec(r).wt)').^2);
end
end
