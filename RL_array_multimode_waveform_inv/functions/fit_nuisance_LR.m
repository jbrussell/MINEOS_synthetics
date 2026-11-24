function [Ltot, chi, linc] = fit_nuisance_LR(rec, Cf, sig, rho, cfg, mc, nm, ngn)
% [Ltot, chi, linc] = fit_nuisance_LR(rec, Cf, sig, rho, cfg, mc, nm, ngn)
%
% Fit the nuisance parameters of every record (tau_n(w), dth_n(w) and the
% smooth spectral amplitudes B_n(w)) for a FIXED structure, given as band
% dispersion Cf = {Love, Rayleigh} [Nf x nm], by damped Gauss-Newton
% (spec_refine.m), and return the linearized-marginal log-likelihood (as
% in the samplers) and the reduced chi-square of each record (residual of
% the best fit over sigma).
%
% Used for the log-likelihood of the true structure (idagrn data have no
% "true" nuisance parameters), to initialize the samplers' linearization
% points, and to check how well a dispersion model can fit the data.
if nargin < 8, ngn = 30; end
nrec = numel(rec);
linc = spec_refine([], rec, Cf, sig, rho, cfg, mc, nm, ngn);
Ltot = 0; chi = zeros(nrec, 1);
for r = 1:nrec
    [~, L, G, xh] = nuis_relin(rec(r), Cf{rec(r).wave}, linc{r}, cfg, sig, rho, mc, nm, false);
    Ltot = Ltot + L;
    res = rec(r).d - G*xh;
    nl_ = rec_nlev(rec(r)); if isempty(nl_), nl_ = 1; end
    chi(r) = mean((res ./ (sig(rec(r).wt)' .* nl_)).^2);
end
end
