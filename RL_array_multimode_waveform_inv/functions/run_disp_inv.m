function out = run_disp_inv(CLd, CRd, sig, cfg, mc, Tdisp, fb, nm, zplot, nit, burn, thin, seed, w)
% out = run_disp_inv(CLd, CRd, sig, cfg, mc, Tdisp, fb, nm, zplot, nit, burn, thin, seed, w)
%
% Dispersion-only rjMcMC (diagnostic): the same model space, prior,
% proposals and MINEOS linearization as the waveform inversion
% (lr_common.m: B-spline Vsv and xi, kx <= kv, bottom tied to the
% reference), with a Gaussian likelihood on phase velocities instead of
% the waveforms:
%       log L = -1/2 sum ((C(m) - C_data) / (sig * C_data))^2
% over the Love and Rayleigh branches 0..nm-1 at the periods Tdisp where
% the reference has the branch (cfg.need), weighted by w = [wL wR]
% (default [1 1]; 0 drops a wave type). The optional priors of mc
% (smoothness, Vsv toward the reference, xi toward xi_ref_mean;
% lr_common 'logprior') enter the acceptance ratio.
%
% CLd, CRd : [numel(Tdisp) x nm] data phase velocities (km/s)
% sig      : relative standard deviation (e.g. 0.0025 = 0.25%), or a cell
%            {sigL, sigR} of [numel(Tdisp) x nm] (or [numel(Tdisp) x 1])
%            relative standard deviations per period (and branch)
% nit, burn, thin, seed : iterations, burn-in, thinning, random seed
% out      : vsv, vsh, xi [numel(zplot) x nsave], kv, kx, L (trace),
%            CL, CR [numel(Tdisp) x nm x nsave], acc (accepted fraction)
if nargin < 14 || isempty(w), w = [1 1]; end
rng(seed);
need = cfg.need;
okL = need & isfinite(CLd); okR = need & isfinite(CRd);
if iscell(sig), sL = sig{1} .* ones(size(CLd)); sR = sig{2} .* ones(size(CRd));
else, sL = sig * ones(size(CLd)); sR = sig * ones(size(CRd)); end
okL = okL & isfinite(sL) & sL > 0; okR = okR & isfinite(sR) & sR > 0;
loglik = @(S) -0.5*(w(1)*sum(((S.CL(okL) - CLd(okL)) ./ (sL(okL).*CLd(okL))).^2) + ...
                    w(2)*sum(((S.CR(okR) - CRd(okR)) ./ (sR(okR).*CRd(okR))).^2));
while true
    m = lr_common('draw_prior', cfg, mc);
    [S, ok] = lr_common('disp', struct(), m, cfg, Tdisp, fb, nm, 'both');
    if ok, break; end
end
L = loglik(S); lp = lr_common('logprior', S.m, cfg, mc);   % smoothness / Vsv / xi priors (mc; 0 if none)
nsave = floor((nit - burn)/thin); k = 0;
out.vsv = zeros(numel(zplot), nsave); out.vsh = out.vsv; out.xi = out.vsv;
out.kv = zeros(nsave, 1); out.kx = out.kv; out.L = zeros(nit, 1);
out.CL = zeros(numel(Tdisp), nm, nsave); out.CR = out.CL; nacc = 0; out.xs = NaN(nsave, 1); out.ss = NaN(nsave, 2);
for it = 1:nit
    [mn, okp] = lr_common('propose', S.m, cfg, mc, 1);
    if okp
        [Sn, ok] = lr_common('disp', S, mn, cfg, Tdisp, fb, nm, 'both');
        if ok
            Ln = loglik(Sn); lpn = lr_common('logprior', mn, cfg, mc);
            if log(rand) < Ln - L + lpn - lp, S = Sn; L = Ln; lp = lpn; nacc = nacc + 1; end
        end
    end
    if ((isfield(mc, 'xi_ref_hier') && mc.xi_ref_hier) || (isfield(mc, 'smooth_hier') && mc.smooth_hier)), [S.m, lp] = lr_common('hyper', S.m, cfg, mc); end
    out.L(it) = L;
    if it > burn && mod(it - burn, thin) == 0
        k = k + 1;
        v = lr_common('profiles', S, zplot, cfg);
        out.vsv(:, k) = v(:, 1); out.vsh(:, k) = v(:, 2); out.xi(:, k) = v(:, 3);
        out.kv(k) = numel(S.m.V.zc); out.kx(k) = numel(S.m.X.zc);
        if isfield(S.m, 'xs'), out.xs(k) = S.m.xs; end
        if isfield(S.m, 'ss'), out.ss(k, :) = S.m.ss; end
        out.CL(:, :, k) = S.CL; out.CR(:, :, k) = S.CR;
    end
end
out.acc = nacc/nit;
end
