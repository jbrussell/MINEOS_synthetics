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
% (default [1 1]; 0 drops a wave type).
%
% CLd, CRd : [numel(Tdisp) x nm] data phase velocities (km/s)
% sig      : relative standard deviation (e.g. 0.0025 = 0.25%)
% nit, burn, thin, seed : iterations, burn-in, thinning, random seed
% out      : vsv, vsh, xi [numel(zplot) x nsave], kv, kx, L (trace),
%            CL, CR [numel(Tdisp) x nm x nsave], acc (accepted fraction)
if nargin < 14 || isempty(w), w = [1 1]; end
rng(seed);
need = cfg.need;
okL = need & isfinite(CLd); okR = need & isfinite(CRd);
loglik = @(S) -0.5*(w(1)*sum(((S.CL(okL) - CLd(okL)) ./ (sig*CLd(okL))).^2) + ...
                    w(2)*sum(((S.CR(okR) - CRd(okR)) ./ (sig*CRd(okR))).^2));
while true
    m = lr_common('draw_prior', cfg, mc);
    [S, ok] = lr_common('disp', struct(), m, cfg, Tdisp, fb, nm, 'both');
    if ok, break; end
end
L = loglik(S);
nsave = floor((nit - burn)/thin); k = 0;
out.vsv = zeros(numel(zplot), nsave); out.vsh = out.vsv; out.xi = out.vsv;
out.kv = zeros(nsave, 1); out.kx = out.kv; out.L = zeros(nit, 1);
out.CL = zeros(numel(Tdisp), nm, nsave); out.CR = out.CL; nacc = 0;
for it = 1:nit
    [mn, okp] = lr_common('propose', S.m, cfg, mc, 1);
    if okp
        [Sn, ok] = lr_common('disp', S, mn, cfg, Tdisp, fb, nm, 'both');
        if ok
            Ln = loglik(Sn);
            if log(rand) < Ln - L, S = Sn; L = Ln; nacc = nacc + 1; end
        end
    end
    out.L(it) = L;
    if it > burn && mod(it - burn, thin) == 0
        k = k + 1;
        v = lr_common('profiles', S, zplot, cfg);
        out.vsv(:, k) = v(:, 1); out.vsh(:, k) = v(:, 2); out.xi(:, k) = v(:, 3);
        out.kv(k) = numel(S.m.V.zc); out.kx(k) = numel(S.m.X.zc);
        out.CL(:, :, k) = S.CL; out.CR(:, :, k) = S.CR;
    end
end
out.acc = nacc/nit;
end
