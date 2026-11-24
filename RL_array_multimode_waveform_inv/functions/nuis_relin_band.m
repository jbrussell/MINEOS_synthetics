function lr = nuis_relin_band(ev, Cf, lr, cfg, sig, rho, mc, nm, LP)
% lr = nuis_relin_band(ev, Cf, lr, cfg, sig, rho, mc, nm, LP)
%
% One Gauss-Newton step of the nuisance linearization point (as
% nuis_relin.m with update = true) on LOW-PASSED data: the windowed data
% ev.d and every column of the linearized operator G (event_G_spec.m) are
% put back on their [N x Nst] time grid (zero outside the windows),
% filtered with the zero-phase low-pass LP ([N x 1], on the FFT frequency
% grid) and windowed again, so data and predictions see the same filter
% and the same window edges. Used for the frequency continuation of
% spec_refine.m (long periods first, against cycle skipping); the noise
% model (sig, rho) is that of the unfiltered data, which only weights the
% Gauss-Newton step.
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); nb = 2*nm*K;
G = event_G_spec(ev, Cf, lr, cfg, sig(2*ev.wave-1:2*ev.wave), mc.sb);
Nst = size(cfg.xy, 1);
F = @(X) filt_(X, ev.idx, cfg.N, Nst, LP);
G = F(G); d = F(ev.d);
[mu, sb2, it, ith] = spec_layout(lr, nm, K, Ka, mc.sb, mc.taup, mc.dthp);
[~, xh] = event_loglike_ar(d, G, ev.wt, ev.first, sig, rho, mu, sb2);
lr.b0   = reshape(xh(1:2:nb) + 1i*xh(2:2:nb), K, nm).';
lr.tau0 = lr.tau0 + max(min(reshape(xh(it), size(it)), mc.max_dtau), -mc.max_dtau);
lr.dth0 = lr.dth0 + max(min(reshape(xh(ith), size(ith)), mc.max_dth), -mc.max_dth);
end

function Y = filt_(X, idx, N, Nst, LP)
% columns of X: samples at idx of [N x Nst] traces -> low-passed, same samples
Y = zeros(size(X));
Z = zeros(N, Nst);
for c = 1:size(X, 2)
    Z(:) = 0; Z(idx) = X(:, c);
    W = real(ifft(fft(Z) .* LP));
    Y(:, c) = W(idx);
end
end
