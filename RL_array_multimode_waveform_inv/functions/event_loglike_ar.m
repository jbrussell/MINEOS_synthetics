function [L, bhat] = event_loglike_ar(d, G, wt, first, sig, rho, mu_b, sb2)
% EVENT_LOGLIKE_AR  Marginal log-likelihood for one record with correlated
% noise (Bodin et al., 2012, exponential correlation law): within each
% contiguous window segment the noise covariance is
%
%       C_ij = sigma_w^2 * r_w^|i-j|          (an AR(1) process)
%
% with a noise level sigma_w and correlation r_w for each window type w
% (Love/Rayleigh x overtone/fundamental). The inverse of an AR(1)
% covariance is tridiagonal, so the data and every column of G are
% "whitened" segment by segment,
%
%       y_1 = x_1,  y_i = (x_i - r x_(i-1)) / sqrt(1 - r^2),
%
% after which the noise is white with variance sigma_w^2 and
% log det C = sum_w [ n_w log sigma_w^2 + (n_w - n_seg,w) log(1 - r_w^2) ].
% The complex modal amplitudes are then marginalized as in event_loglike.
%
%   first : logical, true at the first sample of every contiguous segment
%   sig, rho : 1 x 4 noise levels and correlations

rw = reshape(rho(wt), [], 1);
nf = ~first(:);
s  = sqrt(1 - rw(nf).^2);
ip = find(nf);
dw = d;  dw(ip) = (d(ip) - rw(ip).*d(ip-1)) ./ s;
Gw = G;  Gw(ip,:) = (G(ip,:) - rw(ip).*G(ip-1,:)) ./ s;
if nargout > 1
    [L, bhat] = event_loglike(dw, Gw, wt, sig, mu_b, sb2);
else
    L = event_loglike(dw, Gw, wt, sig, mu_b, sb2);
end
L = L - sum(log(s));                 % -0.5 * log det of the correlation part
end
