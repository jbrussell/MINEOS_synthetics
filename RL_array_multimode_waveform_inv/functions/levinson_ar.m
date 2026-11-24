function [phi, v] = levinson_ar(acf, p, nugget)
% [phi, v] = levinson_ar(acf, p, nugget)
%
% Levinson-Durbin recursion: the stationary AR(p) process whose
% autocorrelation matches acf(1:p+1) = rho(0..p) exactly (the
% maximum-entropy extension of the measured autocorrelation).
%
% acf    : autocorrelation at lags 0..>=p (acf(1) = rho(0) > 0; normalized
%          here)
% p      : order
% nugget : white-noise fraction added to rho(0) before the fit (default
%          0.05). Band-passed noise has almost no power outside the band,
%          which makes high-order fits nearly singular (tiny innovation
%          variances, unstable whitening); the nugget keeps them stable.
%
% phi{k} : 1 x k prediction coefficients of the order-k model (k = 0..p,
%          phi{1} = [] for k = 0): x_i ~ sum_m phi{k+1}(m) x_(i-m)
% v(k+1) : innovation variance of the order-k model relative to the
%          process variance (v(1) = 1)
% The order is lowered if a reflection coefficient reaches |kappa| >= 0.999
% (an acf that is not positive definite): phi and v then stop there.
if nargin < 3 || isempty(nugget), nugget = 0.05; end
r = acf(:)' / acf(1);
r(1) = 1 + nugget; r = r / r(1);
phi = cell(1, p + 1); v = zeros(1, p + 1);
phi{1} = []; v(1) = 1;
a = [];
for k = 1:p
    kap = (r(k+1) - sum(a .* r(k:-1:2))) / v(k);
    if abs(kap) >= 0.999 || ~isfinite(kap)       % not positive definite: stop here
        phi = phi(1:k); v = v(1:k); return
    end
    a = [a - kap*fliplr(a), kap];
    phi{k+1} = a; v(k+1) = v(k)*(1 - kap^2);
end
end
