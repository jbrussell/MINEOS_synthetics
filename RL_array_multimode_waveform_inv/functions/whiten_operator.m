function [W, wld] = whiten_operator(first, jst, models)
% [W, wld] = whiten_operator(first, jst, models)
%
% Sparse whitening operator of a record whose noise in each contiguous
% segment (first = true at segment starts) is a stationary AR(p) process
% with unit variance (the noise std sigma, and per-sample levels, are
% applied separately in event_loglike_ar.m):
%       y_i = (x_i - sum_m phi_k(m) x_(i-m)) / sqrt(v_k),  k = min(q, p),
% q = position of sample i in its segment (0-based); this is the exact
% Cholesky whitening of the stationary process (the first samples of a
% segment use the lower-order models of the Levinson recursion).
%
% first  : [n x 1] logical, segment starts
% jst    : [n x 1] model index of each sample (e.g. its station), into models
% models : struct array with fields phi, v (levinson_ar.m)
% W      : [n x n] sparse, lower triangular, banded
% wld    : log det W = -0.5 * sum(log v_k) (added to the log-likelihood)
n = numel(first); first = first(:); jst = jst(:);
seg = cumsum(first);
fs = find(first); pos = (1:n)' - fs(seg);      % position within the segment (0-based)
pmax = max(arrayfun(@(M) numel(M.v) - 1, models));
I = zeros(n*(pmax + 1), 1); J = I; V = I; c = 0; wld = 0;
for i = 1:n
    M = models(jst(i));
    k = min(pos(i), numel(M.v) - 1);
    sv = sqrt(M.v(k+1));
    c = c + 1; I(c) = i; J(c) = i; V(c) = 1/sv;
    for m = 1:k
        c = c + 1; I(c) = i; J(c) = i - m; V(c) = -M.phi{k+1}(m)/sv;
    end
    wld = wld - log(sv);
end
W = sparse(I(1:c), J(1:c), V(1:c), n, n);
end
