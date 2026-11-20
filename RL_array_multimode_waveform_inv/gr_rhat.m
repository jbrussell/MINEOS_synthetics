function R = gr_rhat(X)
% GR_RHAT  Gelman-Rubin potential scale reduction factor; X is n x m (m chains).
[n, m] = size(X);
W = mean(var(X, 0, 1));
B = n * var(mean(X, 1), 0);
if W < 1e-12, R = NaN; return; end
R = sqrt(((n-1)/n*W + B/n) / W);
end
