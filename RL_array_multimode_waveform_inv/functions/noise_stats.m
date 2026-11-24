function [sd, ac, ps] = noise_stats(x, dt, fc, maxlag, fg)
% [sd, ac, ps] = noise_stats(x, dt, fc, maxlag, fg)
% Statistics of a noise segment x (sample interval dt) band-passed as the
% data (Tukey 0.1 taper, one-sided band-pass taper on fc = [f1 f2 f3 f4]),
% on the central 80% of the samples:
%   sd : standard deviation
%   ac : autocorrelation at lags 0..maxlag samples (biased estimator; NaN
%        beyond the available lags)
%   ps : two-sided power spectral density (units^2/Hz; periodogram, Tukey
%        0.2 taper, 3-bin running mean) at the frequencies fg (Hz)
% Used by read_event_data.m (inversion) and A4_plot_event_data.m (figures).
x = x(:) - mean(x); n = numel(x);
f = (0:n-1)'/(n*dt); ib = find(f > fc(1) & f < fc(4));
X = fft(x .* tukeywin_(n, 0.1)); Y = zeros(n, 1); Y(ib) = X(ib) .* bp_taper(f(ib), fc);
u = 2*real(ifft(Y)); k = round(0.1*n) + 1 : round(0.9*n);
sd = std(u(k));
uk = u(k) - mean(u(k)); nk = numel(uk); ac = NaN(1, maxlag + 1);
for L = 0:min(maxlag, nk - 1), ac(L+1) = sum(uk(1:nk-L) .* uk(1+L:nk)) / nk; end
ac = ac / ac(1);
w = tukeywin_(nk, 0.2); P = abs(fft(uk .* w)).^2 * dt / sum(w.^2);   % two-sided PSD estimate
P = conv(P, ones(3, 1)/3, 'same'); fk = (0:nk-1)'/(nk*dt);
h = fk < 1/(2*dt);
ps = interp1(fk(h), P(h), fg(:)', 'linear', NaN);
end
