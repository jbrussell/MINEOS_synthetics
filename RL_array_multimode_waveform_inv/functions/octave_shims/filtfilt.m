function y = filtfilt(b, a, x)
% y = filtfilt(b, a, x)  (fallback when the signal toolbox/package is
% missing) Zero-phase forward-backward filtering along the first
% non-singleton dimension, with odd reflection padding of 3*(order) samples
% and steady-state initial conditions, as in MATLAB's filtfilt.
b = b(:).'; a = a(:).';
n = max(numel(a), numel(b));
b = [b, zeros(1, n - numel(b))] / a(1); a = [a, zeros(1, n - numel(a))] / a(1);
row = isrow(x); if row, x = x(:); end
nf = 3*(n - 1);
if size(x, 1) <= nf, error('filtfilt: data must be longer than %d samples', nf); end
if n > 1
    M = eye(n - 1) - [-a(2:n).', [eye(n - 2); zeros(1, n - 2)]];
    zi = M \ (b(2:n).' - a(2:n).'*b(1));
else
    zi = zeros(0, 1);
end
y = zeros(size(x));
for c = 1:size(x, 2)
    xc = x(:, c);
    xp = [2*xc(1) - xc(nf+1:-1:2); xc; 2*xc(end) - xc(end-1:-1:end-nf)];
    v = filter(b, a, xp, zi*xp(1));
    v = flipud(v);
    v = filter(b, a, v, zi*v(1));
    v = flipud(v);
    y(:, c) = v(nf+1:end-nf);
end
if row, y = y.'; end
end
