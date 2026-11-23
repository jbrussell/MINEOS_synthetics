function w = tukeywin(N, r)
% w = tukeywin(N, r)  (fallback when the signal toolbox/package is missing)
% Tukey (tapered cosine) window of N points, taper fraction r (default 0.5),
% as in MATLAB's tukeywin.
if nargin < 2, r = 0.5; end
x = linspace(0, 1, N)'; w = ones(N, 1);
if r <= 0, return; end
if r >= 1, w = 0.5*(1 - cos(2*pi*x)); return; end
i1 = x < r/2; i3 = x >= 1 - r/2;
w(i1) = 0.5*(1 + cos(2*pi/r*(x(i1) - r/2)));
w(i3) = 0.5*(1 + cos(2*pi/r*(x(i3) - 1 + r/2)));
end
