function w = tukeywin_(n, a)
% w = tukeywin_(n, a)  Tukey (tapered cosine) window of length n with taper
% fraction a (as tukeywin in the Signal Processing Toolbox).
x = (0:n-1)'/(n-1);
w = ones(n, 1);
i1 = x < a/2;      w(i1) = 0.5*(1 + cos(2*pi/a*(x(i1) - a/2)));
i2 = x >= 1 - a/2; w(i2) = 0.5*(1 + cos(2*pi/a*(x(i2) - 1 + a/2)));
end
