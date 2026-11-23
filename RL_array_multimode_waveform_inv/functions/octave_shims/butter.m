function [b, a] = butter(n, Wn, ftype)
% [b, a] = butter(n, Wn, ftype)  (fallback when the signal toolbox/package
% is missing) Digital Butterworth filter of order n: low-pass (scalar Wn),
% band-pass (Wn = [w1 w2]) or high-pass (ftype 'high'); Wn normalized to the
% Nyquist frequency. Analog prototype, frequency transformation and the
% bilinear transform with pre-warping, as in MATLAB's butter.
if nargin < 3, ftype = ''; end
p0 = exp(1i*pi*(2*(1:n) + n - 1)/(2*n)).';     % unit-cutoff prototype poles
W = tan(pi*Wn/2);                              % pre-warped (bilinear with s = (z-1)/(z+1))
if numel(Wn) == 2                              % band-pass
    Bw = W(2) - W(1); W0 = sqrt(W(1)*W(2));
    pa = [p0*Bw/2 + sqrt((p0*Bw/2).^2 - W0^2); p0*Bw/2 - sqrt((p0*Bw/2).^2 - W0^2)];
    za = zeros(n, 1); ka = Bw^n;
elseif strcmpi(ftype, 'high')
    pa = W ./ p0; za = zeros(n, 1); ka = 1;
else                                           % low-pass
    pa = W * p0; za = zeros(0, 1); ka = W^n;
end
pd = (1 + pa) ./ (1 - pa);
zd = [(1 + za) ./ (1 - za); -ones(numel(pa) - numel(za), 1)];
kd = ka * real(prod(1 - za) / prod(1 - pa));
b = real(kd * poly(zd)); a = real(poly(pd));
end
