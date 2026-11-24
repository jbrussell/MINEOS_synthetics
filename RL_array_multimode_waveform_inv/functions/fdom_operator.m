function Wst = fdom_operator(first, wt, jst, dt, fband, S, taper_frac)
% Wst = fdom_operator(first, wt, jst, dt, fband, S, taper_frac)
%
% Frequency-domain "whitening" operator of a record (opt.noise_model =
% 'fdom'): in every contiguous segment (one station's window; first = true
% at segment starts) the residual is tapered (Tukey, taper_frac, default
% 0.2) and transformed to its DFT coefficients at the in-band frequencies
% f_k = k/(n dt), fband(1) < f_k < fband(2); the real and imaginary parts
% are scaled by 1/sqrt(c_k), with
%       c_k = 0.5 * sum(w.^2) * S_j(f_k) / dt
% the variance of each part for stationary noise of two-sided power
% spectral density S_j (station of the segment). For stationary noise the
% transformed coefficients have unit variance and are nearly independent
% across frequency (their correlation comes from the taper).
%
% first, wt, jst : [n x 1] segment starts, window type, station of each sample
% dt     : sample interval (s)
% fband  : [fmin fmax] (Hz)
% S      : function handle S(j, f): PSD of station j at frequencies f (data
%          units^2 per Hz, two-sided)
% Wst    : struct with M (sparse [m x n] operator), wt ([m x 1] window type
%          of each row), ld (log det term -0.5*sum(log c_k), added to the
%          log-likelihood) and the row frequencies f
if nargin < 7 || isempty(taper_frac), taper_frac = 0.2; end
first = first(:); n = numel(first);
st = find(first); en = [st(2:end) - 1; n];
I = []; J = []; V = []; wtr = []; fr = []; ld = 0; m = 0;
for s = 1:numel(st)
    ii = (st(s):en(s))'; ns = numel(ii); j = jst(ii(1));
    w = tukeywin_(ns, taper_frac);
    k = find((0:ns-1)'/(ns*dt) > fband(1) & (0:ns-1)'/(ns*dt) < fband(2) & (0:ns-1)' < ns/2) - 1;
    if isempty(k), continue; end
    f = k/(ns*dt);
    c = 0.5 * sum(w.^2) * S(j, f) / dt;               % variance of the real / imaginary part
    t = (0:ns-1);
    E = exp(-2i*pi*k*t/ns) .* w';                      % [nk x ns]
    R = [real(E) ./ sqrt(c); imag(E) ./ sqrt(c)];      % [2 nk x ns]
    [rr, cc] = find(R ~= 0);
    I = [I; m + rr]; J = [J; ii(cc)]; V = [V; R(sub2ind(size(R), rr, cc))]; %#ok<AGROW>
    m = m + 2*numel(k);
    wtr = [wtr; wt(ii(1))*ones(2*numel(k), 1)]; fr = [fr; f; f]; %#ok<AGROW>
    ld = ld - sum(log(c));                             % -0.5*log(c) for each of the 2 parts
end
Wst = struct('M', sparse(I, J, V, m, n), 'wt', wtr, 'ld', ld, 'f', fr);
end
