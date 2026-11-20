function G = event_G_lin(ev, Cf, lin, cfg, sig, sb)
% EVENT_G_LIN  Linearized operator for one event in which the per-mode time
% shifts and arrival-angle anomalies are linear parameters, so they can be
% marginalized analytically together with the complex amplitudes.
%
%   d ~ G * [Re b1; Im b1; ...; Re bN; Im bN; dtau_1..N; dth_1..N]
%
% The b-columns are the exact modal traces at the linearization point
% (tau0, dth0). The dtau/dth columns are the derivatives of the predicted
% traces at amplitudes b0 (Gauss-Newton linearization):
%   d/dtau_n : multiply mode-n spectrum by -i w
%   d/dth_n  : multiply by -i w / c_n * d(dDelta_j)/d(theta)
%
%   lin.tau0, lin.dth0 : 1 x nm (s, deg);  lin.b0 : 1 x nm complex
%   sig, sb (optional) : if given, b0 is re-estimated for this structure by
%   damped least squares on the b-columns (noise sig, prior std sb) before
%   the derivative columns are formed. This keeps the linearization
%   accurate for structures far from the pilot (e.g. in hot PT chains).

nm  = size(Cf, 2);
Nst = size(cfg.xy, 1);
w   = cfg.w(:);
M   = numel(ev.idx);
G   = zeros(M, 4*nm);
S = cell(nm, 1); cc = cell(nm, 1); dDa = cell(nm, 1);
for n = 1:nm
    az = ev.az_prop + lin.dth0(n);
    dD     = cfg.xy * [sind(az); cosd(az)];
    dDa{n} = cfg.xy * [cosd(az); -sind(az)] * pi/180;   % d(dD)/d(az in deg)
    c  = Cf(:, n);
    ok = ~isnan(c) & ~isnan(ev.T0(:, n));
    c(~ok) = 1; cc{n} = c;
    src = ev.A(:, n) .* ok .* exp(-1i * w .* (ev.T0(:, n) - ev.tred + lin.tau0(n)));
    S{n} = repmat(src, 1, Nst) .* exp(-1i * (w ./ c) * dD');   % Nf x Nst
    X = zeros(cfg.N, Nst); X(cfg.ib, :) = S{n};
    a = ifft(X);
    gt = 2*real(a); ht = -2*imag(a);
    G(:, 2*n-1) = gt(ev.idx);
    G(:, 2*n)   = ht(ev.idx);
end
b0 = lin.b0;
if nargin >= 6                                          % re-estimate b0
    Gb = G(:, 1:2*nm); lam = (mean(sig)/sb)^2;
    x = (Gb'*Gb + lam*eye(2*nm)) \ (Gb'*ev.d + lam*repmat([1; 0], nm, 1));
    b0 = (x(1:2:end) + 1i*x(2:2:end)).';
end
for n = 1:nm
    X = zeros(cfg.N, Nst);
    X(cfg.ib, :) = b0(n) * (-1i*w) .* S{n};                      % d/dtau
    at = 2*real(ifft(X));
    X(cfg.ib, :) = b0(n) * ((-1i*w ./ cc{n}) * dDa{n}') .* S{n}; % d/dtheta
    aa = 2*real(ifft(X));
    G(:, 2*nm + n) = at(ev.idx);
    G(:, 3*nm + n) = aa(ev.idx);
end
end
