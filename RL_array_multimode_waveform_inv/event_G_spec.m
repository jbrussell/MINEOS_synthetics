function G = event_G_spec(ev, Cf, lin, cfg, sig, sb, basis_only)
% EVENT_G_SPEC  Linearized operator for one record with
%   (1) a smooth complex spectral correction B_n(w) per mode (source
%       excitation and path dispersion; common to all stations), and
%   (2) a FREQUENCY-DEPENDENT arrival-angle anomaly dth_n(w) per mode,
% applied to reference modal spectra Uref (e.g. idagrn6_mask branches):
%
%   U_jn(w) = B_n(w) Uref_jn(w) exp(-i w tau_n)
%             exp(-i w [dD_j(az + dth_n(w))/c_n(w;m) - dD_j(az)/cref_n(w)])
%
%   B_n(w)   = sum_k beta_nk Phi_k(w)     (Phi = cfg.Phi, Nf x K)
%   dth_n(w) = sum_a  th_na  Psi_a(w)     (Psi = cfg.Psi, Nf x Ka; see spec_basis)
%
% beta is exactly linear and marginalized; tau and th are linearized about
% lin (Gauss-Newton), as in event_G_lin. With K = Ka = 1 and Uref, cref
% from the analytic source term, this reproduces event_G_lin exactly.
%
% lin.tau0 : 1 x nm (s)    lin.dth0 : nm x Ka (deg)    lin.b0 : nm x K complex
% ev.Uref  : 1 x nm cell of Nf x Nst      ev.cref : Nf x nm
% Optional sig, sb: re-estimate beta for this structure by damped LS
% before forming the derivative columns (as event_G_lin). Pass [] to skip.
% Optional basis_only = true: return only the beta columns, i.e. the exact
% (non-linearized) modal waveform basis for the given tau0 and dth0; used by
% the ordinary sampler, which samples tau and dth directly.
%
% Columns (see spec_layout): [beta (2*nm*K) | dtau (nm) | dth (nm*Ka, mode-major)]

nm  = size(Cf, 2);
Nst = size(cfg.xy, 1);
w   = cfg.w(:);
Phi = cfg.Phi; Psi = cfg.Psi; K = size(Phi, 2); Ka = size(Psi, 2);
nb  = 2*nm*K;
G   = zeros(numel(ev.idx), nb + nm + nm*Ka);
x1 = cfg.xy(:,1)'; y1 = cfg.xy(:,2)';
dD0 = sind(ev.az_prop)*x1 + cosd(ev.az_prop)*y1;          % 1 x Nst
S = cell(nm, 1); cc = cell(nm, 1); dDa = cell(nm, 1);
for n = 1:nm
    az = ev.az_prop + Psi * lin.dth0(n, :).';              % Nf x 1 (deg)
    dD     = sind(az)*x1 + cosd(az)*y1;                     % Nf x Nst
    dDa{n} = (cosd(az)*x1 - sind(az)*y1) * pi/180;          % d dD / d az(deg)
    c  = Cf(:, n); c0 = ev.cref(:, n);
    ok = ~isnan(c) & ~isnan(c0);
    c(~ok) = 1; c0(~ok) = 1; cc{n} = c;
    dphi = (w ./ c) .* dD - (w ./ c0) * dD0;
    S{n} = ev.Uref{n} .* ok .* exp(-1i * w * lin.tau0(n)) .* exp(-1i * dphi);
    for k = 1:K
        X = zeros(cfg.N, Nst); X(cfg.ib, :) = Phi(:, k) .* S{n};
        a = ifft(X);
        col = 2*((n-1)*K + k) - 1;
        gt = 2*real(a); ht = -2*imag(a);
        G(:, col) = gt(ev.idx); G(:, col+1) = ht(ev.idx);
    end
end
if nargin >= 7 && basis_only, G = G(:, 1:nb); return; end
b0 = lin.b0;
if nargin >= 6 && ~isempty(sig)
    Gb = G(:, 1:nb); lam = (mean(sig)/sb)^2;
    x = (Gb'*Gb + lam*eye(nb)) \ (Gb'*ev.d + lam*repmat([1; 0], nm*K, 1));
    b0 = reshape(x(1:2:end) + 1i*x(2:2:end), K, nm).';
end
for n = 1:nm
    BS = (Phi * b0(n, :).') .* S{n};                        % predicted spectra
    X = zeros(cfg.N, Nst);
    X(cfg.ib, :) = (-1i*w) .* BS;                           % d/dtau_n
    at = 2*real(ifft(X));
    G(:, nb + n) = at(ev.idx);
    for a = 1:Ka                                            % d/dth_na
        X(cfg.ib, :) = (-1i*w ./ cc{n} .* Psi(:, a)) .* dDa{n} .* BS;
        aa = 2*real(ifft(X));
        G(:, nb + nm + (n-1)*Ka + a) = aa(ev.idx);
    end
end
end
