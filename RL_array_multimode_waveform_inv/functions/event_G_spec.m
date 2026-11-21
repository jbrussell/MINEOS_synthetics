function G = event_G_spec(ev, Cf, lin, cfg, sig, sb, basis_only)
% EVENT_G_SPEC  Linearized operator for one record with
%   (1) a smooth complex spectral correction B_n(w) per mode (source
%       excitation and path dispersion; common to all stations),
%   (2) a FREQUENCY-DEPENDENT time shift tau_n(w) per mode (path delay
%       relative to the reference spectra; common to all stations), and
%   (3) a FREQUENCY-DEPENDENT arrival-angle anomaly dth_n(w) per mode,
% applied to reference modal spectra Uref (idagrn6_mask branches):
%
%   U_jn(w) = B_n(w) Uref_jn(w) exp(-i w tau_n(w))
%             exp(-i w [dD_j(az + dth_n(w))/c_n(w;m) - dD_j(az)/cref_n(w)])
%
%   B_n(w)   = sum_k beta_nk Phi_k(w)     (Phi  = cfg.Phi,  Nf x K)
%   tau_n(w) = sum_a  ta_na Psit_a(w)     (Psit = cfg.Psit, Nf x Kt; see spec_basis)
%   dth_n(w) = sum_a  th_na  Psi_a(w)     (Psi  = cfg.Psi,  Nf x Ka)
%
% beta is exactly linear and marginalized; tau and th are linearized about
% lin (Gauss-Newton).
%
% lin.tau0 : nm x Kt (s; a 1 x nm row is a constant tau_n)
% lin.dth0 : nm x Ka (deg)    lin.b0 : nm x K complex
% ev.Uref  : 1 x nm cell of Nf x Nst      ev.cref : Nf x nm
% Optional sig, sb: re-estimate beta for this structure by damped LS
% before forming the derivative columns. Pass [] to skip.
% Optional basis_only = true: return only the beta columns, i.e. the exact
% (non-linearized) modal waveform basis for the given tau0 and dth0; used
% when tau and dth are sampled directly.
%
% Columns (see spec_layout):
%   [beta (2*nm*K) | dtau (nm*Kt, mode-major) | dth (nm*Ka, mode-major)]

nm  = size(Cf, 2);
Nst = size(cfg.xy, 1);
w   = cfg.w(:);
Phi = cfg.Phi; Psi = cfg.Psi; K = size(Phi, 2); Ka = size(Psi, 2);
T0 = lin.tau0;
if size(T0, 1) == 1 && size(T0, 2) == nm && nm > 1, T0 = T0(:); end   % constant tau_n
Kt = size(T0, 2);
if isfield(cfg, 'Psit'), Psit = cfg.Psit(:, 1:Kt); else, Psit = ones(numel(w), 1); end
nb  = 2*nm*K;
G   = zeros(numel(ev.idx), nb + nm*Kt + nm*Ka);
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
    S{n} = ev.Uref{n} .* ok .* exp(-1i * w .* (Psit * T0(n, :).')) .* exp(-1i * dphi);
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
    for a = 1:Kt                                            % d/dtau_na
        X(cfg.ib, :) = (-1i*w .* Psit(:, a)) .* BS;
        at = 2*real(ifft(X));
        G(:, nb + (n-1)*Kt + a) = at(ev.idx);
    end
    for a = 1:Ka                                            % d/dth_na
        X(cfg.ib, :) = (-1i*w ./ cc{n} .* Psi(:, a)) .* dDa{n} .* BS;
        aa = 2*real(ifft(X));
        G(:, nb + nm*Kt + (n-1)*Ka + a) = aa(ev.idx);
    end
end
end
