function [rec, pc] = path_precorrect(rec, cfg, nm, opt_pc)
% [rec, pc] = path_precorrect(rec, cfg, nm, opt_pc)
%
% Correct the reference modal spectra Uref for the source-to-array path,
% using the data. With a reference model that differs from the true path
% structure (here PREM vs PA5), the predicted fundamental-mode Rayleigh
% phase at the array is off by up to ~100 s at 16 s over 9500 km (many
% cycles), far more than the nuisance model of event_G_spec (a time shift
% and a 3-node smooth spectral correction per mode) can absorb. The other
% path delays (Love, overtones) are below ~30 s and are left to the
% nuisance model.
%
% Fundamental mode, for the wave types in opt_pc.waves (default 2 =
% Rayleigh): the data and the reference prediction are windowed with the
% fundamental-mode window and cross-correlated, stacked over the stations,
%     C(w) = sum_j D_j(w) conj(P_j(w)),   A(w) = sum_j |P_j(w)|^2,
% and smoothed over 5 frequencies. The phase of C is unwrapped by tracking
% from the longest period with a stacked coherence above 0.8 (where the
% delay is a fraction of a period) toward short periods, choosing at each
% frequency the 2pi branch closest to a linear extrapolation of the phase.
% A smooth delay tau(w) and log amplitude ratio log(|C|/A) (opt_pc.Kt
% Legendre terms in normalized log w, weighted by the coherence) are fitted
% and folded into Uref_0. This is the array-average (array-center)
% correction: the moveout across the array, which constrains the local
% structure, is not changed.
%
% Quality control: where the coherence is low, the unwrapping can skip
% cycles. A record is flagged (pc.ok(r) = false) when the unwrapped phase
% jumps by more than opt_pc.maxstep (default 0.8 rad) between neighbouring
% frequencies inside the pass band; the caller drops flagged records.
%
% Tests with PA5 data and a PREM reference (Rayleigh, 15 events): the
% fitted delay is within 1-5 s of the true path delay (up to 107 s) for
% the records that pass; the 3 flagged records (low signal-to-noise) were
% off by 2-100 s. Love and overtone corrections estimated the same way were
% less accurate than no correction (the Love fundamental and overtone
% windows overlap, and branches above 2 contaminate the overtone window).
%
% pc.tau{r}      [Nf x nm] delay (s) applied to each mode of record r
% pc.tau_meas{r} [Nf x 1] measured (unwrapped) fundamental-mode delay
% pc.coh{r}      [Nf x 1] stacked coherence of the fundamental window
% pc.amp{r}      [Nf x 1] amplitude factor applied to mode 0
% pc.ok(r)       false if the record failed the unwrapping check
%
% Bulk delay, for the wave types in opt_pc.waves_bulk (default 1 = Love):
% the Love path delays are smaller (< ~20 s here) and vary little with
% frequency, but they are large enough to put the per-mode time shifts of
% the inversion on the wrong cycle. A single delay for all branches, the
% lag (within +-opt_pc.tmax s) of the peak of the stacked cross-correlation
% of the data with the total reference prediction (both windows), is
% folded into Uref.
%
% jbrussell 11/2026

if nargin < 4, opt_pc = struct(); end
Kt = optget(opt_pc, 'Kt', 4);
waves = optget(opt_pc, 'waves', 2);
waves_bulk = optget(opt_pc, 'waves_bulk', 1);
tmax = optget(opt_pc, 'tmax', 40);
lags = (-tmax:0.25:tmax)';
maxstep = optget(opt_pc, 'maxstep', 0.8);
w = cfg.w(:); Nf = numel(w); N = cfg.N; Nst = size(cfg.xy, 1);
x = log(w); xi = 2*(x - min(x))/(max(x) - min(x)) - 1;
Pl = [ones(Nf,1), xi, 1.5*xi.^2 - 0.5, 2.5*xi.^3 - 1.5*xi, (35*xi.^4 - 30*xi.^2 + 3)/8];
Pl = Pl(:, 1:Kt);
fb = w/(2*pi);
band = fb >= cfg.fc(2) & fb <= cfg.fc(3);
pc.ok = true(1, numel(rec));
for r = 1:numel(rec)
    ev = rec(r);
    pc.tau{r} = zeros(Nf, nm); pc.amp{r} = ones(Nf, 1);
    pc.tau_meas{r} = NaN(Nf, 1); pc.coh{r} = NaN(Nf, 1);
    if any(ev.wave == waves_bulk)
        D = zeros(N, Nst); D(ev.idx) = ev.d;
        X = zeros(N, Nst); for n = 1:nm, X(cfg.ib, :) = X(cfg.ib, :) + ev.Uref{n}; end
        P = 2*real(ifft(X));
        M = zeros(N, Nst); M(ev.idx) = 1;
        Dw = fft(D .* M); Pw = fft(P .* M);
        C = sum(Dw(cfg.ib, :) .* conj(Pw(cfg.ib, :)), 2);
        [~, il] = max(real(exp(1i*lags*w') * C));   % D = P exp(-i w tau): peak at lag = tau
        for n = 1:nm, rec(r).Uref{n} = ev.Uref{n} .* exp(-1i*w*lags(il)); end
        pc.tau{r}(:) = lags(il);
        continue
    end
    if ~any(ev.wave == waves), continue; end
    D = zeros(N, Nst); D(ev.idx) = ev.d;
    X = zeros(N, Nst); X(cfg.ib, :) = ev.Uref{1};
    P = 2*real(ifft(X));
    M = zeros(N, Nst); M(ev.idx(ev.wt == 2*ev.wave)) = 1;
    Dw = fft(D .* M); Pw = fft(P .* M); Dw = Dw(cfg.ib, :); Pw = Pw(cfg.ib, :);
    C = sum(Dw .* conj(Pw), 2);
    A = sum(abs(Pw).^2, 2);
    coh = abs(C) ./ sqrt(A .* sum(abs(Dw).^2, 2) + eps);
    C = conv(C, ones(5,1)/5, 'same');
    ph = angle(C);
    k0 = find(coh > 0.8, 1, 'first');
    if isempty(k0), [~, k0] = max(coh); end
    phu = ph;
    for k = k0+1:Nf                                   % toward short periods
        if k >= k0 + 2, pr = 2*phu(k-1) - phu(k-2); else, pr = phu(k-1); end
        phu(k) = ph(k) + 2*pi*round((pr - ph(k))/(2*pi));
    end
    for k = k0-1:-1:1                                 % toward long periods
        if k <= k0 - 2, pr = 2*phu(k+1) - phu(k+2); else, pr = phu(k+1); end
        phu(k) = ph(k) + 2*pi*round((pr - ph(k))/(2*pi));
    end
    st = abs(diff(phu)); st = st(band(2:end) & band(1:end-1));
    pc.ok(r) = max(st) <= maxstep;
    tm = -phu ./ w;
    Wt = coh.^2 .* w;                                 % phase precision
    a = (Pl .* Wt) \ (tm .* Wt);
    tau = Pl*a;
    la = log(abs(C) ./ A); ok = isfinite(la);
    b = (Pl(ok,1:3) .* coh(ok).^2) \ (la(ok) .* coh(ok).^2);
    amp = exp(max(min(Pl(:,1:3)*b, log(3)), -log(3)));
    rec(r).Uref{1} = ev.Uref{1} .* amp .* exp(-1i*w.*tau);
    pc.tau{r}(:, 1) = tau; pc.tau_meas{r} = tm; pc.coh{r} = coh; pc.amp{r} = amp;
end
end
