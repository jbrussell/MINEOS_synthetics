function linc = spec_refine(linc, rec, Cf, sig, rho, cfg, mc, nm, ngn)
% linc = spec_refine(linc, rec, Cf, sig, rho, cfg, mc, nm, ngn)
% Refine the linearization points of the nuisance parameters (tau_n(w),
% dth_n(w), B_n(w)) of every record for a FIXED structure Cf = {Love,
% Rayleigh} band dispersion by damped Gauss-Newton steps (nuis_relin).
% linc : {nrec} structs tau0 (nm x Kt), dth0 (nm x Ka), b0 (nm x K); pass
%        [] to start from tau = dth = 0, B = 1.
% ngn  : number of Gauss-Newton steps (default 30).
%
% Continuation: when starting from [] with a frequency-dependent time shift
% (Kt > 1), the first ngn steps keep the higher tau_n(w) terms at zero
% (constant tau_n), and ngn more steps then free them. Freeing all terms
% at once lets the overlapping branches trade delays and Gauss-Newton ends
% in poorer local minima (with PA5 data and the step-1 reference: median
% chi-square of the Love fundamental-mode windows 2.0 vs 1.13 with the
% continuation, at the true structure).
%
% Frequency continuation (from [], mc.fcont_stages = ns > 1; default 1 = off,
% because in a test with PA5 data it fit worse: log L of the true structure
% 205.9k vs 215.4k. The Legendre tau_n(w) terms fitted in the long-period
% stages extrapolate freely to the short periods; a basis on the current
% band or per-frequency unwrapping would be needed):
% against cycle skipping, the fit starts at the long periods, where the
% path delays are a small fraction of a cycle, and the band is widened
% toward the short periods in ns stages (low-pass corners log-spaced over
% the band, the last stage the full band; nuis_relin_band.m), each stage
% starting from the time shifts tau_n(w) of the previous one. Without it,
% Gauss-Newton is dominated by the short periods (the most cycles) and can
% lock onto a time shift a cycle off there but not at the long periods
% (not degenerate, unlike a whole-cycle shift k/f). The first stage keeps
% the higher tau_n(w) terms at zero as above.
if nargin < 9 || isempty(ngn), ngn = 30; end
nrec = numel(rec);
Kt = size(cfg.Psit, 2);
stages = {mc}; LPs = {[]};
if isempty(linc)
    linc = nuis_zero(nrec, nm, cfg);
    mc1 = mc;
    if Kt > 1, mc1.taup = [mc.taup(1), 1e-3*ones(1, Kt - 1)]; end   % higher terms held at 0
    ns = 1; if isfield(mc, 'fcont_stages') && ~isempty(mc.fcont_stages), ns = mc.fcont_stages; end
    if ns > 1
        f = cfg.w(:)/(2*pi); f1 = min(f); f2 = max(f);
        fk = min((0:cfg.N-1)', cfg.N - (0:cfg.N-1)') / (cfg.N*cfg.dt);   % FFT frequency grid
        stages = {}; LPs = {};
        for j = 1:ns-1
            fc = f1*(f2/f1)^(j/ns);                                       % low-pass corner
            LP = cos(min(max((fk - fc)/(0.3*fc), 0), 1)*pi/2).^2;         % cos^2 taper to 1.3 fc
            stages{end+1} = mc; LPs{end+1} = LP; %#ok<AGROW>
        end
        stages{1} = mc1;
        stages{end+1} = mc; LPs{end+1} = [];                              % full band
    elseif Kt > 1
        stages = {mc1, mc}; LPs = {[], []};
    end
end
for k = 1:numel(stages)
    for r = 1:nrec
        for itr = 1:ngn
            if isempty(LPs{k})
                linc{r} = nuis_relin(rec(r), Cf{rec(r).wave}, linc{r}, cfg, sig, rho, stages{k}, nm, true);
            else
                linc{r} = nuis_relin_band(rec(r), Cf{rec(r).wave}, linc{r}, cfg, sig, rho, stages{k}, nm, LPs{k});
            end
        end
    end
end
end
