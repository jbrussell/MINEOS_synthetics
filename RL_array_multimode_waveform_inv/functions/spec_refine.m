function [linc, hist, pick] = spec_refine(linc, rec, Cf, sig, rho, cfg, mc, nm, ngn)
% [linc, hist, pick] = spec_refine(linc, rec, Cf, sig, rho, cfg, mc, nm, ngn)
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
%
% Band-restricted time-shift basis (mc.fcont_basis = 'banded'; default
% 'global' = the Legendre basis of the whole band in every stage): in a low-passed
% stage with corner fc, tau_n(w) is described by Legendre polynomials in
% normalized log f on [f_min, fc] only, held constant (its value at fc) at
% shorter periods, where the stage has no data. With the full-band basis the
% terms fitted at long periods extrapolate freely to the short periods and
% start the next stage cycles off there. Between stages tau_n(w) is carried
% over by least squares on the next stage's band; the last stage uses
% cfg.Psit.
%
% Test (PA5, matched 8 branches, nuisance fitted at the true dispersion; log L
% of the true structure): direct start 215367.5, continuation with the
% global basis 205886.7, with the banded basis 197717.2. The banded basis fixes the records where the full-band terms extrapolate
% (Rayleigh record 17: 7284 vs -6145, old 7316) but fails where the Love
% fundamental mode and first overtone overlap at long periods: the first
% stage already lands a cycle off or trades delays between them, and the
% banded basis carries that through (Love record 11: -1752; global 6689,
% direct 7257).
%
% Start (from [], mc.nuis_start):
%   'best' (default): the three starts below, and for every record the one
%          with the highest log-likelihood (records are independent, so the
%          choice is exact); pick (optional output) = 1, 2, 3 per record.
%          Test (PA5, matched 8 branches): at the true dispersion log L
%          220383 (direct 215368, global 205887, banded 197717; records won
%          8/16/6); at the reference dispersion (the chains' lin0) 219068
%          (direct 216177, global 217464, banded 189666; 9/13/8). One
%          cycle-jump pass (cycle_jump.m) then adds 817 and 275.
%   'direct'              : the whole band from the start: constant tau_n,
%                           then all terms
%   'continuation_global' : frequency continuation (mc.fcont_stages stages,
%                           default 4), tau basis over the whole band
%   'continuation_banded' : frequency continuation, tau basis restricted to
%                           each stage's band
% Without mc.nuis_start: mc.fcont_stages (> 1: continuation) and
% mc.fcont_basis decide, as before.
%
% hist (optional output): per stage, the linearization points (lin), the
% time-shift basis (Psit) and the low-pass (LP, [] for the full band).

if nargin < 9 || isempty(ngn), ngn = 30; end
nrec = numel(rec); pick = [];
if isempty(linc) && isfield(mc, 'nuis_start') && ~isempty(mc.nuis_start)
    ns_ = 4; if isfield(mc, 'fcont_stages') && ~isempty(mc.fcont_stages) && mc.fcont_stages > 1, ns_ = mc.fcont_stages; end
    switch mc.nuis_start
        case 'best'
            nm_ = {'direct', 'continuation_global', 'continuation_banded'}; Ls = zeros(nrec, 3); Lc = cell(1, 3);
            for v = 1:3
                m_ = mc; m_.nuis_start = nm_{v};
                Lc{v} = spec_refine([], rec, Cf, sig, rho, cfg, m_, nm, ngn);
                for r = 1:nrec
                    [~, Ls(r, v)] = nuis_relin(rec(r), Cf{rec(r).wave}, Lc{v}{r}, cfg, sig, rho, mc, nm, false);
                end
            end
            [~, pick] = max(Ls, [], 2);
            linc = cell(nrec, 1);
            for r = 1:nrec, linc{r} = Lc{pick(r)}{r}; end
            hist = struct('lin', {}, 'Psit', {}, 'LP', {});
            fprintf('nuisance start, best per record: direct %d, continuation_global %d, continuation_banded %d records; log L %.1f (direct %.1f, global %.1f, banded %.1f)\n', ...
                sum(pick == 1), sum(pick == 2), sum(pick == 3), sum(max(Ls, [], 2)), sum(Ls(:, 1)), sum(Ls(:, 2)), sum(Ls(:, 3)));
            return
        case 'direct'
            mc.fcont_stages = 1;
        case 'continuation_global'
            mc.fcont_stages = ns_; mc.fcont_basis = 'global';
        case 'continuation_banded'
            mc.fcont_stages = ns_; mc.fcont_basis = 'banded';
        otherwise
            error('spec_refine: unknown mc.nuis_start ''%s''', mc.nuis_start);
    end
end
Kt = size(cfg.Psit, 2);
stages = {mc}; LPs = {[]}; Pst = {cfg.Psit};
if isempty(linc)
    linc = nuis_zero(nrec, nm, cfg);
    mc1 = mc;
    if Kt > 1, mc1.taup = [mc.taup(1), 1e-3*ones(1, Kt - 1)]; end   % higher terms held at 0
    ns = 1; if isfield(mc, 'fcont_stages') && ~isempty(mc.fcont_stages), ns = mc.fcont_stages; end
    bb = isfield(mc, 'fcont_basis') && strcmp(mc.fcont_basis, 'banded');
    if ns > 1
        f = cfg.w(:)/(2*pi); f1 = min(f); f2 = max(f);
        fk = min((0:cfg.N-1)', cfg.N - (0:cfg.N-1)') / (cfg.N*cfg.dt);   % FFT frequency grid
        stages = {}; LPs = {}; Pst = {};
        for j = 1:ns-1
            fc = f1*(f2/f1)^(j/ns);                                       % low-pass corner
            LP = cos(min(max((fk - fc)/(0.3*fc), 0), 1)*pi/2).^2;         % cos^2 taper to 1.3 fc
            stages{end+1} = mc; LPs{end+1} = LP; %#ok<AGROW>
            if bb, Pst{end+1} = band_legendre(f, f1, fc, Kt); else, Pst{end+1} = cfg.Psit; end %#ok<AGROW>
        end
        stages{1} = mc1;
        stages{end+1} = mc; LPs{end+1} = []; Pst{end+1} = cfg.Psit;      % full band
    elseif Kt > 1
        stages = {mc1, mc}; LPs = {[], []}; Pst = {cfg.Psit, cfg.Psit};
    end
end
hist = struct('lin', {}, 'Psit', {}, 'LP', {});
for k = 1:numel(stages)
    cfk = cfg; cfk.Psit = Pst{k};
    if k > 1 && ~isequal(Pst{k}, Pst{k-1})          % carry tau_n(w) over to this stage's basis
        f = cfg.w(:)/(2*pi); inb = true(size(f));
        if ~isempty(LPs{k}), inb = f <= max(f(Pst{k}(:, 2) < 1 - 1e-12)) + eps; end
        for r = 1:nrec
            tv = Pst{k-1} * linc{r}.tau0.';                % Nf x nm (s)
            linc{r}.tau0 = (Pst{k}(inb, :) \ tv(inb, :)).';
        end
    end
    for r = 1:nrec
        for itr = 1:ngn
            if isempty(LPs{k})
                linc{r} = nuis_relin(rec(r), Cf{rec(r).wave}, linc{r}, cfk, sig, rho, stages{k}, nm, true);
            else
                linc{r} = nuis_relin_band(rec(r), Cf{rec(r).wave}, linc{r}, cfk, sig, rho, stages{k}, nm, LPs{k});
            end
        end
    end
    if nargout > 1, hist(k).lin = linc; hist(k).Psit = Pst{k}; hist(k).LP = LPs{k}; end
end
end

function P = band_legendre(f, f1, fc, Kt)
% Legendre polynomials in normalized log f on [f1, fc], constant above fc
x = (log(min(f(:), fc)) - log(f1)) / (log(fc) - log(f1));
xi = 2*x - 1;
A = [ones(size(xi)), xi, 1.5*xi.^2 - 0.5, 2.5*xi.^3 - 1.5*xi];
P = A(:, 1:Kt);
end
