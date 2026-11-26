function P = write_posterior_summary(svc, zplot, cfg, mc, tru, refp, Tdisp, tag, outfile)
% P = write_posterior_summary(svc, zplot, cfg, mc, tru, refp, Tdisp, tag, outfile)
%
% Final posterior of a step-2 run: the T = 1 samples of ALL chains/PT
% ensembles pooled (equal weight per saved sample) and their statistics,
% saved to outfile (posterior_summary<tag>.mat; called by
% sub_plot_results.m). Read it with  P = load(outfile); P = P.P;
%
% svc : {nens x 1} saved T = 1 samples per ensemble (chain<tag>_XX.mat, sv)
%
% P.z                 depth below the seafloor (km, zplot), nz x 1
% P.zinv              [zmin zmax] (km) of the inverted part of the model
%                     (above zmax the model is sampled; below, the reference)
% P.<p>, p = vsv, vsh (km/s), xi (= (Vsh/Vsv)^2), and eta when it is
%   sampled, per depth (nz x 1, or nz x 2 for the intervals):
%   .mean .median .std
%   .ci95  [2.5 97.5] percentiles  (middle 95 %)
%   .ci68  [16 84] percentiles     (middle 68 %)
%   .rhat  Gelman-Rubin R-hat across the ensembles (NaN where the
%          parameter does not vary, e.g. below zmax)
%   .pdf, .bins  per-depth histogram (nz x nbins, rows sum to 1) and the
%          bin centres (the probability fields of plot_posterior_pdf.m)
%   .true, .ref  true model (synthetic tests; NaN otherwise) and the
%          reference (step-1) model
% P.ncp.vsv, P.ncp.xi  number of interior control points (Vsv and xi
%   B-spline knots: sv.kv, sv.kx): .samples (ns x 1), .mean .median .mode
%   .ci95 .ci68 .rhat, .values/.prob (posterior probability of each count)
% P.disp.T (s), P.disp.L / .R  Love / Rayleigh phase velocities of the
%   sampled models (km/s), nT x nmodes: .mean .median .ci95 (nT x nmodes x 2)
%   .ci68 .rhat
% P.noise.sig, .rho  noise standard deviations and correlations (ns x 4;
%   L ovt, L fund, R ovt, R fund); P.smooth_sig (ns x 2) and P.xi_sig
%   (ns x 1) when sampled (hierarchical priors)
% P.samples  pooled samples: .vsv .vsh .xi (.eta) (nz x ns), .ensemble
%   (1 x ns, ensemble of each sample), .logL (1 x ns, NaN if unknown)
% P.ns, P.nens, P.ns_per_ensemble, P.tag, P.created
nc = numel(svc);
cat2 = @(f) cell2mat(cellfun(@(x) x.(f), svc', 'uniformoutput', false));   % pooled along dim 2
catr = @(f) cell2mat(cellfun(@(x) x.(f), svc, 'uniformoutput', false));    % pooled along dim 1
nsi = cellfun(@(x) size(x.vsv, 2), svc);
ns = sum(nsi);

P = struct();
P.tag = tag; P.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
P.nens = nc; P.ns = ns; P.ns_per_ensemble = nsi(:)';
P.z = zplot(:);
P.zinv = [cfg.zmin cfg.zmax];
P.description = ['Pooled T = 1 posterior samples of all ensembles (equal weight) and their ' ...
                 'statistics; see write_posterior_summary.m'];
names = {'vsv', 'vsh', 'xi'};
if isfield(svc{1}, 'eta') && isfield(mc, 'eta_free') && mc.eta_free, names{end+1} = 'eta'; end
bw = struct('vsv', 0.01, 'vsh', 0.01, 'xi', 0.0025, 'eta', 0.0025);   % histogram bin widths
ens = cell2mat(arrayfun(@(i) i*ones(1, nsi(i)), 1:nc, 'uniformoutput', false));
P.samples = struct('ensemble', ens);
for q = 1:numel(names)
    p = names{q};
    X = cat2(p);
    S = struct();
    S.mean = mean(X, 2); S.median = pctl(X, 50); S.std = std(X, 0, 2);
    S.ci95 = [pctl(X, 2.5) pctl(X, 97.5)];
    S.ci68 = [pctl(X, 16) pctl(X, 84)];
    S.rhat = NaN(numel(zplot), 1);
    if nc > 1 && all(nsi == nsi(1))
        for iz = 1:numel(zplot)
            S.rhat(iz) = gr_rhat(cell2mat(cellfun(@(x) x.(p)(iz, :)', svc', 'uniformoutput', false)));
        end
    end
    lo = floor(min(S.ci95(:, 1) - 3*S.std)/bw.(p))*bw.(p);
    hi = ceil(max(S.ci95(:, 2) + 3*S.std)/bw.(p))*bw.(p);
    lo = min(lo, min(X(:))); hi = max(hi, max(X(:)));
    S.bins = (floor(lo/bw.(p))*bw.(p):bw.(p):ceil(hi/bw.(p))*bw.(p))';
    S.pdf = zeros(numel(zplot), numel(S.bins));
    for iz = 1:numel(zplot)
        ib = min(max(round((X(iz, :) - S.bins(1))/bw.(p)) + 1, 1), numel(S.bins));
        S.pdf(iz, :) = accumarray(ib(:), 1, [numel(S.bins) 1])' / ns;
    end
    S.true = NaN(numel(zplot), 1); S.ref = NaN(numel(zplot), 1);
    if isfield(tru, p), S.true = tru.(p)(:); end
    if isfield(refp, p), S.ref = refp.(p)(:); end
    P.(p) = S;
    P.samples.(p) = X;
end

% number of interior control points
P.ncp = struct();
fk = {'vsv', 'kv'; 'xi', 'kx'};
for q = 1:2
    k = catr(fk{q, 2}); k = k(:);
    S = struct('samples', k, 'mean', mean(k), 'median', median(k), 'mode', mode(k), ...
               'ci95', [pctl(k', 2.5) pctl(k', 97.5)], 'ci68', [pctl(k', 16) pctl(k', 84)], 'rhat', NaN);
    if nc > 1 && all(nsi == nsi(1))
        S.rhat = gr_rhat(reshape(k, nsi(1), nc));
    end
    S.values = (min(k):max(k))';
    S.prob = arrayfun(@(v) mean(k == v), S.values);
    P.ncp.(fk{q, 1}) = S;
end

% dispersion of the sampled models
P.disp = struct('T', Tdisp(:));
wl = 'LR';
for wv = 1:2
    C = cat(3, svc{1}.(['C' wl(wv)]));
    for ic = 2:nc, C = cat(3, C, svc{ic}.(['C' wl(wv)])); end
    [nT, nmo, ~] = size(C);
    Cr = reshape(C, nT*nmo, ns);
    S = struct('mean', reshape(mean(Cr, 2), nT, nmo), 'median', reshape(pctl(Cr, 50), nT, nmo), ...
               'ci95', reshape([pctl(Cr, 2.5) pctl(Cr, 97.5)], nT, nmo, 2), ...
               'ci68', reshape([pctl(Cr, 16) pctl(Cr, 84)], nT, nmo, 2), 'rhat', NaN(nT, nmo));
    if nc > 1 && all(nsi == nsi(1))
        for j = 1:nT*nmo
            S.rhat(j) = gr_rhat(reshape(Cr(j, :), nsi(1), nc));
        end
    end
    if isfield(cfg, 'need'), S.need = cfg.need; end   % (mode, period) pairs constrained by the data
    P.disp.(wl(wv)) = S;
end

% noise and hierarchical prior parameters
P.noise = struct('sig', catr('sig'), 'rho', catr('rho'), 'order', {{'L ovt', 'L fund', 'R ovt', 'R fund'}});
if isfield(svc{1}, 'ss'), ss = catr('ss'); if any(isfinite(ss(:))), P.smooth_sig = ss; end, end
if isfield(svc{1}, 'xs')
    xs = cell2mat(cellfun(@(x) x.xs(:), svc, 'uniformoutput', false));
    if any(isfinite(xs)), P.xi_sig = xs; end
end

% log-likelihood of each saved sample (sv.L is the trace of every iteration)
P.samples.logL = NaN(1, ns);
if isfield(mc, 'burn') && isfield(mc, 'thin')
    j0 = 0;
    for ic = 1:nc
        L = svc{ic}.L(:)'; it = mc.burn + mc.thin*(1:nsi(ic));
        if numel(L) >= it(end), P.samples.logL(j0 + (1:nsi(ic))) = L(it); end
        j0 = j0 + nsi(ic);
    end
end

save(outfile, 'P', '-v7');
end

function q = pctl(X, p)
% percentile of each row of X across columns (linear interpolation)
X = sort(X, 2); n = size(X, 2);
pos = 1 + (n - 1)*p/100;
lo = floor(pos); hi = min(lo + 1, n); fr = pos - lo;
q = (1 - fr)*X(:, lo) + fr*X(:, hi);
end
