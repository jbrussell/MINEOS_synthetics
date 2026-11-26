function fig = plot_posterior_disp(P, outfile, nfine, mode)
% fig = plot_posterior_disp(P, outfile, nfine, mode)
%
% Posterior probability fields of the Love (top row) and Rayleigh (bottom
% row) phase velocities of every mode branch (columns), in the style of the
% Vsv/Vsh/xi fields of plot_posterior_pdf.m (P from
% write_posterior_summary.m, which keeps the dispersion of every pooled
% T = 1 sample in P.samples.CL / .CR).
%
% mode = 'anomaly' (default): plotted as the anomaly relative to the
% reference (step-1) dispersion, 100 (c / c_ref - 1) %, so the posterior
% spread is visible for the fundamental mode and the overtones alike (the
% posterior median is the baseline when there is no reference);
% 'absolute': the phase velocity c (km/s), with the reference curve. Each sampled dispersion curve is
% interpolated linearly in log period onto nfine periods (default 80)
% between the shortest and longest period at which the data constrain the
% branch (P.disp.<L|R>.need); at each period the histogram of the samples,
% smoothed by a Gaussian of 1.5 bins, is normalized to its maximum. Lines:
% 68 % and 95 % intervals, mean, median, true dispersion (synthetic tests)
% and the reference; dots: the periods of the dispersion grid.
% Called by sub_plot_results.m (fig11_posterior_dispersion.png, anomaly;
% fig11b_posterior_dispersion_abs.png, absolute).
if nargin < 3 || isempty(nfine), nfine = 80; end
if nargin < 4 || isempty(mode), mode = 'anomaly'; end
abs_ = strcmpi(mode, 'absolute');
wl = 'LR'; wn = {'Love', 'Rayleigh'};
T = P.disp.T(:);
nmo = size(P.samples.CL, 2);
cm = pdf_cmap(256); cb_ = [0 0.3 0.95];
fig = numbered_figure(1110 + abs_, 'position', [50 50 max(900, 380*nmo) 900], 'color', 'w');
for wv = 1:2
    D = P.disp.(wl(wv)); C = P.samples.(['C' wl(wv)]);         % nT x nmo x ns
    for n = 1:nmo
        ax = subplot(2, nmo, (wv - 1)*nmo + n); hold on; box on;
        ip = find(D.need(:, n));
        if numel(ip) < 2
            axis off; text(0.5, 0.5, sprintf('%s mode %d: not constrained', wn{wv}, n - 1), ...
                'horizontalalignment', 'center'); continue
        end
        Tn = T(ip); X = squeeze(C(ip, n, :));                  % nper x ns
        if size(X, 1) == 1, X = X'; end
        cref = D.ref(ip, n);
        if any(~isfinite(cref)), cref = median(X, 2); rlab = 'posterior median'; else, rlab = 'reference (step 1)'; end
        if abs_                                                 % absolute (km/s)
            A = X; At = D.true(ip, n); Aref = cref;
        else                                                    % anomaly (%)
            A = 100*(X ./ cref - 1); At = 100*(D.true(ip, n) ./ cref - 1); Aref = 0*cref;
        end
        % fine log-period grid; each sample's curve interpolated
        lf = linspace(log(Tn(1)), log(Tn(end)), nfine)';
        Af = interp1(log(Tn), A, lf, 'linear');
        % bins in the anomaly
        lo = min(prctile_rows(Af, 0.5)); hi = max(prctile_rows(Af, 99.5));
        if any(isfinite(At)), lo = min(lo, min(At)); hi = max(hi, max(At)); end
        lo = min(lo, min(Aref)); hi = max(hi, max(Aref));
        pad = 0.08*max(hi - lo, 0.05*(1 + 0.01*abs_)); lo = lo - pad; hi = hi + pad;
        nb = 120 + 130*abs_; bw = (hi - lo)/nb; bins = lo + bw*((1:nb)' - 0.5);
        H = zeros(nb, nfine);
        for j = 1:nfine
            ib = min(max(floor((Af(j, :) - lo)/bw) + 1, 1), nb);
            H(:, j) = accumarray(ib(:), 1, [nb 1]);
        end
        g = exp(-0.5*((-5:5)'/1.5).^2); g = g/sum(g);
        H = conv2(g, 1, H, 'same');
        H = H ./ max(max(H, [], 1), eps);
        Tf = exp(lf);
        % imagesc needs a uniform x axis: draw in log10(T), label in seconds
        imagesc(log10(Tf), bins, H); colormap(ax, cm); caxis([0 1]); set(gca, 'ydir', 'normal');
        x = log10(Tn);
        h = []; l = {};
        h(end+1) = plot(x, Aref, '--', 'color', [0.6 0 0.7], 'linewidth', 1.3); l{end+1} = rlab;
        if any(isfinite(At))
            h(end+1) = plot(x, At, '-', 'color', [0 0.65 0], 'linewidth', 2.2); l{end+1} = 'true';
        end
        h(end+1) = plot(x, mean(A, 2), 'k-', 'linewidth', 1.8); l{end+1} = 'mean';
        h(end+1) = plot(x, prctile_rows(A, 50), 'k:', 'linewidth', 1.6); l{end+1} = 'median';
        h(end+1) = plot(x, prctile_rows(A, 16), '-', 'color', cb_, 'linewidth', 1.1); l{end+1} = '68 %';
        plot(x, prctile_rows(A, 84), '-', 'color', cb_, 'linewidth', 1.1);
        h(end+1) = plot(x, prctile_rows(A, 2.5), '--', 'color', cb_, 'linewidth', 1.1); l{end+1} = '95 %';
        plot(x, prctile_rows(A, 97.5), '--', 'color', cb_, 'linewidth', 1.1);
        plot(x, mean(A, 2), 'k.', 'markersize', 10);
        tk = [10 15 20 25 30 40 50 60 70 80 100 120 150 200];
        tk = tk(tk >= Tn(1)*0.98 & tk <= Tn(end)*1.02);
        set(gca, 'xtick', log10(tk), 'xticklabel', arrayfun(@num2str, tk, 'uniformoutput', false), ...
            'layer', 'top', 'fontsize', 10);
        xlim(log10(Tn([1 end]))' + [-0.01 0.01]); ylim([lo hi]);
        title(sprintf('%s mode %d', wn{wv}, n - 1));
        if n == 1, if abs_, ylabel('c (km/s)'); else, ylabel('c / c_{ref} - 1 (%)'); end, end
        if wv == 2, xlabel('Period (s)'); end
        if wv == 1 && n == 1, legend(h, l, 'location', 'southwest', 'fontsize', 8); end
        if n == nmo && wv == 1
            cbh = colorbar; ylabel(cbh, 'probability (normalized at each period)');
        end
    end
end
annotation('textbox', [0 0.955 1 0.04], 'string', sprintf('Run %s: phase-velocity posterior (%s), %d ensembles, %d samples (T = 1, pooled)', ...
    P.tag, lower(mode), P.nens, P.ns), 'edgecolor', 'none', 'horizontalalignment', 'center', 'fontsize', 11, 'interpreter', 'none');
if nargin >= 2 && ~isempty(outfile)
    save_figure(fig, outfile, 150);
end
end

function q = prctile_rows(X, p)
X = sort(X, 2); n = size(X, 2);
pos = 1 + (n - 1)*p/100;
lo = floor(pos); hi = min(lo + 1, n); fr = pos - lo;
q = (1 - fr)*X(:, lo) + fr*X(:, hi);
end

function c = pdf_cmap(n)
% white -> yellow -> orange -> red -> dark red (as plot_posterior_pdf.m)
k = [0 0.15 0.4 0.7 1];
r = [1 1.00 1.00 0.85 0.45];
g = [1 0.95 0.65 0.20 0.00];
b = [1 0.60 0.15 0.10 0.05];
x = linspace(0, 1, n)';
c = [interp1(k, r, x) interp1(k, g, x) interp1(k, b, x)];
end
