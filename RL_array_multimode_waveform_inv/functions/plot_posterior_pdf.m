function fig = plot_posterior_pdf(P, outfile, zmax_plot)
% fig = plot_posterior_pdf(P, outfile, zmax_plot)
%
% Posterior probability fields of Vsv, Vsh and xi with depth, in the style
% of Eilon et al. (2018): at each depth the marginal posterior histogram of
% the pooled samples (P from write_posterior_summary.m), normalized to its
% maximum at that depth (so narrow and broad parts of the profile are both
% visible), with the middle 68 % and 95 % intervals, the posterior mean and
% median, the reference (step-1) model and, for synthetic tests, the true
% model. The inverted depth range (P.zinv) is marked. Called by
% sub_plot_results.m (fig10_posterior_pdf.png).
%
% zmax_plot : bottom of the plot (km below the seafloor; default 400)
if nargin < 3 || isempty(zmax_plot), zmax_plot = 400; end
names = {'vsv', 'vsh', 'xi'};
labs = {'V_{SV} (km/s)', 'V_{SH} (km/s)', '\xi = (V_{SH}/V_{SV})^2'};
cm = pdf_cmap(256);
z = P.z; iz = z <= zmax_plot;
fig = numbered_figure(1100, 'position', [50 50 1500 700], 'color', 'w');
for q = 1:3
    S = P.(names{q});
    ax = subplot(1, 3, q); hold on; box on;
    % field: the histogram smoothed with a Gaussian kernel (1.5 bins along
    % the parameter axis, 1 grid step in depth) so a few hundred samples give
    % a continuous density, then normalized to its maximum at each depth
    H = smooth2(S.pdf, 1.5, 1);
    H = H ./ max(max(H, [], 2), eps);
    % x range: where the field is non-negligible within the inverted depth
    % range (the water and crust above it are fixed and would set the scale)
    zr = z >= P.zinv(1) & z <= min(zmax_plot, P.zinv(2) + 20);
    keep = any(H(zr, :) > 0.02, 1);
    xl = S.bins([find(keep, 1) find(keep, 1, 'last')])';
    xr = [S.ref(zr); S.true(zr)]; xr = xr(isfinite(xr));
    if ~isempty(xr), xl = [min(xl(1), min(xr)) max(xl(2), max(xr))]; end
    pad = 0.05*diff(xl); xl = xl + [-pad pad];
    imagesc(S.bins, z, H);
    colormap(ax, cm); caxis([0 1]);
    cb_ = [0 0.3 0.95];                                       % interval/mean colour
    h = []; l = {};
    if any(isfinite(S.true))
        h(end+1) = plot(S.true, z, '-', 'color', [0 0.65 0], 'linewidth', 2.2); l{end+1} = 'true';
    end
    h(end+1) = plot(S.ref, z, '--', 'color', [0.6 0 0.7], 'linewidth', 1.3); l{end+1} = 'reference (step 1)';
    h(end+1) = plot(S.mean, z, '-', 'color', 'k', 'linewidth', 1.8); l{end+1} = 'mean';
    h(end+1) = plot(S.median, z, ':', 'color', 'k', 'linewidth', 1.6); l{end+1} = 'median';
    h(end+1) = plot(S.ci68(:, 1), z, '-', 'color', cb_, 'linewidth', 1.1); l{end+1} = '68 %';
    plot(S.ci68(:, 2), z, '-', 'color', cb_, 'linewidth', 1.1);
    h(end+1) = plot(S.ci95(:, 1), z, '--', 'color', cb_, 'linewidth', 1.1); l{end+1} = '95 %';
    plot(S.ci95(:, 2), z, '--', 'color', cb_, 'linewidth', 1.1);
    for zz = P.zinv(:)', plot(xl, zz*[1 1], 'k-.', 'linewidth', 0.5); end
    set(gca, 'ydir', 'reverse', 'layer', 'top', 'fontsize', 11);
    xlim(xl); ylim([0 zmax_plot]); xlabel(labs{q});
    if q == 1, ylabel('Depth below seafloor (km)'); end
    title(sprintf('Posterior %s', labs{q}));
    if q == 1, legend(h, l, 'location', 'southwest'); end
    if q == 3
        cb = colorbar; ylabel(cb, 'probability (normalized at each depth)');
    end
end
annotation('textbox', [0 0.955 1 0.04], 'string', sprintf('Run %s: %d ensembles, %d samples (T = 1, pooled); dash-dot: inverted depth range', ...
    P.tag, P.nens, P.ns), 'edgecolor', 'none', 'horizontalalignment', 'center', 'fontsize', 11, 'interpreter', 'none');
if nargin >= 2 && ~isempty(outfile)
    save_figure(fig, outfile, 150);
end
end

function H = smooth2(H, sx, sz)
% Gaussian smoothing of the rows (sx bins) and columns (sz depth steps)
gx = exp(-0.5*((-ceil(3*sx):ceil(3*sx))/sx).^2); gx = gx/sum(gx);
gz = exp(-0.5*((-ceil(3*sz):ceil(3*sz))'/sz).^2); gz = gz/sum(gz);
H = conv2(gz, gx, H, 'same');
end

function c = pdf_cmap(n)
% white -> yellow -> orange -> red -> dark red
k = [0 0.15 0.4 0.7 1];
r = [1 1.00 1.00 0.85 0.45];
g = [1 0.95 0.65 0.20 0.00];
b = [1 0.60 0.15 0.10 0.05];
x = linspace(0, 1, n)';
c = [interp1(k, r, x) interp1(k, g, x) interp1(k, b, x)];
end
