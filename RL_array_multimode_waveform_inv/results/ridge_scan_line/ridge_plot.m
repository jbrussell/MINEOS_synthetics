% RIDGE_PLOT  figure of ridge_scan.mat: log L along the lines, relative to the truth
addpath('functions');
R = load('ridge_scan.mat'); t = R.t; L0 = R.L(abs(t) < 1e-9, 1);   % truth
set(0, 'defaultfigurevisible', 'off');
fig = figure('position', [50 50 1500 900]);
grp = {1:4, 5:8, 9:14, 15:18};
tit = {'truth (t=0) -> 2b ensembles (t=1)', 'truth (t=0) -> 2f ensembles (t=1)', ...
       '2b ensemble (t=0) -> 2b ensemble (t=1)', '2f mean (t=0) -> 2b ensembles (t=1)'};
col = lines(6);
for g = 1:4
    subplot(2, 3, g); hold on; box on; grid on;
    lg = {};
    for j = 1:numel(grp{g})
        k = grp{g}(j);
        plot(t, R.L(:, k) - L0, '-', 'color', col(j, :), 'linewidth', 1.5);
        lg{end+1} = sprintf('%s - %s', R.names{R.pairs(k, 1)}, R.names{R.pairs(k, 2)}); %#ok<SAGROW>
    end
    plot([0 0], ylim, 'k:'); plot([1 1], ylim, 'k:');
    xlabel('t'); ylabel('log L - log L(truth)'); title(tit{g}); legend(lg, 'location', 'south', 'interpreter', 'none');
end
subplot(2, 3, 5); hold on; box on; grid on;
names = R.names; z = R.zplot;
for i = [1 2:5 10]
    if i == 1, st = 'k-'; lw = 2.5; elseif i == 10, st = 'r-'; lw = 2; else, st = '-'; lw = 1; end
    p = R.prof{i}; h = plot(p(:, 2), z, st, 'linewidth', lw);
    if i > 1 && i < 10, set(h, 'color', [0.4 0.4 0.9]); end
end
set(gca, 'ydir', 'reverse'); ylim([0 350]); xlabel('\xi'); ylabel('depth (km)');
title('mean \xi: truth (black), 2b ensembles (blue), 2f mean (red)');
subplot(2, 3, 6); axis off;
tx = {'Line test, 2b setup (step-1 isotropic reference,', 'tau fixed at the truth, B free and marginalized, dth ~0).', ...
      'C(t) = (1-t) C_A + t C_B between posterior-mean dispersions', '(linear forward: = line between the mean profiles).', ...
      sprintf('sigma, rho: medians of all posterior samples (sigma %.2f-%.2f).', min(R.sig), max(R.sig)), '', ...
      'No barrier on any line: L stays above the lower end point', 'everywhere in 0<t<1 (one connected basin, no separate modes).', ...
      'The truth is not at the top: models with lid xi ~1.03-1.04', 'and xi ~1.04 at 160-250 km fit up to ~240 better.'};
for q = 1:numel(tx), text(0, 1 - 0.08*q, tx{q}, 'fontsize', 10); end
print(fig, '-dpng', '-r110', 'ridge_scan.png');
