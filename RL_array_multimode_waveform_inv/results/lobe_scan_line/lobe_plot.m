R = load('lobe_scan.mat'); e = R.eps*100;
set(0, 'defaultfigurevisible', 'off');
cases = {'oracle_s1', 'oracle_s3', 'b2_s1', 'b2_s3'};
lab = {'oracle (B, tau, dth fixed), sigma_{true}', 'oracle, 3 sigma_{true}', '2b (B free, tau fixed, dth 0), sigma_{true}', '2b, 3 sigma_{true}'};
col = [0 0 0; 0.5 0.5 0.5; 0.85 0.2 0.2; 0.95 0.6 0.5];
wn = {'Love', 'Rayleigh'};
fig = figure('position', [50 50 1800 900]);
for wv = 1:2, for n = 1:4
    subplot(2, 4, (wv-1)*4 + n); hold on; box on;
    for c = 1:4
        L = R.(cases{c})(:, n, wv); i0 = find(abs(e) < 1e-9);
        plot(e, L - L(i0), '-', 'color', col(c,:), 'linewidth', 1.2);
    end
    ylim([-3000 200]); xlim([-3 3]); grid on;
    title(sprintf('%s branch %d', wn{wv}, n-1)); xlabel('shift of c_n (%)'); if n == 1, ylabel('logL - logL(truth)'); end
end, end
legend(lab, 'location', 'south', 'fontsize', 7);
print(fig, '-dpng', '-r110', 'figs_line2b/lobe_scan.png');
% local maxima within +-3 % (ΔlogL > -5000), each case
for c = 1:4
    fprintf('== %s\n', cases{c});
    for wv = 1:2, for n = 1:4
        L = R.(cases{c})(:, n, wv); i0 = find(abs(e) < 1e-9); d = L - L(i0);
        pk = find(d(2:end-1) > d(1:end-2) & d(2:end-1) >= d(3:end)) + 1;
        pk = pk(d(pk) > -5000);
        [~, ib] = max(d);
        fprintf('%s%d: best %+.3f %% (dL %+.0f); maxima: %s\n', wn{wv}(1), n-1, e(ib), d(ib), ...
            strjoin(arrayfun(@(i) sprintf('%+.2f%%(%+.0f)', e(i), d(i)), pk(:)', 'uniformoutput', false), ' '));
    end, end
end
