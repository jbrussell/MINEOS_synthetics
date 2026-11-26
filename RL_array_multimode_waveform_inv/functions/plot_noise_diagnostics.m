function plot_noise_diagnostics(rec, cfg, ev, sta, outdir)
% plot_noise_diagnostics(rec, cfg, ev, sta, outdir)
%
% Figures of the noise measured before the first arrival (pre-P noise;
% read_event_data.m, stored in the records by sub_setup_records.m as
% rec.nstd_, rec.npsd_, rec.nacf_; 'aswms' data only):
%   noise_levels.png  per record and station: pre-P noise std / the record's
%                     peak in its windows (map, events x stations), and per
%                     station the median and range over events
%   noise_spectra.png per station (thin) and pooled (median, 10-90%) power
%                     spectral density normalized by the station's variance,
%                     against period
%   noise_acf.png     autocorrelation per station (thin) and pooled, and the
%                     AR(1) model r^k with the pooled lag-1 r
%   noise_absolute.png absolute noise std (data units) per station over the
%                     events (which stations are noisy regardless of event size)
% Love (top/left) and Rayleigh (bottom/right) separately.
% Called by sub_plot_results.m (D4) when the records have the fields; also
% usable on its own: S = load('setup<tag>.mat'); plot_noise_diagnostics(S.rec,
% S.cfg, S.ev, S.sta, 'figs<tag>').
if ~isfield(rec, 'nstd_') || all(arrayfun(@(x) isempty(x.nstd_), rec)), return; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
wn = {'Love', 'Rayleigh'}; Nst = numel(sta);
fg = linspace(cfg.fc(1), cfg.fc(4), 60);
lag = (0:size(rec(1).nacf_, 2) - 1) * cfg.dt;
col = lines(Nst);

% ---- 1: levels relative to the record peak ----
fig = numbered_figure(1501, 'position', [0 0 2300 1000]);
for wv = 1:2
    rr = find([rec.wave] == wv);
    Lv = NaN(numel(rr), Nst);
    for q = 1:numel(rr), Lv(q, :) = rec(rr(q)).nstd_ / rec(rr(q)).scale; end
    subplot(2, 2, 2*wv - 1);
    imagesc(1:Nst, 1:numel(rr), log10(Lv)); axis tight; c = colorbar;
    ylabel(c, 'log_{10}(noise std / record peak)');
    set(gca, 'xtick', 1:Nst, 'xticklabel', sta, 'ytick', 1:numel(rr), 'yticklabel', {ev([rec(rr).event]).id});
    title(sprintf('%s: pre-P noise relative to each record''s peak', wn{wv}));
    subplot(2, 2, 2*wv); hold on; box on;
    for j = 1:Nst
        v = Lv(isfinite(Lv(:, j)), j);
        if isempty(v), continue; end
        plot([j j], [min(v) max(v)], '-', 'color', [0.6 0.6 0.6], 'linewidth', 2);
        plot(j, median(v), 'o', 'color', col(j, :), 'markerfacecolor', col(j, :));
    end
    set(gca, 'yscale', 'log', 'xtick', 1:Nst, 'xticklabel', sta); xlim([0 Nst + 1]);
    ylabel('noise std / record peak'); title(sprintf('%s: per station, median and range over %d events (1/SNR)', wn{wv}, numel(rr)));
end
print(fig, '-dpng', '-r110', fullfile(outdir, 'noise_levels.png')); if fig_close_mode, close(fig); end

% ---- 2: spectra normalized by each station's variance ----
fig = numbered_figure(1502, 'position', [0 0 1400 600]);
for wv = 1:2
    subplot(1, 2, wv); hold on; box on;
    A = [];
    for r = find([rec.wave] == wv)
        for j = 1:Nst
            p = rec(r).npsd_(j, :) / rec(r).nstd_(j)^2;
            if all(isfinite(p)), plot(1./fg, p, '-', 'color', [col(j, :) * 0.5 + 0.5], 'linewidth', 0.3); A = [A; p]; end %#ok<AGROW>
        end
    end
    if ~isempty(A)
        q = sort(A, 1); n = size(q, 1);
        plot(1./fg, q(max(1, round(0.1*n)), :), 'k--', 1./fg, q(round(0.9*n), :), 'k--', 'linewidth', 1.2);
        plot(1./fg, median(A, 1), 'k-', 'linewidth', 2.5);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'xtick', [17 20 30 50 80 112]); xlim([1/cfg.fc(4) 1/cfg.fc(1)]);
    xlabel('Period (s)'); ylabel('PSD / variance (1/Hz)');
    title(sprintf('%s pre-P noise spectra: stations (thin), pooled median and 10-90%% (black)', wn{wv}));
end
print(fig, '-dpng', '-r110', fullfile(outdir, 'noise_spectra.png')); if fig_close_mode, close(fig); end

% ---- 3: autocorrelation ----
fig = numbered_figure(1503, 'position', [0 0 1400 600]);
for wv = 1:2
    subplot(1, 2, wv); hold on; box on;
    A = [];
    for r = find([rec.wave] == wv)
        for j = 1:Nst
            a = rec(r).nacf_(j, :);
            if all(isfinite(a)), plot(lag, a, '-', 'color', [col(j, :) * 0.5 + 0.5], 'linewidth', 0.3); A = [A; a]; end %#ok<AGROW>
        end
    end
    h = [];
    if ~isempty(A)
        m = median(A, 1);
        h(1) = plot(lag, m, 'k-', 'linewidth', 2.5);
        h(2) = plot(lag, m(2).^(0:numel(lag)-1), 'r--', 'linewidth', 1.8);
        legend(h, {'pooled (median)', sprintf('AR(1) model, r = %.2f', m(2))}, 'location', 'northeast');
    end
    plot(lag([1 end]), [0 0], 'k:');
    xlabel('Lag (s)'); ylabel('autocorrelation'); ylim([-1 1]);
    title(sprintf('%s pre-P noise autocorrelation: stations (thin), pooled (black)', wn{wv}));
end
print(fig, '-dpng', '-r110', fullfile(outdir, 'noise_acf.png')); if fig_close_mode, close(fig); end

% ---- 4: absolute levels per station ----
fig = numbered_figure(1504, 'position', [0 0 1400 600]);
for wv = 1:2
    subplot(1, 2, wv); hold on; box on;
    rr = find([rec.wave] == wv); A = NaN(numel(rr), Nst);
    for q = 1:numel(rr), A(q, :) = rec(rr(q)).nstd_; end
    for j = 1:Nst
        v = A(isfinite(A(:, j)), j);
        if isempty(v), continue; end
        plot(j + 0.15*randn(size(v)), v, '.', 'color', col(j, :), 'markersize', 10);
        plot(j + [-0.3 0.3], median(v)*[1 1], 'k-', 'linewidth', 2);
    end
    set(gca, 'yscale', 'log', 'xtick', 1:Nst, 'xticklabel', sta); xlim([0 Nst + 1]);
    ylabel('pre-P noise std (data units)'); title(sprintf('%s: absolute noise per station (dots: events, bar: median)', wn{wv}));
end
print(fig, '-dpng', '-r110', fullfile(outdir, 'noise_absolute.png')); if fig_close_mode, close(fig); end
fprintf('noise diagnostics (pre-P noise) in %s/noise_*.png\n', outdir);
end
