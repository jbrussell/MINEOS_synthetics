% A4_PLOT_EVENT_DATA  Data preparation, step 4 (synthetic and real data):
% figures of the event data after A1-A3 (ASWMS-Q event files), with noise
% diagnostics per event, before the inversion (D1, D2, D3).
%
% For each component of opt.aswms_components and each event file
% <opt.aswms_dir>eventmat/<event id>_<component>.mat:
%   figs_events/<component>/<event id>.png
%     (a) record section: traces band-passed to the inversion band (opt.fc),
%         each normalized by its peak, against epicentral distance; the A3
%         window dist/v1 + t1 .. dist/v2 + t2 (blue) and the end of the
%         pre-P noise window dist/13.5 km/s - 60 s (red); stations that A3
%         flagged (stadata.isgood < 0) in grey
%     (b) SNR per station: peak of the band-passed trace in the window /
%         std of its pre-P noise (noise_stats.m, as used by the inversion)
%     (c) pre-P noise spectra per station, normalized by the variance, and
%         their median
%     (d) pre-P noise autocorrelation per station, their median and the
%         AR(1) model r^k with the median lag-1 r
%   figs_events/<component>/summary.png
%     SNR map (events x stations, log10) and per-station SNR median and
%     range over the events
% Times are relative to event.otime (SAC origin). Output in opt.aswms_dir.
% Previous: A3_run_autowinpick.m (optional A3b). Next: D1 (step 1), D2.

parameters_RL;
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
wd = opt.aswms_dir; if wd(end) ~= '/' && wd(end) ~= filesep, wd = [wd filesep]; end
fc = optget(opt, 'fc', [1/112 1/100 1/20 1/17]);
maxlag = optget(opt, 'noise_maxlag', 12);
fg = linspace(fc(1), fc(4), 60);
for icomp = 1:numel(opt.aswms_components)
    comp = opt.aswms_components{icomp};
    od = [wd 'figs_events' filesep comp filesep];
    if ~exist(od, 'dir'), mkdir(od); end
    files = dir([wd 'eventmat' filesep '*_' comp '.mat']);
    allsta = {}; SNR = {}; evids = {};
    for ie = 1:numel(files)
        E = load([wd 'eventmat' filesep files(ie).name]); event = E.event;
        ns = numel(event.stadata);
        dist = [event.stadata.dist]; [~, so] = sort(dist);
        wp = []; if isfield(event, 'winpara') && numel(event.winpara) == 4, wp = event.winpara; end
        snr = NaN(1, ns); sd = NaN(1, ns); AC = NaN(ns, maxlag + 1); PS = NaN(ns, numel(fg));
        U = cell(1, ns); T = U; good = true(1, ns);
        for k = 1:ns
            s = event.stadata(k);
            if isfield(s, 'isgood') && ~isempty(s.isgood) && s.isgood < 0, good(k) = false; end
            x = s.data(:) - mean(s.data); n = numel(x); dt = s.delta;
            t = (s.otime - event.otime) + (0:n-1)'*dt;
            f = (0:n-1)'/(n*dt); ib = find(f > fc(1) & f < fc(4));
            X = fft(x .* tukeywin_(n, 0.05)); Y = zeros(n, 1); Y(ib) = X(ib) .* bp_taper(f(ib), fc);
            u = 2*real(ifft(Y)); U{k} = u; T{k} = t;
            ip = find(t < s.dist/13.5 - 60);
            if numel(ip) >= optget(opt, 'noise_min_n', 40)
                [sd(k), AC(k, :), PS(k, :)] = noise_stats(s.data(ip), dt, fc, maxlag, fg);
            end
            if ~isempty(wp)
                iw = t >= s.dist/wp(1) + wp(2) & t <= s.dist/wp(3) + wp(4);
            else
                iw = t >= s.dist/5 & t <= s.dist/2.5;     % (no A3 window: 5-2.5 km/s)
            end
            if any(iw), snr(k) = max(abs(u(iw))) / sd(k); end
        end
        % ---- per-event figure ----
        fig = figure('position', [0 0 1600 1100]);
        subplot(2, 2, 1); hold on; box on;
        dd = max(dist) - min(dist); sc = 0.04*max(dd, 100);
        for k = so
            u = U{k}; c = [0 0 0]; if ~good(k), c = [0.6 0.6 0.6]; end
            plot(T{k}, dist(k) + sc*u/max(abs(u)), '-', 'color', c, 'linewidth', 0.5);
        end
        yl = [min(dist) - 2*sc, max(dist) + 2*sc]; dl = linspace(yl(1), yl(2), 20);
        if ~isempty(wp)
            plot(dl/wp(1) + wp(2), dl, 'b-', dl/wp(3) + wp(4), dl, 'b-', 'linewidth', 1.5);
        end
        plot(dl/13.5 - 60, dl, 'r--', 'linewidth', 1.2);
        tmax = max(cellfun(@(t) t(end), T));
        xlim([0, min(tmax, max(dist)/2 + 600)]); ylim(yl);
        xlabel('Time after origin (s)'); ylabel('Epicentral distance (km)');
        title(sprintf('%s %s: band-passed %.0f-%.0f s (blue: A3 window, red: end of pre-P noise)', event.id, comp, 1/fc(4), 1/fc(1)), 'interpreter', 'none');
        subplot(2, 2, 2); hold on; box on;
        stn = {event.stadata.stnm};
        for q = 1:ns
            k = so(q); c = [0.2 0.4 0.8]; if ~good(k), c = [0.6 0.6 0.6]; end
            if isfinite(snr(k)), plot([q q], [1 snr(k)], '-', 'color', c, 'linewidth', 3); plot(q, snr(k), 'o', 'color', c, 'markerfacecolor', c); end
        end
        set(gca, 'yscale', 'log', 'xtick', 1:ns, 'xticklabel', stn(so), 'fontsize', 7); xlim([0 ns + 1]);
        if any(isfinite(snr)), ylim([1, 10^ceil(log10(max(snr)))]); end
        ylabel('SNR (window peak / pre-P noise std)'); title('SNR per station (sorted by distance)');
        subplot(2, 2, 3); hold on; box on;
        okp = all(isfinite(PS), 2);
        for k = find(okp)', plot(1./fg, PS(k, :) / sd(k)^2, '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.5); end
        if any(okp), plot(1./fg, median(PS(okp, :) ./ sd(okp)'.^2, 1), 'k-', 'linewidth', 2); end
        set(gca, 'xscale', 'log', 'yscale', 'log', 'xtick', [17 20 30 50 80 112]); xlim([1/fc(4) 1/fc(1)]);
        xlabel('Period (s)'); ylabel('PSD / variance (1/Hz)'); title('Pre-P noise spectra: stations (grey), median (black)');
        subplot(2, 2, 4); hold on; box on;
        oka = all(isfinite(AC), 2); lag = (0:maxlag)*event.stadata(1).delta; h = [];
        for k = find(oka)', plot(lag, AC(k, :), '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.5); end
        if any(oka)
            m = median(AC(oka, :), 1);
            h(1) = plot(lag, m, 'k-', 'linewidth', 2); h(2) = plot(lag, m(2).^(0:maxlag), 'r--', 'linewidth', 1.5);
            legend(h, {'median', sprintf('AR(1), r = %.2f', m(2))}, 'location', 'northeast');
        end
        plot(lag([1 end]), [0 0], 'k:'); ylim([-1 1]);
        xlabel('Lag (s)'); ylabel('autocorrelation'); title('Pre-P noise autocorrelation');
        print(fig, '-dpng', '-r100', [od event.id '.png']); close(fig);
        % ---- for the summary ----
        evids{end+1} = event.id; row = containers.Map(stn, num2cell(snr)); %#ok<SAGROW>
        allsta = union(allsta, stn); SNR{end+1} = row; %#ok<SAGROW>
        fprintf('%s %s: %d stations, SNR median %.1f (min %.1f), %d without pre-P noise\n', event.id, comp, ns, ...
            median(snr(isfinite(snr))), min(snr(isfinite(snr))), sum(~isfinite(sd)));
    end
    if isempty(evids), continue; end
    % ---- summary over events ----
    M = NaN(numel(evids), numel(allsta));
    for ie = 1:numel(evids)
        for j = 1:numel(allsta)
            if isKey(SNR{ie}, allsta{j}), M(ie, j) = SNR{ie}(allsta{j}); end
        end
    end
    fig = figure('position', [0 0 1500 700]);
    subplot(1, 2, 1);
    imagesc(1:numel(allsta), 1:numel(evids), log10(M)); c = colorbar; ylabel(c, 'log_{10} SNR');
    set(gca, 'xtick', 1:numel(allsta), 'xticklabel', allsta, 'ytick', 1:numel(evids), 'yticklabel', evids, 'fontsize', 6);
    title(sprintf('%s: SNR per event and station', comp));
    subplot(1, 2, 2); hold on; box on;
    for j = 1:numel(allsta)
        v = M(isfinite(M(:, j)), j);
        if isempty(v), continue; end
        plot([j j], [min(v) max(v)], '-', 'color', [0.6 0.6 0.6], 'linewidth', 2);
        plot(j, median(v), 'o', 'color', [0.2 0.4 0.8], 'markerfacecolor', [0.2 0.4 0.8]);
    end
    set(gca, 'yscale', 'log', 'xtick', 1:numel(allsta), 'xticklabel', allsta, 'fontsize', 7); xlim([0 numel(allsta) + 1]);
    ylabel('SNR'); title(sprintf('%s: per-station SNR, median and range over %d events', comp, numel(evids)));
    print(fig, '-dpng', '-r110', [od 'summary.png']); close(fig);
    fprintf('%s: event figures and summary in %s\n', comp, od);
end
