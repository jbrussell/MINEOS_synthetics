function plot_record_sections(P, outdir, fc, win)
% plot_record_sections(P, outdir, fc, win)
%
% Record sections of the step-2 prep structure P (prep_mineos_idagrn.m):
% the data (P.u_true_*: the event files of A1-A3 in D2, or the full PA5
% mode sum, branches 0..P.nbranch_data-1), the step-1 starting model (if
% present: P.u_start_*, opt.start_synth; same branches as the data) and the
% final step-1 model (the reference: the sum of its single-branch
% synthetics P.u_ref_*, branches 0..P.nm-1), transverse (Love, T) and vertical (Rayleigh, Z)
% displacement, band-passed to the inversion band fc = [f1 f2 f3 f4] Hz
% (cosine taper; default [1/112 1/100 1/20 1/17]).
%
% win (optional): the automatic waveform windows of A3 (ASWMS-Q,
%   A3_run_autowinpick.m), win{wv}{e} for wave type wv (1 Love, 2 Rayleigh)
%   and event e: a struct with winpara [v1 t1 v2 t2], dist (km, per
%   station of P.sta) and toff (event.otime - centroid time, s), as
%   read_event_data.m returns them; the window of a station is
%   dist/v1 + t1 + toff .. dist/v2 + t2 + toff on the time axis of P.t.
%   Empty (or missing for an event): no window is drawn.
%
% Figures in outdir:
%   record_section_all_<T|Z>.png    all events and stations vs epicentral
%                                   distance, each trace normalized; the
%                                   window bounds of every event
%   record_section_events_<T|Z>.png one panel per event (stations vs
%                                   distance), normalized per event; data
%                                   black, step-1 start green, final step-1
%                                   model red, window blue
if nargin < 3 || isempty(fc), fc = [1/112 1/100 1/20 1/17]; end
if nargin < 4, win = {}; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
t = P.t(:); Nt = numel(t); dt = t(2) - t(1);
f = (0:Nt-1)'/(Nt*dt);
H = bp_taper(min(f, 1/dt - f), fc);                % two-sided band-pass
tap = tukeywin_(Nt, 0.05);
bp = @(u) real(ifft(fft(u .* tap) .* H));
nev = numel(P.ev); nst = numel(P.sta.lat);
km2deg = 180/(6371*pi);
comp = {'T', 'Z'}; ty = {'T', 'S'}; wname = {'Love (transverse)', 'Rayleigh (vertical)'};
ecol = lines(nev);
for wv = 1:2
    D = zeros(nev, nst); Ut = cell(nev, 1); Us = cell(nev, 1); U0 = cell(nev, 1);
    W1 = NaN(nev, nst); W2 = NaN(nev, nst);       % window start/end per station (s)
    for e = 1:nev
        D(e, :) = gc_dist_az(P.sta.lat, P.sta.lon, P.ev(e).lat*ones(nst, 1), P.ev(e).lon*ones(nst, 1))';
        Ut{e} = bp(P.(['u_true_' ty{wv}]){e});
        fr = ['u_ref_' ty{wv}];
        if isfield(P, fr)
            us = 0; for n = 1:size(P.(fr), 2), us = us + P.(fr){e, n}; end
            Us{e} = bp(us);
        end
        if isfield(P, ['u_start_' ty{wv}]), U0{e} = bp(P.(['u_start_' ty{wv}]){e}); end
        w = [];
        if numel(win) >= wv && numel(win{wv}) >= e, w = win{wv}{e}; end
        if ~isempty(w) && numel(w.winpara) == 4
            wp = w.winpara;
            W1(e, :) = w.dist(:)'/wp(1) + wp(2) + w.toff;
            W2(e, :) = w.dist(:)'/wp(3) + wp(4) + w.toff;
        end
    end
    hasw = any(isfinite(W1(:)));
    tit_ = sprintf('%s, %g-%g s: data (branches 0-%d), step-1 start (green) and final step-1 model (red, branches 0-%d)', ...
        wname{wv}, round(1/fc(4)), round(1/fc(1)), P.nbranch_data - 1, P.nm - 1);
    if hasw, wl_ = 'blue: automatic window (A3)'; else, wl_ = 'no automatic windows (run A3 before D2)'; end
    % ---- all events in one record section ----
    fig = numbered_figure(2000 + 2*wv - 1, 'position', [50 50 1100 1300]); hold on; box on;
    sc = 0.45;                                      % half-amplitude of a trace (deg)
    for e = 1:nev
        for j = 1:nst
            u = Ut{e}(:, j); u = u / max(abs(u));
            plot(t, D(e, j) + sc*u, '-', 'color', ecol(e, :), 'linewidth', 0.5);
        end
        if any(isfinite(W1(e, :)))                  % window bounds of the event (stations by distance)
            [ds, is] = sort(D(e, :));
            plot(W1(e, is), ds, 'k-', W2(e, is), ds, 'k-', 'linewidth', 1);
        end
    end
    dd = [floor(min(D(:))) - 1, ceil(max(D(:))) + 1];
    set(gca, 'ydir', 'normal'); xlim([0 t(end)]); ylim(dd);
    xlabel('Time since the centroid (s)'); ylabel('Epicentral distance (deg)');
    if hasw, wl2_ = 'black: automatic window of each event (A3)'; else, wl2_ = wl_; end
    title(sprintf('%s, data, %g-%g s, traces normalized; %s', wname{wv}, round(1/fc(4)), round(1/fc(1)), wl2_));
    save_figure(fig, fullfile(outdir, ['record_section_all_' comp{wv} '.png']), 130);
    if fig_close_mode, close(fig); end
    % ---- one panel per event ----
    fig = numbered_figure(2000 + 2*wv, 'position', [50 50 1900 1250]);
    nr = 3; nc = ceil(nev/nr);
    for e = 1:nev
        ir = ceil(e/nc); icol = e - (ir - 1)*nc;     % tight panel layout
        axes('position', [0.04 + (icol - 1)*0.192, 0.07 + (nr - ir)*0.305, 0.17, 0.255]); hold on; box on;
        de = D(e, :); rng_ = max(de) - min(de);
        sp = max(rng_, 0.5) / nst * 1.6;            % half-amplitude of a trace (deg)
        sn = max(abs(Ut{e}(:)));                    % one normalization per event
        for j = 1:nst
            if ~isempty(U0{e}), plot(t, de(j) + sp*U0{e}(:, j)/sn, '-', 'color', [0 0.6 0], 'linewidth', 0.5); end
            if ~isempty(Us{e}), plot(t, de(j) + sp*Us{e}(:, j)/sn, '-', 'color', [0.85 0.2 0.2], 'linewidth', 0.5); end
            plot(t, de(j) + sp*Ut{e}(:, j)/sn, 'k-', 'linewidth', 0.6);
        end
        dl = [min(de) - 2*sp, max(de) + 2*sp];
        if any(isfinite(W1(e, :)))
            [ds, is] = sort(de);
            plot(W1(e, is), ds, 'b-', W2(e, is), ds, 'b-', 'linewidth', 1.2);
            t1 = max(0, min(W1(e, :)) - 300); t2 = min(t(end), max(W2(e, :)) + 300);
        else
            t1 = max(0, min(de)/km2deg/5 - 300); t2 = min(t(end), max(de)/km2deg/3.5 + 400);
        end
        xlim([t1 t2]); ylim(dl);
        title(sprintf('%s  %.0f-%.0f deg', P.ev(e).id, min(de), max(de)), 'interpreter', 'none');
        if mod(e - 1, nc) == 0, ylabel('\Delta (deg)'); end
        if e > (nr - 1)*nc, xlabel('Time since the centroid (s)'); end
    end
    st_ = sprintf('%s; data black; %s', tit_, wl_);
    if exist('OCTAVE_VERSION', 'builtin')
        axes('position', [0 0.965 1 0.03], 'visible', 'off');
        text(0.5, 0.5, st_, 'horizontalalignment', 'center');
    else
        sgtitle(st_);
    end
    save_figure(fig, fullfile(outdir, ['record_section_events_' comp{wv} '.png']), 130);
    if fig_close_mode, close(fig); end
end
end
