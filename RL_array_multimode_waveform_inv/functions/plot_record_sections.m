function plot_record_sections(P, outdir, fc, U_start)
% plot_record_sections(P, outdir, fc, U_start)
%
% Record sections of the full mode-sum idagrn6 synthetics of the step-2
% prep structure P (prep_mineos_idagrn.m): the PA5 "data" (all branches
% 0..P.nbranch_data-1) and, if present, the step-1 starting model (same
% branches; P.u_start_*), transverse
% (Love, T) and vertical (Rayleigh, Z) displacement, band-passed to the
% inversion band fc = [f1 f2 f3 f4] Hz (cosine taper; default
% [1/112 1/100 1/20 1/17]). The lines mark the overtone and fundamental-
% mode windows of the inversion (group velocities of the reference model,
% as in sub_setup_records.m): Delta/U_ov (overtone window start), Delta/U_b
% (overtone/fundamental boundary; with opt.windows = 'single' the two are
% one window) and Delta/U_f + 60 s (end). U_start (optional, km/s, scalar
% or [Love Rayleigh]) replaces U_ov as in opt.U_start.
%
% Figures in outdir:
%   record_section_all_<T|Z>.png    all events and stations vs epicentral
%                                   distance, each trace normalized
%   record_section_events_<T|Z>.png one panel per event (16 stations vs
%                                   distance), normalized per event; PA5
%                                   black, step-1 start red
if nargin < 3 || isempty(fc), fc = [1/112 1/100 1/20 1/17]; end
if nargin < 4, U_start = []; end
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
    U = P.(['grv_ref_' ty{wv}]);
    Ub = max(U(:,1)) + 0.2; Uov = max(max(U(:, 2:end))) + 0.6; Uf = min(U(:,1)) - 0.35;
    if ~isempty(U_start), Uov = U_start(min(wv, numel(U_start))); end
    D = zeros(nev, nst); Ut = cell(nev, 1); Up = cell(nev, 1);
    for e = 1:nev
        D(e, :) = gc_dist_az(P.sta.lat, P.sta.lon, P.ev(e).lat*ones(nst, 1), P.ev(e).lon*ones(nst, 1))';
        Ut{e} = bp(P.(['u_true_' ty{wv}]){e});
        if isfield(P, ['u_start_' ty{wv}]), Up{e} = bp(P.(['u_start_' ty{wv}]){e}); end
    end
    % ---- all events in one record section ----
    fig = figure('position', [50 50 1100 1300]); hold on; box on;
    sc = 0.45;                                      % half-amplitude of a trace (deg)
    for e = 1:nev
        for j = 1:nst
            u = Ut{e}(:, j); u = u / max(abs(u));
            plot(t, D(e, j) + sc*u, '-', 'color', ecol(e, :), 'linewidth', 0.5);
        end
    end
    dd = linspace(floor(min(D(:))) - 1, ceil(max(D(:))) + 1, 50);
    hl = plot(dd/km2deg/Uov, dd, 'k--', dd/km2deg/Ub, dd, 'k-', dd/km2deg/Uf + 60, dd, 'k:', 'linewidth', 1.2);
    legend(hl, {sprintf('window start, U_{ov} = %.2f km/s', Uov), ...
                sprintf('overtone/fundamental boundary, U_b = %.2f km/s', Ub), ...
                sprintf('fundamental window end, U_f = %.2f km/s (+60 s)', Uf)}, 'location', 'southeast');
    set(gca, 'ydir', 'normal'); xlim([0 t(end)]); ylim(dd([1 end]));
    xlabel('Time since origin (s)'); ylabel('Epicentral distance (deg)');
    title(sprintf('%s, PA5 full mode sum (branches 0-%d), %g-%g s, traces normalized', ...
        wname{wv}, P.nbranch_data - 1, round(1/fc(4)), round(1/fc(1))));
    print(fig, '-dpng', '-r130', fullfile(outdir, ['record_section_all_' comp{wv} '.png']));
    close(fig);
    % ---- one panel per event ----
    fig = figure('position', [50 50 1900 1250]);
    nr = 3; nc = ceil(nev/nr);
    for e = 1:nev
        ir = ceil(e/nc); icol = e - (ir - 1)*nc;     % tight panel layout
        axes('position', [0.04 + (icol - 1)*0.192, 0.07 + (nr - ir)*0.305, 0.17, 0.255]); hold on; box on;
        de = D(e, :); rng_ = max(de) - min(de);
        sp = max(rng_, 0.5) / nst * 1.6;            % half-amplitude of a trace (deg)
        sn = max(abs(Ut{e}(:)));                    % one normalization per event
        for j = 1:nst
            if ~isempty(Up{e}), plot(t, de(j) + sp*Up{e}(:, j)/sn, '-', 'color', [0.85 0.2 0.2], 'linewidth', 0.5); end
            plot(t, de(j) + sp*Ut{e}(:, j)/sn, 'k-', 'linewidth', 0.6);
        end
        dl = [min(de) - 2*sp, max(de) + 2*sp];
        plot(dl/km2deg/Uov, dl, 'b--', dl/km2deg/Ub, dl, 'b-', dl/km2deg/Uf + 60, dl, 'b:');
        t1 = max(0, min(de)/km2deg/Uov - 300); t2 = min(t(end), max(de)/km2deg/Uf + 400);
        xlim([t1 t2]); ylim(dl);
        title(sprintf('%s  %.0f-%.0f deg', P.ev(e).id, min(de), max(de)), 'fontsize', 8, 'interpreter', 'none');
        if mod(e - 1, nc) == 0, ylabel('\Delta (deg)'); end
        if e > (nr - 1)*nc, xlabel('Time (s)'); end
        set(gca, 'fontsize', 7);
    end
    if exist('OCTAVE_VERSION', 'builtin')
        axes('position', [0 0.965 1 0.03], 'visible', 'off');
        text(0.5, 0.5, sprintf('%s full mode sums, %g-%g s: PA5 (black, branches 0-%d), step-1 start (red); blue: window bounds', ...
            wname{wv}, round(1/fc(4)), round(1/fc(1)), P.nbranch_data - 1), 'horizontalalignment', 'center');
    else
        sgtitle(sprintf('%s full mode sums, %g-%g s: PA5 (black, branches 0-%d), step-1 start (red); blue: window bounds', ...
            wname{wv}, round(1/fc(4)), round(1/fc(1)), P.nbranch_data - 1));
    end
    print(fig, '-dpng', '-r130', fullfile(outdir, ['record_section_events_' comp{wv} '.png']));
    close(fig);
end
end
