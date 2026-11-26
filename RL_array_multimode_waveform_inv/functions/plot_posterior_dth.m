function [fig, As] = plot_posterior_dth(rec, cfg, svc, nm, events, outdir, tag, linc, dth_true)
% [fig, As] = plot_posterior_dth(rec, cfg, svc, nm, events, outdir, tag, linc, dth_true)
%
% Posterior of the arrival-angle anomalies dth_n(w) (deg; event_G_spec.m,
% Legendre terms in log w, cfg.Psi) of every record and mode branch from
% the pooled T = 1 samples (sv.dthc: for each sample the conditional
% posterior mean of the coefficients given the structure, as saved by
% run_pt_LR.m, so the spread across samples leaves out the conditional
% uncertainty of dth itself and is a lower bound).
%
% Figure fig14_arrival_angles.png: columns Love, Rayleigh; rows branches;
% dth_n against period. Thin lines: posterior median of every record
% (colour by event); thick lines with shaded 95 % bands: the records of the
% events in events (as fig4); dotted: the nuisance fit at the true
% structure (linc of plot_nuisance_truth.m, synthetic tests, else []);
% dashed: the true anomaly (dth_true, 1 x nrec, deg; synthetic tests with
% opt.synth_arrival_shift, else []). A second panel row summarizes all
% records: the band-mean (constant term) posterior median and 95 % range of
% every record against the true value (or against record number when the
% truth is unknown).
%
% As : struct saved to <outdir>/dth_posterior<tag>.mat: T, w and per record
%      and branch (nrec x nm x Nf) .med, .lo, .hi (deg), the constant term
%      .c_med, .c_lo, .c_hi (nrec x nm), .true (nrec x 1) and .fit
%      (nrec x nm x Nf, from linc), .id, .wave, .event.
if nargin < 8, linc = []; end
if nargin < 9, dth_true = []; end
nrec = numel(rec); Nf = numel(cfg.w); T = 2*pi./cfg.w(:); Ka = size(cfg.Psi, 2);
Dc = cat(4, svc{1}.dthc); for ic = 2:numel(svc), Dc = cat(4, Dc, svc{ic}.dthc); end   % nrec x nm x Ka x ns
ns = size(Dc, 4);
D = zeros(nrec, nm, Nf, ns);
for r = 1:nrec
    for n = 1:nm
        D(r, n, :, :) = reshape(cfg.Psi * reshape(Dc(r, n, 1:Ka, :), Ka, ns), 1, 1, Nf, ns);
    end
end
q = @(X, p) prc4(X, p);
c = reshape(Dc(:, :, 1, :), nrec, nm, ns);
As = struct('T', T, 'w', cfg.w(:), 'tag', tag, 'nsamples', ns, ...
            'med', q(D, 50), 'lo', q(D, 2.5), 'hi', q(D, 97.5), ...
            'c_med', prc3(c, 50), 'c_lo', prc3(c, 2.5), 'c_hi', prc3(c, 97.5), ...
            'true', NaN(nrec, 1), 'fit', NaN(nrec, nm, Nf), ...
            'id', {{rec.id}}, 'wave', [rec.wave], 'event', [rec.event]);
if ~isempty(dth_true), As.true = dth_true(:); end
if ~isempty(linc)
    for r = 1:nrec, As.fit(r, :, :) = reshape((cfg.Psi * linc{r}.dth0(:, 1:Ka).').', 1, nm, Nf); end
end
if nargin >= 6 && ~isempty(outdir)
    save(fullfile(outdir, ['dth_posterior' tag '.mat']), 'As', '-v7');
end

wn = {'Love', 'Rayleigh'};
ev = [rec.event]; nev = max(ev); ecol = hsv(max(nev, 2))*0.75;
hcol = [0.85 0.1 0.1; 0 0.35 0.9; 0 0.6 0; 0.6 0 0.7];
tk = [15 20 25 30 40 50 60 80 100 120]; tk = tk(tk >= min(T) & tk <= max(T));
fig = numbered_figure(1400, 'position', [50 50 1300 300*(nm + 1)], 'color', 'w');
for wv = 1:2
    rr = find([rec.wave] == wv);
    for n = 1:nm
        subplot(nm + 1, 2, 2*(n - 1) + wv); hold on; box on;
        plot(T([end 1]), [0 0], 'k:');
        if isempty(rr), axis off; continue; end
        for r = rr, plot(T, squeeze(As.med(r, n, :)), '-', 'color', 0.45*ecol(ev(r), :) + 0.55, 'linewidth', 0.6); end
        h = []; l = {};
        for ie = 1:numel(events)
            r = rr(ev(rr) == events(ie));
            if isempty(r), continue; end
            r = r(1); cc = hcol(1 + mod(ie - 1, size(hcol, 1)), :);
            fill([T; flipud(T)], [squeeze(As.lo(r, n, :)); flipud(squeeze(As.hi(r, n, :)))], cc, ...
                 'facealpha', 0.2, 'edgecolor', 'none');
            h(end+1) = plot(T, squeeze(As.med(r, n, :)), '-', 'color', cc, 'linewidth', 2); %#ok<AGROW>
            l{end+1} = sprintf('%s posterior (95%%)', rec(r).id); %#ok<AGROW>
            if any(isfinite(As.fit(r, n, :)))
                h(end+1) = plot(T, squeeze(As.fit(r, n, :)), ':', 'color', cc, 'linewidth', 1.6); %#ok<AGROW>
                l{end+1} = 'fit at true structure'; %#ok<AGROW>
            end
            if isfinite(As.true(r))
                h(end+1) = plot(T([end 1]), As.true(r)*[1 1], '--', 'color', cc, 'linewidth', 1.4); %#ok<AGROW>
                l{end+1} = 'true'; %#ok<AGROW>
            end
        end
        set(gca, 'xscale', 'log', 'xtick', tk, 'xticklabel', arrayfun(@num2str, tk, 'uniformoutput', false), 'fontsize', 10); xlim([min(T) max(T)]);
        ylabel('d\theta_n (deg)'); title(sprintf('%s mode %d: arrival-angle anomaly', wn{wv}, n - 1));
        if n == nm, xlabel('Period (s)'); end
        if n == 1 && ~isempty(h), legend(h, l, 'location', 'best', 'fontsize', 7); end
    end
    % summary row: constant term of every record and branch
    subplot(nm + 1, 2, 2*nm + wv); hold on; box on;
    if isempty(rr), axis off; continue; end
    mk = 'osd^v><ph';
    tr = all(isfinite(As.true(rr)));
    for n = 1:nm
        if tr, x = As.true(rr); else, x = (1:numel(rr))'; end
        x = x + 0.06*(n - (nm + 1)/2)*(tr*0.5 + ~tr);
        cc = min(1, [0 0.3 0.8]*(1 - 0.15*(n - 1)) + 0.1*(n - 1));
        plot([x x]', [As.c_lo(rr, n) As.c_hi(rr, n)]', '-', 'color', cc);   % 95 % range
        plot(x, As.c_med(rr, n), mk(1 + mod(n - 1, numel(mk))), 'color', cc, 'markerfacecolor', cc);
    end
    if tr
        lo_ = As.c_lo(rr, :); hi_ = As.c_hi(rr, :);
        lim = [min([As.true(rr); lo_(:)]) max([As.true(rr); hi_(:)])] + [-0.5 0.5];
        plot(lim, lim, 'k--'); xlim(lim); ylim(lim); axis square;
        xlabel('true d\theta (deg)'); err = As.c_med(rr, :) - As.true(rr);
        title(sprintf('%s: band-mean d\\theta_n, rms error %.2f deg', wn{wv}, sqrt(mean(err(:).^2))));
    else
        xlabel('record'); title(sprintf('%s: band-mean d\\theta_n (median, 95%%)', wn{wv}));
    end
    ylabel('posterior d\theta (deg)');
end
annotation('textbox', [0 0.965 1 0.03], 'string', sprintf(['Run %s: arrival-angle anomalies d\\theta_n(w) ' ...
    '(thin: posterior median of every record; %d pooled samples; bands leave out the conditional uncertainty)'], strrep(tag, '_', '-'), ns), ...
    'edgecolor', 'none', 'horizontalalignment', 'center', 'fontsize', 11);
if nargin >= 6 && ~isempty(outdir)
    save_figure(fig, fullfile(outdir, 'fig14_arrival_angles.png'), 130);
end
end

function Q = prc4(X, p)
X = sort(X, 4); n = size(X, 4);
pos = 1 + (n - 1)*p/100; lo = floor(pos); hi = min(lo + 1, n); fr = pos - lo;
Q = (1 - fr)*X(:, :, :, lo) + fr*X(:, :, :, hi);
end

function Q = prc3(X, p)
X = sort(X, 3); n = size(X, 3);
pos = 1 + (n - 1)*p/100; lo = floor(pos); hi = min(lo + 1, n); fr = pos - lo;
Q = (1 - fr)*X(:, :, lo) + fr*X(:, :, hi);
end
