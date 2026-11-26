function [figs, W] = plot_posterior_waveforms(rec, cfg, mc, svc, Tdisp, fb, nm, sta, events, pick, outdir, tag, nsub)
% [figs, W] = plot_posterior_waveforms(rec, cfg, mc, svc, Tdisp, fb, nm, sta, events, pick, outdir, tag, nsub)
%
% Posterior probability fields of the predicted waveforms, in the style of
% plot_posterior_pdf.m: for nsub pooled T = 1 samples (default 200, spread
% evenly over the samples of all ensembles) the noise-free waveforms that
% the sample predicts are computed with the inversion's forward operator
% (event_G_spec.m): the sample's dispersion (sv.CL / CR, to the processing
% band with disp_to_band.m), its time shifts tau_n(w) and arrival angles
% dth_n(w) (sv.tauc, sv.dthc: their conditional posterior means given the
% structure, as saved by run_pt_LR.m) and, for each sample, the
% best-fitting spectral amplitudes B_n(w) (the conditional posterior mean
% of event_loglike_ar.m with the sample's noise parameters). At each time
% the histogram of the predictions (smoothed by 1.5 amplitude bins) is
% normalized to its maximum; lines: data, median and 95 % interval of the
% predictions. One figure per event of events (rows: stations pick,
% columns: Love, Rayleigh), fig12_posterior_waveforms_event<e>.png in outdir,
% and the same relative to the median prediction (prediction - median and
% data - median, i.e. the residual of the median prediction, on a scale
% set by them), fig12b_posterior_waveform_spread_event<e>.png: the
% posterior spread of the predictions compared with the misfit.
%
% W (optional output): per event and wave type the predictions at the
% picked stations, W{ie, wv}.t (time), .d (data), .pred (nt x nsub), .sta.
% Called by sub_plot_results.m.
if nargin < 13 || isempty(nsub), nsub = 300; end
nc = numel(svc); nsi = cellfun(@(x) size(x.vsv, 2), svc);
ns = sum(nsi);
js = unique(round(linspace(1, ns, min(nsub, ns))));               % pooled sample indices
ic_ = zeros(1, ns); is_ = zeros(1, ns); j0 = 0;
for ic = 1:nc, ic_(j0 + (1:nsi(ic))) = ic; is_(j0 + (1:nsi(ic))) = 1:nsi(ic); j0 = j0 + nsi(ic); end
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2);
mu_b = repmat([1; 0], nm*K, 1); sb2 = mc.sb^2*ones(2*nm*K, 1);
wn = {'Love', 'Rayleigh'};
cm = pdf_cmap(256); cb_ = [0 0.3 0.95];
sw_ = isfield(mc, 'single_window') && mc.single_window;
figs = []; W = cell(numel(events), 2);
for ie = 1:numel(events)
    e0 = events(ie);
    for wv = 1:2
        r = find([rec.event]' == e0 & [rec.wave]' == wv);
        if isempty(r), continue; end
        r = r(1); R = rec(r);
        stn = ceil(R.idx/cfg.N); tt = mod(R.idx - 1, cfg.N)*cfg.dt + R.tred;
        sel = ismember(stn, pick);
        Pr = zeros(nnz(sel), numel(js));
        for k = 1:numel(js)
            sv = svc{ic_(js(k))}; i = is_(js(k));
            C = {sv.CL(:, :, i), sv.CR(:, :, i)};
            Cf = disp_to_band(C{wv}, Tdisp, fb);
            lin = struct('tau0', reshape(sv.tauc(r, :, 1:Kt, i), nm, Kt), ...
                         'dth0', reshape(sv.dthc(r, :, 1:Ka, i), nm, Ka), 'b0', ones(nm, K));
            G = event_G_spec(R, Cf, lin, cfg, [], [], true);
            [~, b] = event_loglike_ar(R.d, G, R.wt, R.first, sv.sig(i, :), sv.rho(i, :), mu_b, sb2, ...
                                      rec_nlev(R), rec_W(R), rec_wld(R));
            p = G*b; Pr(:, k) = p(sel);
        end
        W{ie, wv} = struct('t', tt(sel), 'd', R.d(sel), 'pred', Pr, 'sta', stn(sel), 'wt', R.wt(sel), 'id', R.id, 'Del0', R.Del0);
    end
    for dev = [false true]
        fig = numbered_figure(1200 + 50*dev + e0, 'position', [50 50 1500 950], 'color', 'w'); figs(end+1) = fig; %#ok<AGROW>
        for wv = 1:2
            if isempty(W{ie, wv}), continue; end
            draw_record(W{ie, wv}, wv, pick, sta, dev, sw_, cm, cb_, cfg.dt, wn);
        end
        if dev, ttl = 'predictions and data minus the median prediction'; fn = 'fig12b_posterior_waveform_spread_event%d.png';
        else, ttl = 'predicted waveforms (noise-free)'; fn = 'fig12_posterior_waveforms_event%d.png'; end
        annotation('textbox', [0 0.955 1 0.04], 'string', sprintf('Run %s: posterior %s, %d of %d pooled samples, %d ensembles', ...
            tag, ttl, numel(js), ns, nc), 'edgecolor', 'none', 'horizontalalignment', 'center', 'fontsize', 11, 'interpreter', 'none');
        if nargin >= 11 && ~isempty(outdir)
            save_figure(fig, fullfile(outdir, sprintf(fn, e0)), 130);
        end
    end
end
end



function draw_record(Wr, wv, pick, sta, dev, sw_, cm, cb_, dt, wn)
% one column (wave type) of the figure: stations pick in rows
tt = Wr.t; stn = Wr.sta; d = Wr.d; Pr = Wr.pred; wt = Wr.wt;
for q = 1:numel(pick)
    ax = subplot(numel(pick), 2, 2*(q - 1) + wv); hold on; box on;
    s = stn == pick(q);
    if ~any(s), axis off; continue; end
    ts = tt(s); ds = d(s); X = Pr(s, :);
    [ts, o] = sort(ts); ds = ds(o); X = X(o, :);
    ws = wt(s); ws = ws(o);
    if dev                                           % relative to the median prediction
        md = prc(X, 50); X = X - md; ds = ds - md;
        ymax = 1.15*max([max(abs([prc(X, 0.5) prc(X, 99.5)]), [], 2); abs(ds)]);
    else
        ymax = 1.15*max([abs(ds); max(abs(X), [], 2)]);
    end
    % regular time grid (gaps between windows left blank)
    tg = (ts(1):dt:ts(end))'; ig = round((ts - ts(1))/dt) + 1;
    nb = 160; bw = 2*ymax/nb; bins = -ymax + bw*((1:nb)' - 0.5);
    H = NaN(nb, numel(tg));
    Hs = zeros(nb, numel(ts));
    for j = 1:numel(ts)
        ib = min(max(floor((X(j, :) + ymax)/bw) + 1, 1), nb);
        Hs(:, j) = accumarray(ib(:), 1, [nb 1]);
    end
    g = exp(-0.5*((-5:5)'/1.5).^2); g = g/sum(g);
    Hs = conv2(g, 1, Hs, 'same'); Hs = Hs ./ max(max(Hs, [], 1), eps);
    H(:, ig) = Hs;
    hi_ = imagesc(tg, bins, H); colormap(ax, cm); caxis([0 1]);
    set(hi_, 'alphadata', double(~isnan(H)));
    set(gca, 'ydir', 'normal', 'layer', 'top', 'color', [0.93 0.93 0.93]);
    % window boundaries
    io = ws == 2*wv - 1; iF = ws == 2*wv;
    for w_ = {io, iF}
        if any(w_{1}), plot(min(ts(w_{1}))*[1 1], ymax*[-1 1], 'k:'); plot(max(ts(w_{1}))*[1 1], ymax*[-1 1], 'k:'); end
    end
    gp = NaN(numel(tg), 1);
    hd = plot(tg, fillg(gp, ig, ds), 'k-', 'linewidth', 1.4);
    hm = plot(tg, fillg(gp, ig, prc(X, 50)), '-', 'color', cb_, 'linewidth', 1.1);
    h95 = plot(tg, fillg(gp, ig, prc(X, 2.5)), '--', 'color', cb_, 'linewidth', 0.8);
    plot(tg, fillg(gp, ig, prc(X, 97.5)), '--', 'color', cb_, 'linewidth', 0.8);
    xlim(tg([1 end])); ylim(ymax*[-1 1]); ylabel(sta{pick(q)});
    if q == 1
        title(sprintf('%s, %s (\\Delta_0 = %.0f km)', wn{wv}, Wr.id, Wr.Del0));
        if wv == 1
            if dev, ll = {'data - median', 'median - median (0)', '95 % of predictions'};
            else, ll = {'data', 'median prediction', '95 % of predictions'}; end
            legend([hd hm h95], ll, ...
                'location', 'northwest', 'fontsize', 7);
        end
    end
    if q == numel(pick)
        if sw_, xlabel('Time since origin (s); dotted: surface-wave window');
        else, xlabel('Time since origin (s); dotted: overtone and fundamental windows'); end
    end
    if q == 1 && wv == 2
        cbh = colorbar; ylabel(cbh, 'probability (normalized at each time)');
    end
end
end

function y = fillg(y, ig, v)
y(ig) = v;
end

function q = prc(X, p)
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
