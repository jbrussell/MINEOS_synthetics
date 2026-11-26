function [figs, Bs] = plot_posterior_B(rec, cfg, mc, svc, Tdisp, fb, nm, events, outdir, tag, linc, nsub)
% [figs, Bs] = plot_posterior_B(rec, cfg, mc, svc, Tdisp, fb, nm, events, outdir, tag, linc, nsub)
%
% Posterior of the spectral corrections B_n(w) (source excitation and the
% part of the path effect common to the array; event_G_spec.m) of every
% record and mode branch. B_n is integrated out in the inversion, so for
% nsub pooled T = 1 samples (default 150, spread evenly over all
% ensembles) the conditional posterior mean of B_n given the sample
% (its dispersion, tau_n(w), dth_n(w) and noise parameters) is computed as
% in plot_posterior_waveforms.m, and B_n(w) = sum_k beta_nk Phi_k(w)
% evaluated on the frequency grid of the records (cfg.w).
%
% One figure per wave type (fig13_B_<Love|Rayleigh>.png): rows = branches,
% columns = |B_n(w)| and arg B_n(w) (deg) against period. Thin lines: the
% posterior median of every record (colour by event); thick lines with
% shaded 95 % bands: the records of the events in events (as fig4); dashed
% (synthetic tests, linc = the nuisance fit at the true structure of
% plot_nuisance_truth.m, else []): B_n of that fit for the same records.
% The phase of each sample is unwrapped along frequency and shifted by
% whole turns so that it lies in (-180, 180] at the band centre.
%
% Bs : struct saved to <outdir>/B_posterior<tag>.mat: T (s), w, and per
%      record r and branch n (nrec x nm x Nf) the median, 2.5 and 97.5
%      percentiles of |B| (.amp, .amp_lo, .amp_hi) and of arg B in degrees
%      (.ph, .ph_lo, .ph_hi), .ph_true/.amp_true from linc (NaN otherwise),
%      .id, .wave, .event.
if nargin < 11, linc = []; end
if nargin < 12 || isempty(nsub), nsub = 150; end
nc = numel(svc); nsi = cellfun(@(x) size(x.vsv, 2), svc); ns = sum(nsi);
js = unique(round(linspace(1, ns, min(nsub, ns))));
ic_ = zeros(1, ns); is_ = zeros(1, ns); j0 = 0;
for ic = 1:nc, ic_(j0 + (1:nsi(ic))) = ic; is_(j0 + (1:nsi(ic))) = 1:nsi(ic); j0 = j0 + nsi(ic); end
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2);
mu_b = repmat([1; 0], nm*K, 1); sb2 = mc.sb^2*ones(2*nm*K, 1);
nrec = numel(rec); Nf = numel(cfg.w); T = 2*pi./cfg.w(:);
[~, imid] = min(abs(log(cfg.w) - mean(log(cfg.w([1 end])))));
Bc = complex(zeros(nrec, nm, Nf, numel(js)));
for k = 1:numel(js)
    sv = svc{ic_(js(k))}; i = is_(js(k));
    C = {sv.CL(:, :, i), sv.CR(:, :, i)};
    Cfs = {disp_to_band(C{1}, Tdisp, fb), disp_to_band(C{2}, Tdisp, fb)};
    for r = 1:nrec
        R = rec(r);
        lin = struct('tau0', reshape(sv.tauc(r, :, 1:Kt, i), nm, Kt), ...
                     'dth0', reshape(sv.dthc(r, :, 1:Ka, i), nm, Ka), 'b0', ones(nm, K));
        G = event_G_spec(R, Cfs{R.wave}, lin, cfg, [], [], true);
        [~, b] = event_loglike_ar(R.d, G, R.wt, R.first, sv.sig(i, :), sv.rho(i, :), mu_b, sb2, ...
                                  rec_nlev(R), rec_W(R), rec_wld(R));
        bb = reshape(b(1:2:end) + 1i*b(2:2:end), K, nm);        % K x nm
        Bc(r, :, :, k) = reshape((cfg.Phi * bb).', 1, nm, Nf);
    end
end
% amplitude and unwrapped phase (deg) statistics
A = abs(Bc);
Ph = angle(Bc);
Ph = unwrap(Ph, [], 3);
Ph = Ph - 2*pi*round(Ph(:, :, imid, :)/(2*pi));
Ph = Ph*180/pi;
q = @(X, p) prc4(X, p);
Bs = struct('T', T, 'w', cfg.w(:), 'tag', tag, 'nsamples', numel(js), ...
            'amp', q(A, 50), 'amp_lo', q(A, 2.5), 'amp_hi', q(A, 97.5), ...
            'ph', q(Ph, 50), 'ph_lo', q(Ph, 2.5), 'ph_hi', q(Ph, 97.5), ...
            'amp_true', NaN(nrec, nm, Nf), 'ph_true', NaN(nrec, nm, Nf), ...
            'id', {{rec.id}}, 'wave', [rec.wave], 'event', [rec.event]);
if ~isempty(linc)
    for r = 1:nrec
        Bt = (cfg.Phi * linc{r}.b0.').';                         % nm x Nf
        pt = unwrap(angle(Bt), [], 2); pt = pt - 2*pi*round(pt(:, imid)/(2*pi));
        Bs.amp_true(r, :, :) = reshape(abs(Bt), 1, nm, Nf);
        Bs.ph_true(r, :, :) = reshape(pt*180/pi, 1, nm, Nf);
    end
end
if nargin >= 9 && ~isempty(outdir)
    save(fullfile(outdir, ['B_posterior' tag '.mat']), 'Bs', '-v7');
end

% figures
wn = {'Love', 'Rayleigh'};
ev = [rec.event]; nev = max(ev); ecol = hsv(max(nev, 2))*0.75;
hcol = [0.85 0.1 0.1; 0 0.35 0.9; 0 0.6 0; 0.6 0 0.7];
figs = [];
for wv = 1:2
    rr = find([rec.wave] == wv);
    if isempty(rr), continue; end
    fig = numbered_figure(1300 + wv, 'position', [50 50 1300 max(560, 320*nm)], 'color', 'w'); figs(end+1) = fig; %#ok<AGROW>
    for n = 1:nm
        for c = 1:2
            subplot(nm, 2, 2*(n - 1) + c); hold on; box on;
            if c == 1, M = Bs.amp; L = Bs.amp_lo; U = Bs.amp_hi; Tr = Bs.amp_true; yl = '|B_n|'; y0 = 1;
            else,      M = Bs.ph;  L = Bs.ph_lo;  U = Bs.ph_hi;  Tr = Bs.ph_true;  yl = 'arg B_n (deg)'; y0 = 0; end
            plot(T([end 1]), y0*[1 1], 'k:');
            for r = rr
                plot(T, squeeze(M(r, n, :)), '-', 'color', 0.45*ecol(ev(r), :) + 0.55, 'linewidth', 0.6);
            end
            h = []; l = {};
            for ie = 1:numel(events)
                r = rr(ev(rr) == events(ie));
                if isempty(r), continue; end
                r = r(1); cc = hcol(1 + mod(ie - 1, size(hcol, 1)), :);
                fill([T; flipud(T)], [squeeze(L(r, n, :)); flipud(squeeze(U(r, n, :)))], cc, ...
                     'facealpha', 0.2, 'edgecolor', 'none');
                h(end+1) = plot(T, squeeze(M(r, n, :)), '-', 'color', cc, 'linewidth', 2); %#ok<AGROW>
                l{end+1} = sprintf('%s (posterior, 95%%)', rec(r).id); %#ok<AGROW>
                if any(isfinite(Tr(r, n, :)))
                    h(end+1) = plot(T, squeeze(Tr(r, n, :)), '--', 'color', cc, 'linewidth', 1.4); %#ok<AGROW>
                    l{end+1} = sprintf('%s (fit at true structure)', rec(r).id); %#ok<AGROW>
                end
            end
            set(gca, 'xscale', 'log', 'fontsize', 10);
            xlim([min(T) max(T)]);
            tk = [15 20 25 30 40 50 60 80 100 120]; tk = tk(tk >= min(T) & tk <= max(T));
            set(gca, 'xtick', tk, 'xticklabel', arrayfun(@num2str, tk, 'uniformoutput', false));
            ylabel(yl); title(sprintf('%s mode %d: %s', wn{wv}, n - 1, yl));
            if n == nm, xlabel('Period (s)'); end
            if n == 1 && c == 1 && ~isempty(h), legend(h, l, 'location', 'best', 'fontsize', 7); end
        end
    end
    annotation('textbox', [0 0.965 1 0.03], 'string', sprintf(['Run %s: spectral corrections B_n(w), %s ' ...
        '(thin: posterior median of every record; %d of %d pooled samples)'], tag, wn{wv}, numel(js), ns), ...
        'edgecolor', 'none', 'horizontalalignment', 'center', 'fontsize', 10, 'interpreter', 'none');
    if nargin >= 9 && ~isempty(outdir)
        save_figure(fig, fullfile(outdir, sprintf('fig13_B_%s.png', wn{wv})), 130);
    end
end
end

function Q = prc4(X, p)
% percentile along dimension 4
X = sort(X, 4); n = size(X, 4);
pos = 1 + (n - 1)*p/100;
lo = floor(pos); hi = min(lo + 1, n); fr = pos - lo;
Q = (1 - fr)*X(:, :, :, lo) + fr*X(:, :, :, hi);
end
