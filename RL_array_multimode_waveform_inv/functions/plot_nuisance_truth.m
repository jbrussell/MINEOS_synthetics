function txt = plot_nuisance_truth(rec, cfg, mc, pd, svc, stc, Ctrue_f, sig_true, rho_emp, nm, dth_true, outdir, cachef)
% txt = plot_nuisance_truth(rec, cfg, mc, pd, svc, stc, Ctrue_f, sig_true, rho_emp, nm, dth_true, outdir, cachef)
%
% Synthetic test only: true vs recovered nuisance parameters of every record
% and branch (fig9_nuisance_true_vs_recovered.png in outdir), and a text
% summary (txt, for summary.txt). Called by sub_plot_results.m.
%
% "True" values (idagrn data have no exact nuisance parameters, so two
% versions are shown where they exist):
%   tau_n(w) coefficients (Legendre terms of cfg.Psit): least-squares
%       projection of the true path delay pd.true_delay (MINEOS true -
%       reference phase delay at the array centre) onto the tau basis
%       (filled markers, x axis); grey crosses: the values fitted to the
%       data at the TRUE structure (fit_nuisance_LR / spec_refine), which
%       also absorb the delay outside the basis and the amplitude trade-offs.
%   dth_n(w) coefficients (cfg.Psi): constant term = the simulated arrival-
%       angle anomaly of the event (dth_true, deg, from A0 with
%       opt.synth_arrival_shift; 0 otherwise), higher terms 0; grey crosses:
%       fitted at the true structure.
%   B_n(w) (complex, at the K nodes of cfg.Phi): no exact value (it absorbs
%       the excitation and amplitude differences between the true and
%       reference models); x axis = the best fit at the true structure.
% Recovered: tau and dth = posterior mean with 5-95% bars, all chains
% (sv.tauc, sv.dthc); B = best fit (posterior mean of the Gaussian
% marginal) for the final state of each chain (structure and mean tau, dth
% of its last samples, as fig4), averaged over the chains.
% Cycle equivalence: tau_n(w) -> tau_n(w) + k 2 pi/w (a whole cycle at every
% frequency, k integer) leaves exp(-i w tau) and so the waveforms unchanged
% (up to what the tau basis cannot represent, which B absorbs), so the time
% shifts are only defined up to k cycles. For each record and branch the k
% (|k| <= 3) that brings the recovered tau_n(w) closest to the true delay
% over the band is removed from the recovered coefficients (and,
% separately, from the fit at the true structure) before the comparison;
% records with k ~= 0 are drawn as open markers and counted in the summary.
% dth_true : [1 x nrec] (deg) or [] (zeros).
% cachef   : file for the true-structure fit (it takes a few minutes);
%            reused when present.
nrec = numel(rec); nc = numel(svc);
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2);
if isempty(dth_true), dth_true = zeros(1, nrec); end
mcol = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; lines(max(nm-3, 1))];
mk = 'o^';                                    % Love, Rayleigh

% ---- truth: projection of the true path delay onto the tau basis ----
Ttrue = NaN(nrec, nm, Kt);
for r = 1:nrec
    for n = 1:nm
        td = pd.true_delay{r}(:, n); ok = isfinite(td);
        if sum(ok) > Kt, Ttrue(r, n, :) = cfg.Psit(ok, :) \ td(ok); end
    end
end
Atrue = zeros(nrec, nm, Ka); Atrue(:, :, 1) = repmat(dth_true(:), 1, nm);

% ---- fit at the true structure (cached) ----
if ~isempty(cachef) && exist(cachef, 'file') == 2
    C = load(cachef); linc = C.linc;
else
    linc = spec_refine([], rec, Ctrue_f, sig_true, rho_emp, cfg, mc, nm, 30);
    for r = 1:nrec      % final B_n at the converged point (posterior mean)
        linc{r} = nuis_relin(rec(r), Ctrue_f{rec(r).wave}, linc{r}, cfg, sig_true, rho_emp, mc, nm, true);
    end
    if ~isempty(cachef), save(cachef, 'linc', '-v7'); end
end
Tfit = NaN(nrec, nm, Kt); Afit = NaN(nrec, nm, Ka); Bfit = NaN(nrec, nm, K);
for r = 1:nrec
    Tfit(r, :, :) = reshape(linc{r}.tau0, 1, nm, Kt);
    Afit(r, :, :) = reshape(linc{r}.dth0, 1, nm, Ka);
    Bfit(r, :, :) = reshape(linc{r}.b0, 1, nm, K);
end

% ---- posterior ----
Tc = cat(4, svc{1}.tauc); Ac = cat(4, svc{1}.dthc);
for ic = 2:nc, Tc = cat(4, Tc, svc{ic}.tauc); Ac = cat(4, Ac, svc{ic}.dthc); end
pr = @(X, p) prctile_4(X, p);
Tm = mean(Tc, 4); T5 = pr(Tc, 5); T95 = pr(Tc, 95); Ts = std(Tc, 0, 4);
Am = mean(Ac, 4); A5 = pr(Ac, 5); A95 = pr(Ac, 95); As = std(Ac, 0, 4);
% cycle equivalence: remove k whole cycles (k 2 pi/w, projected on the tau basis)
w_ = cfg.w(:); Pc = cfg.Psit \ (2*pi./w_);          % Kt x 1: one cycle on the tau basis
kpost = zeros(nrec, nm); kfit = kpost;
for r = 1:nrec
    for n = 1:nm
        td = pd.true_delay{r}(:, n); ok = isfinite(td);
        if ~any(ok), continue; end
        tp = cfg.Psit * squeeze(Tm(r, n, :)); tf = cfg.Psit * squeeze(Tfit(r, n, :));
        ks = -3:3;
        [~, i1] = min(arrayfun(@(k) sum((tp(ok) - td(ok) - k*2*pi./w_(ok)).^2), ks)); kpost(r, n) = ks(i1);
        [~, i2] = min(arrayfun(@(k) sum((tf(ok) - td(ok) - k*2*pi./w_(ok)).^2), ks)); kfit(r, n) = ks(i2);
    end
end
sh = bsxfun(@times, kpost, reshape(Pc, 1, 1, Kt)); shf = bsxfun(@times, kfit, reshape(Pc, 1, 1, Kt));
Tm = Tm - sh; T5 = T5 - sh; T95 = T95 - sh; Tfit = Tfit - shf;
Bp = zeros(nrec, nm, K);
mu_b = repmat([1; 0], nm*K, 1); sb2 = mc.sb^2*ones(2*nm*K, 1);
for ic = 1:nc
    s1 = stc{ic};
    for r = 1:nrec
        wv = rec(r).wave;
        G1 = event_G_spec(rec(r), s1.Cf{wv}, struct('tau0', reshape(s1.tauc(r,:,1:Kt), nm, Kt), ...
            'dth0', reshape(s1.dthc(r,:,1:Ka), nm, Ka), 'b0', ones(nm, K)), cfg, [], [], true);
        [~, b1] = event_loglike_ar(rec(r).d, G1, rec(r).wt, rec(r).first, s1.sig, s1.rho, mu_b, sb2, ...
                                   rec_nlev(rec(r)), rec_W(rec(r)), rec_wld(rec(r)));
        bb = reshape(b1(1:2:end) + 1i*b1(2:2:end), K, nm).';      % nm x K
        Bp(r, :, :) = Bp(r, :, :) + reshape(bb, 1, nm, K) / nc;
    end
end

% ---- figure ----
np = Kt + Ka + 2; ncol = 4; nrow = ceil(np/ncol);
fig = numbered_figure(1090, 'position', [50 50 1900 480*nrow]);
wave = [rec.wave];
ip = 0;
sets = {};
for a = 1:Kt, sets{end+1} = {'tau', a}; end %#ok<AGROW>
for a = 1:Ka, sets{end+1} = {'dth', a}; end %#ok<AGROW>
tnames = {'constant', 'linear', 'quadratic', 'cubic'};
for q = 1:numel(sets)
    ip = ip + 1; subplot(nrow, ncol, ip); hold on; box on;
    a = sets{q}{2};
    if strcmp(sets{q}{1}, 'tau')
        X = Ttrue(:, :, a); Y = Tm(:, :, a); L = T5(:, :, a); U = T95(:, :, a); F = Tfit(:, :, a);
        un = 's'; ttl = sprintf('\\tau_n: %s term', tnames{a});
    else
        X = Atrue(:, :, a); Y = Am(:, :, a); L = A5(:, :, a); U = A95(:, :, a); F = Afit(:, :, a);
        un = 'deg'; ttl = sprintf('d\\theta_n: %s term', tnames{a});
    end
    hh = zeros(1, nm);
    for n = 1:nm
        for wv = 1:2
            r = find(wave == wv);
            for j = r(:)', plot(X(j, n)*[1 1], [L(j, n) U(j, n)], '-', 'color', 0.4 + 0.6*mcol(n,:)); end
            if strcmp(sets{q}{1}, 'tau'), sk = kpost(r, n) ~= 0; else, sk = false(size(r(:))); end
            h = plot(X(r(~sk), n), Y(r(~sk), n), mk(wv), 'color', mcol(n,:), 'markerfacecolor', mcol(n,:), 'markersize', 5);
            if any(sk), plot(X(r(sk), n), Y(r(sk), n), mk(wv), 'color', mcol(n,:), 'markersize', 8, 'linewidth', 1.5); end
            if wv == 1, hh(n) = h; end
            plot(X(r, n), F(r, n), 'x', 'color', [0.45 0.45 0.45], 'markersize', 6);
        end
    end
    v = [X(:); F(:)]; v = v(isfinite(v));
    lim = max(1.4*max(abs(v)), 0.5*max(abs(Y(:)))); if ~(lim > 0), lim = max(1, max(abs(Y(:)))); end
    lim = min(max(lim, 2), 60);
    plot(lim*[-1 1], lim*[-1 1], 'k--'); xlim(lim*[-1 1]); ylim(lim*[-1 1]); axis square;
    xlabel(sprintf('true (%s)', un)); ylabel(sprintf('posterior mean, 5-95%% (%s)', un)); title(ttl);
    if q == 1
        legend(hh, arrayfun(@(n) sprintf('branch %d', n), 0:nm-1, 'uniformoutput', false), ...
               'location', 'northwest', 'fontsize', 9);
    end
end
% B: amplitude and phase
Xa = abs(Bfit); Ya = abs(Bp); Xp = angle(Bfit)*180/pi; Yp = angle(Bp)*180/pi;
for q = 1:2
    ip = ip + 1; subplot(nrow, ncol, ip); hold on; box on;
    if q == 1, X = Xa; Y = Ya; else, X = Xp; Y = Yp; end
    for n = 1:nm
        for wv = 1:2
            r = wave == wv;
            x = X(r, n, :); y = Y(r, n, :);
            plot(x(:), y(:), mk(wv), 'color', mcol(n,:), 'markerfacecolor', mcol(n,:), 'markersize', 4);
        end
    end
    if q == 1
        m_ = max([X(:); Y(:)]); m_ = min(max(m_, 1.2), 4); plot([0 m_], [0 m_], 'k--'); xlim([0 m_]); ylim([0 m_]);
        xlabel('|B| fitted at the true structure'); ylabel('|B| recovered'); title('|B_n(\omega)| (K nodes)');
    else
        plot([-180 180], [-180 180], 'k--'); xlim([-180 180]); ylim([-180 180]);
        xlabel('arg B fitted at the true structure (deg)'); ylabel('arg B recovered (deg)'); title('arg B_n(\omega)');
    end
    axis square;
end
annotation('textbox', [0 0.955 1 0.045], 'string', ...
    {'Nuisance parameters, true vs recovered (circles: Love, triangles: Rayleigh; grey x: fitted at the true structure)', ...
     'tau truth: MINEOS path delay projected on the tau basis; open markers: recovered tau shifted by whole cycles k 2\pi/\omega (removed)'}, ...
    'edgecolor', 'none', 'horizontalalignment', 'center', 'fontsize', 11);
save_figure(fig, fullfile(outdir, 'fig9_nuisance_true_vs_recovered.png'), 130);

% ---- text summary ----
txt = sprintf(['nuisance parameters vs truth (posterior mean - true; z = that / posterior std; tau truth = projection of the MINEOS path delay;\n', ...
               '  tau after removing whole-cycle shifts k 2 pi/w, which leave the waveforms unchanged)\n']);
for n = 1:nm
    txt = [txt sprintf('  branch %d: records with the recovered tau shifted by whole cycles: %d of %d (k = %s); fit at the true structure: %d\n', ...
        n-1, sum(kpost(:, n) ~= 0), nrec, mat2str(kpost(kpost(:, n) ~= 0, n)'), sum(kfit(:, n) ~= 0))]; %#ok<AGROW>
end
for q = 1:numel(sets)
    a = sets{q}{2};
    if strcmp(sets{q}{1}, 'tau'), E = Tm(:,:,a) - Ttrue(:,:,a); S = Ts(:,:,a); lo = T5(:,:,a); hi = T95(:,:,a); X = Ttrue(:,:,a); nm_ = 'tau'; un = 's';
    else, E = Am(:,:,a) - Atrue(:,:,a); S = As(:,:,a); lo = A5(:,:,a); hi = A95(:,:,a); X = Atrue(:,:,a); nm_ = 'dth'; un = 'deg'; end
    ok = isfinite(E); E(~ok) = 0; Z = E./S; Z(~ok) = 0; nok = max(sum(ok, 1), 1);
    rm = sqrt(sum(E.^2, 1)./nok); cv = 100*sum(ok & X >= lo & X <= hi, 1)./nok; zr = sqrt(sum(Z.^2, 1)./nok);
    txt = [txt sprintf('  %s %-9s (%s): rms error by branch %s | truth in 5-95%% %s %% | rms z %s\n', nm_, tnames{a}, un, ...
        mat2str(round(rm*100)/100), mat2str(round(cv)), mat2str(round(zr*10)/10))]; %#ok<AGROW>
end
ea = abs(Ya) - abs(Bfit); ep = angle(Bp ./ Bfit)*180/pi;
txt = [txt sprintf('  B (vs fit at the true structure): rms |B| difference by branch %s | rms phase difference %s deg\n', ...
    mat2str(round(squeeze(sqrt(mean(mean(ea.^2, 3), 1)))'*100)/100), mat2str(round(squeeze(sqrt(mean(mean(ep.^2, 3), 1)))')))];
end

function P = prctile_4(X, p)
% percentile p along dimension 4 (linear interpolation)
n = size(X, 4); Xs = sort(X, 4);
k = 1 + (n - 1)*p/100; k0 = floor(k); k1 = min(k0 + 1, n); f = k - k0;
P = (1 - f)*Xs(:, :, :, k0) + f*Xs(:, :, :, k1);
end
