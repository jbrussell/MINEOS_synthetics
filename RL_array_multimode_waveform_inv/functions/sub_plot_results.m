% SUB_PLOT_RESULTS  (drivee script, called by D4_step2_plot_results.m)
% Combine the chains (or PT ensembles) of a run and make the figures and
% summary.txt in its outdir. Set TAG ('_plain' or '_pt') before calling.
% True model: PA5; reference model: refp (the step-1 isotropic model).
if ~exist('tag', 'var'), tag = '_plain'; end
load(['setup' tag '.mat']);
% reference updated at the end of burn-in (D3, opt.update_ref = 'burnin';
% with 'lin' only cfg.lin, the dispersion linearization, is replaced):
% the post-burn-in samples use the updated reference spectra and MINEOS
% linearization, so the final waveforms, the path delays and log L(true)
% are evaluated with it; the "reference" waveforms and curves stay those of
% the step-1 model
upd = [];
if exist(['refupd' tag '.mat'], 'file') == 2
    upd = load(['refupd' tag '.mat']); upd = upd.U;
    upd_plot = ~(isfield(upd, 'reverted') && upd.reverted);
    upd_lin = isfield(upd, 'mode') && strcmp(upd.mode, 'lin');
    if upd_plot && upd_lin       % linearization only: records and path delays unchanged
        Pu_ = load(upd.prepfile); Pu_ = Pu_.P;
        cfg = apply_linearization(cfg, Pu_); clear Pu_
    elseif upd_plot              % (a rejected update leaves the setup's reference)
        Pu_ = load(upd.prepfile); Pu_ = Pu_.P;
        refsyn0_ = {rec.refsyn};
        [rec, cfg, Cu_f_] = apply_reference(rec, cfg, Pu_, Tdisp, fb);
        [rec.refsyn] = refsyn0_{:};
        Ctrue_f_ = {disp_to_band(CtrueL, Tdisp, fb), disp_to_band(CtrueR, Tdisp, fb)};
        for r = 1:numel(rec)
            pd.true_delay{r} = rec(r).Del0 * (1./Ctrue_f_{rec(r).wave} - 1./Cu_f_{rec(r).wave});
        end
        clear Pu_ refsyn0_ Cu_f_ Ctrue_f_
    end
    Ltrue_setup = Ltrue; Ltrue = upd.Ltrue;
    upd.xi = (upd.vsh./upd.vsv).^2;
end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
nc = mc.nchains; nrec = numel(rec); wave = [rec.wave]';
svc = cell(nc, 1); stc = cell(nc, 1); infc = cell(nc, 1);
for ic = 1:nc
    D = load(sprintf('chain%s_%02d.mat', tag, ic)); svc{ic} = D.sv; stc{ic} = D.st;
    if isfield(D, 'info'), infc{ic} = D.info; end
end
cat2 = @(f) cell2mat(cellfun(@(x) x.(f), svc', 'uniformoutput', false));
cat3 = @(f) cat3c(cellfun(@(x) x.(f), svc, 'uniformoutput', false));
V = struct();                    % (V may be left over from D3)
V.vsv = cat2('vsv'); V.vsh = cat2('vsh'); V.xi = cat2('xi');
CLs = cat3('CL'); CRs = cat3('CR');
ns = size(V.vsv, 2);
ccol = lines(max(nc, 7));
mcol = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; lines(max(nm-3, 1))];
in = zplot >= cfg.moho & zplot < cfg.zmax;
names = {'vsv', 'vsh', 'xi'}; labs = {'V_{SV} (km/s)', 'V_{SH} (km/s)', '\xi = (V_{SH}/V_{SV})^2'};
truv = {tru.vsv, tru.vsh, tru.xi};
refv = {refp.vsv, refp.vsh, refp.xi};
eta_ = isfield(mc, 'eta_free') && mc.eta_free;   % eta a model parameter
sw_ = isfield(mc, 'single_window') && mc.single_window;   % one surface-wave window per wave type
if eta_
    V.eta = cat2('eta'); names{4} = 'eta'; labs{4} = '\eta'; truv{4} = tru.eta; refv{4} = refp.eta;
end

% ---------------- summary statistics ----------------
fid = fopen(fullfile(outdir, 'summary.txt'), 'w');
for fo = [1 fid]
    fprintf(fo, 'Run %s: %d chains/ensembles, %d saved samples\n', tag, nc, ns);
    if ~isempty(upd) && upd_plot && upd_lin
        fprintf(fo, 'MINEOS linearization updated at the end of burn-in (%s, mean of %d late burn-in samples; waveform reference unchanged)\n', ...
            upd.refid, upd.nsamples);
        if isfield(upd, 'stages')                  % several updates during the burn-in
            for k_ = 1:numel(upd.stages)
                u_ = upd.stages{k_};
                fprintf(fo, '  update %d (%s): states'' dispersion change rms %s %%, max %s %%; logL change %s\n', k_, u_.refid, ...
                    mat2str(round(100*u_.dCrms')/100), mat2str(round(100*u_.dCmax')/100), mat2str(round(u_.dL')));
            end
        else
            fprintf(fo, '  end-of-burn-in states: dispersion change rms %s %%, max %s %%; logL change %s\n', ...
                mat2str(round(100*upd.dCrms')/100), mat2str(round(100*upd.dCmax')/100), mat2str(round(upd.dL')));
        end
    elseif ~isempty(upd) && upd_plot
        fprintf(fo, 'reference updated at the end of burn-in (%s, mean of %d late burn-in samples); log L(true) %.1f -> %.1f\n', ...
            upd.refid, upd.nsamples, Ltrue_setup, Ltrue);
    elseif ~isempty(upd)
        fprintf(fo, 'reference update at the end of burn-in REJECTED (end-of-burn-in states fit worse, dL = %s): reference not updated\n', ...
            mat2str(round(upd.dL')));
    end
end
Rz = struct();
for q = 1:numel(names)
    X = V.(names{q});
    R = NaN(numel(zplot), 1);
    for iz = find(in)'
        R(iz) = gr_rhat(cell2mat(cellfun(@(x) x.(names{q})(iz,:)', svc', 'uniformoutput', false)));
    end
    Rz.(names{q}) = R;
    lo = prctile_(X, 5); hi = prctile_(X, 95);
    cov = mean(truv{q}(in) >= lo(in) & truv{q}(in) <= hi(in));
    err = mean(abs(mean(X(in,:), 2) - truv{q}(in)));
    M = cell2mat(cellfun(@(x) mean(x.(names{q}), 2), svc', 'uniformoutput', false));
    spr = mean(std(M(in,:), 0, 2));
    ok = ~isnan(R); if ~any(ok), ok(1) = true; R(1) = NaN; end
    for fo = [1 fid]
        fprintf(fo, '%-4s: R-hat median %.2f max %.2f | truth in 5-95%%: %3.0f%% | mean|err| %.3f | chain-mean spread %.3f | 5-95 width %.3f\n', ...
            names{q}, median(R(ok)), max(R(ok)), 100*cov, err, spr, mean(hi(in) - lo(in)));
    end
end
wl = 'LR'; wfull = {'Love', 'Rayleigh'};
for wv = 1:2
    rc = NaN(1, nm);
    for n = 1:nm
        rr = [];
        for ip = find(cfg.need(:, n))'
            fld = ['C' wl(wv)];
            rr(end+1) = gr_rhat(cell2mat(cellfun(@(x) squeeze(x.(fld)(ip,n,:)), svc', 'uniformoutput', false))); %#ok<SAGROW>
        end
        rc(n) = max(rr);
    end
    for fo = [1 fid], fprintf(fo, '%s phase velocity R-hat max per mode: %s\n', ...
            wfull{wv}, mat2str(round(rc*100)/100)); end
end
Lp = cellfun(@(x) median(x.L(mc.burn+1:end)), svc);
for fo = [1 fid]
    fprintf(fo, 'posterior median logL per chain: %s | true model %.1f\n', mat2str(round(Lp')), Ltrue);
    kx = cell2mat(cellfun(@(x) x.kx, svc, 'uniformoutput', false));
    kv = cell2mat(cellfun(@(x) x.kv, svc, 'uniformoutput', false));
    fprintf(fo, 'interior control points: Vsv median %g (5-95%%: %g-%g); xi median %g (5-95%%: %g-%g)\n', ...
        median(kv), round(prctile_(kv', 5)), round(prctile_(kv', 95)), median(kx), round(prctile_(kx', 5)), round(prctile_(kx', 95)));
    pr_ = struct('smooth_sig', [0 0], 'smooth_param', 'vsv_vsh', 'smooth_len', 10, 'vsv_ref_sig', 0, 'vsv_ref_len', 10, ...
                 'ximin', NaN, 'ximax', NaN, 'xi_ref_sig', 0, 'xi_ref_mean', 1, 'xi_ref_len', 10, 'xi_ref_hier', false, ...
                 'smooth_hier', false);
    for f_ = fieldnames(pr_)', if isfield(mc, f_{1}) && ~isempty(mc.(f_{1})), pr_.(f_{1}) = mc.(f_{1}); end, end
    if ~isnumeric(pr_.vsv_ref_sig), pr_.vsv_ref_sig = -1; end   % (-1 printed: depth-dependent, 'step1')
    pf_ = @(f, varargin) pr_.(f);
    fprintf(fo, 'priors: smoothness smooth_sig %s (%s, per %g km); reference Vsv vsv_ref_sig %g km/s (per %g km; 0 = off)\n', ...
        mat2str(pf_('smooth_sig', [0 0])), pf_('smooth_param', 'vsv_vsh'), pf_('smooth_len', 10), pf_('vsv_ref_sig', 0), pf_('vsv_ref_len', 10));
    fprintf(fo, 'priors: xi uniform %g-%g; xi toward %g with xi_ref_sig %g (per %g km; 0 = off)\n', ...
        pf_('ximin'), pf_('ximax'), pf_('xi_ref_mean'), pf_('xi_ref_sig'), pf_('xi_ref_len'));
    if pr_.smooth_hier && isfield(svc{1}, 'ss')
        ss_ = cell2mat(cellfun(@(x) x.ss, svc, 'uniformoutput', false)); ss_ = ss_(all(isfinite(ss_), 2), :);
        if ~isempty(ss_), fprintf(fo, 'priors: smoothing strengths sampled (hierarchical): median %s\n', mat2str(median(ss_, 1), 3)); end
    end
    if pr_.xi_ref_hier && isfield(svc{1}, 'xs')
        xs_ = cell2mat(cellfun(@(x) x.xs(:), svc, 'uniformoutput', false)); xs_ = sort(xs_(isfinite(xs_)));
        if ~isempty(xs_), fprintf(fo, 'priors: xi strength sampled (hierarchical): sigma_xi median %.3f [5-95%%: %.3f-%.3f]\n', ...
            median(xs_), xs_(max(1, round(0.05*numel(xs_)))), xs_(round(0.95*numel(xs_)))); end
    end
end
R4 = cell2mat(cellfun(@(x) x.rho, svc, 'uniformoutput', false));
S4 = cell2mat(cellfun(@(x) x.sig, svc, 'uniformoutput', false));
for fo = [1 fid]
    if sw_, fprintf(fo, 'single surface-wave window per wave type: noise slots (L, L, R, R) tied\n'); end
    if isfield(mc, 'fixed_corr') && mc.fixed_corr
        fprintf(fo, 'noise correlation: fixed %s model from the pre-P noise (not sampled) | empirical lag-1 %s\n', ...
            upper(mc.noise_model), mat2str(round(rho_emp*100)/100));
    else
        fprintf(fo, 'noise correlation r (L ovt, L fund, R ovt, R fund): posterior %s +- %s | empirical %s\n', ...
            mat2str(round(mean(R4)*100)/100), mat2str(round(std(R4)*100)/100), mat2str(round(rho_emp*100)/100));
    end
    fprintf(fo, 'noise sigma / sigma_true: %s\n', mat2str(round(mean(S4)./sig_true*100)/100));
end
% band averages of xi: posterior mean, 5-95%, truth, reference
bands = [7 60; 60 160; 160 250; 250 350];
for fo = [1 fid]
    fprintf(fo, 'xi band averages (depth below seafloor): posterior mean [5-95%%] | true | reference\n');
    for b = 1:size(bands, 1)
        ib = zplot >= bands(b,1) & zplot < bands(b,2);
        xs = mean(V.xi(ib, :), 1);
        fprintf(fo, '  %3d-%3d km: %.3f [%.3f-%.3f] | %.3f | %.3f\n', bands(b,:), mean(xs), ...
            prctile_(xs, 5), prctile_(xs, 95), mean(tru.xi(ib)), mean(refp.xi(ib)));
    end
end
fclose(fid);

%% 1) posteriors of Vsv, Vsh, xi
fig = figure('position', [50 50 1500 650]);
ranges = {3.5:0.02:5.3, 3.5:0.02:5.3, 0.88:0.005:1.22};
for q = 1:3
    subplot(1, 3, q); hold on; box on;
    X = V.(names{q}); vb = ranges{q}; dv = vb(2) - vb(1);
    H = zeros(numel(zplot), numel(vb));
    for iz = 1:numel(zplot)
        ib = min(max(round((X(iz,:) - vb(1))/dv) + 1, 1), numel(vb));
        H(iz,:) = accumarray(ib(:), 1, [numel(vb) 1])' / ns;
    end
    imagesc(vb, zplot, H); colormap(flipud(gray)); caxis([0 0.25]); set(gca, 'ydir', 'reverse');
    h1 = plot(truv{q}, zplot, 'r-', 'linewidth', 2);
    h4 = plot(refv{q}, zplot, '-', 'color', [0 0.6 0], 'linewidth', 1.2);
    h2 = plot(mean(X, 2), zplot, 'b-', 'linewidth', 1.5);
    h3 = plot(prctile_(X, 5), zplot, 'b--', prctile_(X, 95), zplot, 'b--');
    plot(vb([1 end]), cfg.zmax*[1 1], 'k:');
    xlim(vb([1 end])); ylim([0 400]); xlabel(labs{q}); ylabel('Depth (km)');
    title(sprintf('Posterior %s', labs{q}));
    hh = [h1 h4 h2 h3(1)]; ll = {'true (PA5)', 'reference', 'posterior mean', '5-95%'};
    if ~isempty(upd) && upd_plot
        im_ = zplot >= cfg.zmin & zplot <= cfg.zmax;
        hh(end+1) = plot(upd.(names{q})(im_), zplot(im_), '-.', 'color', [0.8 0.4 0], 'linewidth', 1.2); %#ok<SAGROW>
        if upd_lin, ll{end+1} = 'linearization model (end of burn-in)'; %#ok<SAGROW>
        else, ll{end+1} = 'updated reference (end of burn-in)'; end %#ok<SAGROW>
    end
    if q == 1, legend(hh, ll, 'location', 'southwest'); end
end
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig1_posteriors.png'));

%% final posterior: pooled T = 1 samples of all ensembles and their statistics
% (posterior_summary<tag>.mat: mean, median, 68/95 % intervals, R-hat,
% control points, dispersion; write_posterior_summary.m) and the
% probability fields of Vsv, Vsh, xi (fig10_posterior_pdf.png) and of the
% Love and Rayleigh dispersion of all branches (fig11_posterior_dispersion.png,
% relative to the reference; fig11b_posterior_dispersion_abs.png, absolute)
Psum = write_posterior_summary(svc, zplot, cfg, mc, tru, refp, Tdisp, tag, ...
                               fullfile(outdir, ['posterior_summary' tag '.mat']), ...
                               struct('trueL', CtrueL, 'trueR', CtrueR, 'refL', CrefL, 'refR', CrefR));
fig = plot_posterior_pdf(Psum, fullfile(outdir, 'fig10_posterior_pdf.png'));
if fig_close_mode, close(fig); end
fig = plot_posterior_disp(Psum, fullfile(outdir, 'fig11_posterior_dispersion.png'));
if fig_close_mode, close(fig); end
fig = plot_posterior_disp(Psum, fullfile(outdir, 'fig11b_posterior_dispersion_abs.png'), [], 'absolute');
if fig_close_mode, close(fig); end
% posterior predicted waveforms (fig12_posterior_waveforms_event<e>.png) and
% their spread relative to the median prediction, with the data residual
% (fig12b_posterior_waveform_spread_event<e>.png): events and stations of fig4
try
    pk_ = [1 6 11 16]; pk_ = pk_(pk_ <= numel(sta));
    figs_ = plot_posterior_waveforms(rec, cfg, mc, svc, Tdisp, fb, nm, sta, ...
        unique([1 min(8, max([rec.event]))]), pk_, outdir, tag);
    if fig_close_mode, close(figs_); end
    clear figs_ pk_
catch err_
    fprintf('posterior waveform figures not made: %s\n', err_.message);
end
fid = fopen(fullfile(outdir, 'summary.txt'), 'a');
for fo = [1 fid]
    fprintf(fo, 'final posterior (all ensembles pooled, %d samples): %s\n', Psum.ns, ...
        fullfile(outdir, ['posterior_summary' tag '.mat']));
    fprintf(fo, 'interior control points (pooled): Vsv mean %.1f median %g, 95%% %g-%g, R-hat %.2f; xi mean %.1f median %g, 95%% %g-%g, R-hat %.2f\n', ...
        Psum.ncp.vsv.mean, Psum.ncp.vsv.median, Psum.ncp.vsv.ci95, Psum.ncp.vsv.rhat, ...
        Psum.ncp.xi.mean, Psum.ncp.xi.median, Psum.ncp.xi.ci95, Psum.ncp.xi.rhat);
end
fclose(fid);

%% 2) per-chain means and starts, model dimensions
fig = figure('position', [50 50 1900 700]);
xr = {[3.9 5.1], [3.9 5.4], [0.88 1.22]};
for q = 1:3
    subplot(1, 5, q); hold on; box on; hc = zeros(nc, 1);
    for ic = 1:nc
        plot(svc{ic}.st0.prof(:, q), zplot, ':', 'color', ccol(ic,:));
        hc(ic) = plot(mean(svc{ic}.(names{q}), 2), zplot, '-', 'color', ccol(ic,:), 'linewidth', 1.3);
    end
    plot(truv{q}, zplot, 'r-', 'linewidth', 2);
    set(gca, 'ydir', 'reverse'); ylim([0 400]); xlim(xr{q}); xlabel(labs{q});
    if q == 1, ylabel('Depth (km)'); title('Chain means (solid), starts (dotted)'); end
end
for q = 1:2
    subplot(1, 5, 3 + q); hold on; box on;
    f = {'kv', 'kx'}; kmax = [mc.kvmax mc.kxmax];
    for ic = 1:nc
        cnt = accumarray(svc{ic}.(f{q}) + 1, 1, [kmax(q)+1 1]) / numel(svc{ic}.(f{q}));
        plot(0:kmax(q), cnt, '-o', 'color', ccol(ic,:), 'markerfacecolor', ccol(ic,:), 'markersize', 3);
    end
    xlim([-0.5 kmax(q)+0.5]); xlabel('interior control points');
    tl = {'V_{SV} control points', '\xi control points (k_\xi \leq k_V)'}; title(tl{q});
end
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig2_chains_dimensions.png'));

%% 3) dispersion fits: Love and Rayleigh, absolute and minus truth
fig = figure('position', [50 50 max(1200, 520*nm) 1500]);   % (room for 15-pt axes)
Csets = {CLs, CRs}; Ct = {CtrueL, CtrueR}; Cr = {CrefL, CrefR}; wn = {'Love', 'Rayleigh'};
for wv = 1:2
    for n = 1:nm
        ip = find(cfg.need(:, n)); Tn = Tdisp(ip);
        Cs = squeeze(Csets{wv}(ip, n, :));
        lo = prctile_(Cs, 2.5); hi = prctile_(Cs, 97.5); mn = mean(Cs, 2);
        for row = 1:2
            subplot(4, nm, (2*(wv-1) + row - 1)*nm + n); hold on; box on;
            if row == 1, off = zeros(numel(ip), 1); else, off = Ct{wv}(ip, n); end
            fill([Tn fliplr(Tn)], [lo' - off' fliplr(hi' - off')], mcol(n,:), ...
                'facealpha', 0.3, 'edgecolor', 'none');
            hm = plot(Tn, mn - off, '-', 'color', mcol(n,:), 'linewidth', 1.5);
            ht = plot(Tn, Ct{wv}(ip, n) - off, 'ko', 'markerfacecolor', 'k', 'markersize', 3);
            if row == 1
                hr = plot(Tn, Cr{wv}(ip, n), '--', 'color', [0.5 0.5 0.5]);
                ylabel('c (km/s)'); title(sprintf('%s mode %d', wn{wv}, n-1));
                if n == 1 && wv == 1
                    legend([hm ht hr], {'posterior mean (95%)', 'true (PA5)', 'reference'}, ...
                        'location', 'northwest', 'fontsize', 7);
                end
            else
                plot(Tn([1 end]), [0 0], 'k-'); ylabel('c - c_{true} (km/s)');
                title(sprintf('%s mode %d: misfit', wn{wv}, n-1));
                if wv == 2, xlabel('Period (s)'); end
            end
            xlim([Tdisp(1) Tdisp(end)]);
        end
    end
end
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig3_dispersion_fits.png'));

%% 4) waveform fits: data, step-1 model (reference) and final model, with windows
% step 1: sum of the reference branches 0..nm-1 of the step-1 model
% (idagrn6_mask, no nuisance correction); final: the waveform operator for the final state
% of chain/ensemble 1 (structure, posterior-mean tau_n(w) and dth_n(w) of
% the last samples, and the best-fitting B_n(w)). All processed like the data (band-pass, record normalization).
pick = [1 6 11 16];
for e0 = unique([1 min(8, max([rec.event]))])
    fig = figure('position', [50 50 1500 950]);
    for wv = 1:2
        r = find([rec.event]' == e0 & wave == wv);
        if isempty(r), continue; end
        s1 = stc{1};
        K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2);
        mu_b = repmat([1; 0], nm*K, 1); sb2 = mc.sb^2*ones(2*nm*K, 1);
        G1 = event_G_spec(rec(r), s1.Cf{wv}, struct('tau0', reshape(s1.tauc(r,:,1:Kt), nm, Kt), ...
            'dth0', reshape(s1.dthc(r,:,1:Ka), nm, Ka), 'b0', ones(nm, K)), cfg, [], [], true);
        [~, b1] = event_loglike_ar(rec(r).d, G1, rec(r).wt, rec(r).first, s1.sig, s1.rho, mu_b, sb2, rec_nlev(rec(r)), rec_W(rec(r)), rec_wld(rec(r)));
        p1 = G1*b1;
        stn = ceil(rec(r).idx/cfg.N); tt = mod(rec(r).idx - 1, cfg.N)*cfg.dt + rec(r).tred;
        for q = 1:numel(pick)
            subplot(numel(pick), 2, 2*(q-1) + wv); hold on; box on;
            s = stn == pick(q); ymax = 1.15*max(abs(rec(r).d(s)));
            io = s & rec(r).wt == 2*wv-1; iF = s & rec(r).wt == 2*wv;
            if any(io), patch([min(tt(io)) max(tt(io)) max(tt(io)) min(tt(io))], ymax*[-1 -1 1 1], ...
                    [1 0.88 0.75], 'edgecolor', 'none'); end
            if any(iF), patch([min(tt(iF)) max(tt(iF)) max(tt(iF)) min(tt(iF))], ymax*[-1 -1 1 1], ...
                    [0.85 0.92 1], 'edgecolor', 'none'); end
            hd = plot(tt(s), rec(r).d(s), 'k-', 'linewidth', 1.5);
            hr = plot(tt(s), rec(r).refsyn(s), '-', 'color', [0.3 0.3 1]);
            hf = plot(tt(s), p1(s), '-', 'color', [0.85 0.1 0.1], 'linewidth', 1);
            ylim(ymax*[-1 1]); ylabel(sta{pick(q)});
            if q == 1
                title(sprintf('%s, %s (\\Delta_0 = %.0f km)', wn{wv}, rec(r).id, rec(r).Del0));
                if wv == 1
                    legend([hd hr hf], {'data', 'step 1 (reference)', 'final model'}, ...
                        'location', 'northwest', 'fontsize', 7);
                end
            end
            if q == numel(pick)
                if sw_, xlabel('Time since origin (s); shaded: surface-wave window');
                else, xlabel('Time since origin (s); orange: overtone window, blue: fundamental window'); end
            end
        end
    end
    print(fig, '-dpng', '-r130', fullfile(outdir, sprintf('fig4_waveform_fits_event%d.png', e0)));
end

%% 5) path delays and time shifts
% left: path delay of the data relative to the reference spectra
% (MINEOS true - reference phase velocity over the distance to the array
% center), which the nuisance model (tau_n and B_n(w)) has to absorb;
% right: posterior tau_n vs that delay at 40 s
fig = figure('position', [50 50 1500 600]);
subplot(1, 2, 1); hold on; box on;
for r = 1:nrec
    for n = 1:nm, plot(1./fb, pd.true_delay{r}(:, n), '-', 'color', mcol(n,:)); end
end
set(gca, 'xscale', 'log'); xlim([1/fb(end) 1/fb(1)]); xlabel('Period (s)'); ylabel('delay (s)');
title('Path delay of the data relative to Uref (colours: branches)');
subplot(1, 2, 2); hold on; box on;
T_ = cat3('tau'); [~, i40] = min(abs(1./fb - 40)); hh = zeros(1, nm);
for n = 1:nm
    td = cellfun(@(x) x(i40, n), pd.true_delay)';
    he = errorbar(td, mean(T_(:, n, :), 3), 2*std(T_(:, n, :), 0, 3), 'o');
    set(he, 'color', mcol(n,:), 'markerfacecolor', mcol(n,:), 'markersize', 3); hh(n) = he;
end
lim = [-40 40]; plot(lim, lim, 'k--'); xlim(lim); ylim(lim); axis square;
legend(hh, arrayfun(@(n) sprintf('branch %d', n), 0:nm-1, 'uniformoutput', false), 'location', 'northwest', 'fontsize', 7);
xlabel('path delay at 40 s (s)'); ylabel('posterior \tau_n (s)');
title('Time shifts (Love + Rayleigh)');
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig5_path_nuisance.png'));

%% 5b) nuisance parameters: true vs recovered (synthetic test only;
% plot_nuisance_truth.m): tau_n and dth_n Legendre coefficients and B_n(w);
% the true arrival angles come from the A0 truth file (opt.synth_arrival_shift)
if any(isfinite(CtrueL(:)))          % (real data: no true model)
    try
        dth_tr_ = [];
        tf_ = '';
        if isfield(opt, 'aswms_dir') && isfield(opt, 'truthfile'), tf_ = fullfile(opt.aswms_dir, opt.truthfile); end
        if exist(tf_, 'file') == 2
            T_ = load(tf_);
            if isfield(T_.truth, 'dth_true')
                dth_tr_ = zeros(1, nrec);
                for r = 1:nrec
                    e_ = find(strcmp(T_.truth.ev_id, rec(r).id), 1);
                    if ~isempty(e_), dth_tr_(r) = T_.truth.dth_true(e_); end
                end
            end
        end
        Ctrue_f_ = {disp_to_band(CtrueL, Tdisp, fb), disp_to_band(CtrueR, Tdisp, fb)};
        ntxt_ = plot_nuisance_truth(rec, cfg, mc, pd, svc, stc, Ctrue_f_, sig_true, rho_emp, nm, dth_tr_, ...
                                    outdir, ['nuis_true' tag '.mat']);
        fprintf('%s', ntxt_);
        fid_ = fopen(fullfile(outdir, 'summary.txt'), 'a'); fprintf(fid_, '%s', ntxt_); fclose(fid_);
        clear T_ tf_ dth_tr_ Ctrue_f_ ntxt_ fid_ e_
    catch err_
        fprintf('nuisance truth figure not made: %s\n', err_.message);
    end
end

%% 5c) posterior B_n(w) and arrival angles
% posterior spectral corrections B_n(w) (fig13_B_<Love|Rayleigh>.png,
% B_posterior<tag>.mat) and arrival-angle anomalies dth_n(w)
% (fig14_arrival_angles.png, dth_posterior<tag>.mat); synthetic tests: with
% the nuisance fit at the true structure (nuis_true<tag>.mat, section 5b)
% and the true arrival angles of A0
try
    linc_ = []; dtt_ = [];
    if exist(['nuis_true' tag '.mat'], 'file') == 2, L_ = load(['nuis_true' tag '.mat']); linc_ = L_.linc; end
    if isfield(opt, 'aswms_dir') && isfield(opt, 'truthfile') && exist(fullfile(opt.aswms_dir, opt.truthfile), 'file') == 2
        T_ = load(fullfile(opt.aswms_dir, opt.truthfile));
        if isfield(T_.truth, 'dth_true')
            dtt_ = NaN(1, nrec);
            for r = 1:nrec
                e_ = find(strcmp(T_.truth.ev_id, rec(r).id), 1);
                if ~isempty(e_), dtt_(r) = T_.truth.dth_true(e_); end
            end
        end
    end
    ev_ = unique([1 min(8, max([rec.event]))]);
    figs_ = plot_posterior_B(rec, cfg, mc, svc, Tdisp, fb, nm, ev_, outdir, tag, linc_);
    figs_(end+1) = plot_posterior_dth(rec, cfg, svc, nm, ev_, outdir, tag, linc_, dtt_);
    if fig_close_mode, close(figs_); end
    clear figs_ linc_ dtt_ L_ T_ e_ ev_
catch err_
    fprintf('B_n / arrival-angle figures not made: %s\n', err_.message);
end

%% 6) convergence and noise diagnostics
fig = figure('position', [50 50 1900 600]);
subplot(1, 4, 1); hold on; box on;
for ic = 1:nc, plot(svc{ic}.L - Ltrue, 'color', ccol(ic,:)); end
Lall = cell2mat(cellfun(@(x) x.L(round(mc.burn/2):end), svc, 'uniformoutput', false)) - Ltrue;
ylim([min(Lall) - 20, max(max(Lall), 0) + 20]);
plot([1 mc.nit], [0 0], 'r--'); yl = ylim; plot(mc.burn*[1 1], yl, 'k:');
xlabel('iteration'); ylabel('log L - log L(true)'); title('Log likelihood (red: true model)');
wl4 = {'L ovt', 'L fund', 'R ovt', 'R fund'};
if sw_, wl4 = {'L', 'L (tied)', 'R', 'R (tied)'}; end
subplot(1, 4, 2); hold on; box on;
for i = 1:4, plot(i + 0.08*randn(size(S4,1),1), S4(:,i)./sig_true(i), '.', 'markersize', 3); end
plot([0.5 4.5], [1 1], 'k--'); xlim([0.5 4.5]); set(gca, 'xtick', 1:4, 'xticklabel', wl4);
ylabel('\sigma / \sigma_{true}'); title('Noise level');
subplot(1, 4, 3); hold on; box on;
for i = 1:4
    plot(i + 0.08*randn(size(R4,1),1), R4(:,i), '.', 'markersize', 3);
    plot(i + [-0.3 0.3], rho_emp(i)*[1 1], 'k-', 'linewidth', 2);
end
xlim([0.5 4.5]); ylim([0 1]); set(gca, 'xtick', 1:4, 'xticklabel', wl4);
ylabel('lag-1 correlation r'); title('Noise correlation (black: empirical)');
subplot(1, 4, 4); hold on; box on;
for q = 1:3, plot(Rz.(names{q})(in), zplot(in), '-', 'linewidth', 1.3); end
plot([1.1 1.1], [0 400], 'k--'); set(gca, 'ydir', 'reverse'); ylim([0 400]);
xlabel('R-hat'); ylabel('Depth (km)'); title('Gelman-Rubin R-hat');
legend('V_{SV}', 'V_{SH}', '\xi', 'location', 'southeast');
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig6_diagnostics.png'));

%% 8) individual posterior samples (each is a smooth B-spline profile)
fig = figure('position', [50 50 1500 650]);
pickS = round(linspace(1, ns, min(40, ns)));
for q = 1:3
    subplot(1, 3, q); hold on; box on;
    X = V.(names{q});
    plot(X(:, pickS), zplot, '-', 'color', [0.6 0.75 0.95]);
    plot(truv{q}, zplot, 'r-', 'linewidth', 2);
    plot(mean(X, 2), zplot, 'b-', 'linewidth', 1.5);
    set(gca, 'ydir', 'reverse'); ylim([0 400]); xlim(xr{q}); xlabel(labs{q});
    if q == 1, ylabel('Depth (km)'); end
    title(sprintf('%d posterior samples (blue: mean, red: PA5)', numel(pickS)));
end
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig8_samples.png'));

%% 9) eta posterior (only when eta is a model parameter)
if eta_
    fig = figure('position', [50 50 600 560]); hold on; box on;
    X = V.eta; vb = 0.78:0.005:1.07; dv = vb(2) - vb(1);
    H = zeros(numel(zplot), numel(vb));
    for iz = 1:numel(zplot)
        ib = min(max(round((X(iz,:) - vb(1))/dv) + 1, 1), numel(vb));
        H(iz,:) = accumarray(ib(:), 1, [numel(vb) 1])' / ns;
    end
    imagesc(vb, zplot, H); colormap(flipud(gray)); caxis([0 0.25]); set(gca, 'ydir', 'reverse');
    h1 = plot(tru.eta, zplot, 'r-', 'linewidth', 2);
    h4 = plot(refp.eta, zplot, '-', 'color', [0 0.6 0], 'linewidth', 1.2);
    h2 = plot(mean(X, 2), zplot, 'b-', 'linewidth', 1.5);
    h3 = plot(prctile_(X, 5), zplot, 'b--', prctile_(X, 95), zplot, 'b--');
    xlim(vb([1 end])); ylim([0 400]); xlabel('\eta'); ylabel('Depth (km)'); title('Posterior \eta');
    legend([h1 h4 h2 h3(1)], {'true (PA5)', 'reference', 'posterior mean', '5-95%'}, 'location', 'southwest');
    print(fig, '-dpng', '-r130', fullfile(outdir, 'fig9_eta.png'));
end

%% 7) parallel-tempering diagnostics
if ~isempty(infc{1}) && numel(infc{1}.T) > 1
    fig = figure('position', [50 50 1600 560]);
    subplot(1, 3, 1); hold on; box on;
    for ic = 1:nc
        T = infc{ic}.T; sw = infc{ic}.swap(:,2) ./ max(infc{ic}.swap(:,1), 1);
        plot(sqrt(T(1:end-1).*T(2:end)), 100*sw, '-o', 'color', ccol(ic,:), 'markerfacecolor', ccol(ic,:));
    end
    set(gca, 'xscale', 'log'); ylim([0 100]); xlabel('temperature (pair midpoint)');
    ylabel('swap acceptance (%)'); title('Swap acceptance');
    subplot(1, 3, 2); hold on; box on;
    for ic = 1:nc, plot(infc{ic}.T, infc{ic}.Lmean - Ltrue, '-o', 'color', ccol(ic,:)); end
    set(gca, 'xscale', 'log'); xlabel('temperature'); ylabel('mean log L - log L(true)');
    title('Mean log L per temperature');
    subplot(1, 3, 3); hold on; box on;
    K = numel(infc{1}.T); tc = flipud(gray(K + 2));
    for k = K:-1:1, plot(infc{1}.Ltrace(:, k) - Ltrue, 'color', tc(k+1,:)); end
    plot(infc{1}.Ltrace(:, 1) - Ltrue, 'b'); ylim([-300 30]);
    xlabel('iteration'); ylabel('log L - log L(true)'); title('Ensemble 1 by temperature (blue: T = 1)');
    print(fig, '-dpng', '-r130', fullfile(outdir, 'fig7_pt_diagnostics.png'));
end

%% noise diagnostics of the pre-P noise ('aswms' data; plot_noise_diagnostics.m)
try
    plot_noise_diagnostics(rec, cfg, ev, sta, outdir);
catch err_
    fprintf('noise diagnostics not plotted: %s\n', err_.message);
end
