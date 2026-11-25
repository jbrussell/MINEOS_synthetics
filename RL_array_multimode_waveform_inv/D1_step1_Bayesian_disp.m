% D1_STEP1_BAYESIAN_DISP  Step 1 of the two-step inversion, Bayesian
% dispersion-only version (alternative to D1_step1_rayleigh_isotropic_inv.m).
%
% Transdimensional Bayesian (rjMcMC) inversion of fundamental-mode Rayleigh
% wave phase velocities for Vsv ONLY, with the dispersion-only framework of
% D7 (functions/run_disp_inv.m, lr_common.m): B-spline Vsv profile between
% the Moho (opt.s1b_zmin km below the seafloor) and opt.s1b_zmax, bottom
% tied to the starting model, uniform priors; the linearized MINEOS
% dispersion about the starting model, re-linearized opt.s1b_lin_updates
% times about the posterior mean. Like the linearized D1, the starting
% model's xi, Vpv/Vsv, Vph/Vsh, rho/Vsv and eta are kept exactly
% (cfg.xi_fixed; s1_card_from_vsv.m), so the result has the same form.
% Faster than the waveform version; for users who have measured 1-D
% Rayleigh dispersion (e.g. with ASWMS-Q).
%
% Data (opt.s1b_disp_file):
%   ''   synthetic: MINEOS branch-0 phase velocities of the true model (PA5)
%        at opt.s1_periods, error opt.s1_cstd (fraction of c), plus
%        Gaussian noise of opt.s1b_noise_scale (default 0.5) times that
%        std (opt.s1b_noise, default true; seed
%        opt.s1b_noise_seed); without noise the posterior is centred on an
%        exact fit
%   file .mat with periods (s), c and cstd (km/s), or a text file with
%        three columns: period (s), c, cstd (km/s)
% Starting model: as in D1 (opt.s1_start_model, 'pa5_iso' or 'prem').
%
% Output (same as D1, so D2 can follow either): CARDS/<opt.s1_final>.card
% from the posterior mean Vsv, opt.s1_file (step1_result.mat) with the D1
% fields plus
%   method   'bayes_disp'
%   s1_unc   z (km below the seafloor), vsv_mean, vsv_std, vsv_p05,
%            vsv_p50, vsv_p95 (km/s) on zplot, the posterior samples
%            (vsv_samples, thinned), rhat (per depth), the stage history
% and figs_step1_bayes_disp/ (profile with the posterior band, dispersion
% fit). The posterior std (s1_unc.vsv_std) can be used for the step-2 Vsv
% prior (opt.vsv_ref_sig = 'step1').
%
% Noise level: the data errors are multiplied by a sampled factor
% (opt.s1b_noise_hier, default true; log-uniform hyperprior on
% opt.s1b_noise_range, default [0.1 10]), so the data decide how well they
% can be fit; its posterior is saved (s1_unc.noise_mult) and the scaled
% error bars are drawn next to the assumed ones. Smoothing: hierarchical by
% default (opt.s1b_smooth_hier = true; strength sampled on
% opt.smooth_sig_range). No prior toward a reference Vsv.
%
% Options (opt.s1b_*): disp_file (''), noise (true), noise_seed (1),
% noise_hier (true), noise_range ([0.1 10]),
% zmin (7), zmax (350; km below the seafloor), vsmin/vsmax (3.8/5.0 km/s),
% kvmax (15), smooth_hier (true; false = smooth_sig, default 0 = off, a
% fixed curvature strength on Vsv, lr_common.m),
% chains (4), nit (20000), burn (10000), lin_updates (2).
%
% Next: D2_step2_prep_mineos_idagrn.m

close all;
parameters_RL;
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
mdir = opt.mdir; here = pwd;
cleanupObj = onCleanup(@() cd(here));
od = 'figs_step1_bayes_disp'; if ~exist(od, 'dir'), mkdir(od); end
g = @(f, d) optget(opt, f, d);

%% models (as in D1)
card_prem = make_ref_card(mdir, opt.premref, opt.premid, opt.trueid, opt.zmoho);
card_true = read_model_card(fullfile(mdir, 'CARDS', [opt.trueid '.card']));
if strcmp(opt.s1_start_model, 'pa5_iso')
    card_start = pa5_iso_ref_card(card_true, opt.zmoho, opt.s1_zlin);
else
    card_start = isotropic_card(card_prem);
end
write_MINEOS_mod(card_start, fullfile(mdir, 'CARDS', [opt.s1_start '.card']));
copyfile(fullfile(mdir, 'CARDS', [opt.ref_qmod '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_start '.qmod']));
card_start = read_model_card(fullfile(mdir, 'CARDS', [opt.s1_start '.card']));
zw = max(card_start.z(card_start.vsv == 0 & card_start.z < 20));     % water depth (km)

%% data
fdat = g('s1b_disp_file', '');
if isempty(fdat)
    periods = opt.s1_periods(:);
    cd(mdir); cobs = run_mineos_phv_kernels(opt.trueid, 'R', periods, 0, 0); cd(here);
    cstd = opt.s1_cstd * cobs;
    if g('s1b_noise', true)
        rs_ = randn('state'); randn('state', g('s1b_noise_seed', 1));
        cobs = cobs + optget(opt, 's1b_noise_scale', 0.5) * cstd .* randn(size(cobs)); randn('state', rs_);
    end
    dlabel = sprintf('synthetic (PA5), %.2f%% error', 100*opt.s1_cstd);
    if g('s1b_noise', true), dlabel = [dlabel sprintf(' + noise %.2fx (seed %d)', g('s1b_noise_scale', 0.5), g('s1b_noise_seed', 1))]; end
elseif any(regexp(fdat, '\.mat$'))
    D_ = load(fdat); periods = D_.periods(:); cobs = D_.c(:); cstd = D_.cstd(:); dlabel = fdat;
else
    D_ = load(fdat); periods = D_(:, 1); cobs = D_(:, 2); cstd = D_(:, 3); dlabel = fdat;
end
[periods, io_] = sort(periods); cobs = cobs(io_); cstd = cstd(io_);

%% model space and sampler settings (Vsv only)
zmin = g('s1b_zmin', 7); zmax = g('s1b_zmax', 350); zplot = (0:2:400)';
cfg = struct('zmin', zmin, 'zmax', zmax, 'moho', zmin, 'zw', zw, 'xi_fixed', true);
cp_ = @(f) card_at_depth(card_start, f, zplot + zw)/1000;
cfg.refp = struct('z', zplot, 'vsv', cp_('vsv'), 'vsh', cp_('vsh'), 'eta', card_at_depth(card_start, 'eta', zplot + zw));
mc = struct('kvmin', 0, 'kvmax', g('s1b_kvmax', 15), 'kxmin', 0, 'kxmax', 0, ...
    'vsmin', g('s1b_vsmin', 3.8), 'vsmax', g('s1b_vsmax', 5.0), 'ximin', 0.85, 'ximax', 1.15, 'p_xi', 0, ...
    'step_vs', 0.03, 'step_xi', 0.01, 'step_z', 10, 'eta_free', false, ...
    'smooth_sig', g('s1b_smooth_sig', 0)*[1 1], 'smooth_param', 'vsv_vsh', 'smooth_len', 10, ...
    'smooth_hier', g('s1b_smooth_hier', true), 'smooth_sig_range', g('smooth_sig_range', [5e-5 5e-3]), ...
    'disp_noise_hier', g('s1b_noise_hier', true), 'disp_noise_range', g('s1b_noise_range', [0.1 10]));
nch = g('s1b_chains', 4); nit = g('s1b_nit', 20000); burn = g('s1b_burn', 10000); thin = 20;
nup = g('s1b_lin_updates', 2);
sig = {NaN(numel(periods), 1), cstd ./ cobs};      % relative errors; Love not used (w = [0 1])

%% sampling, re-linearized about the posterior mean
card_lin = card_start; hist_ = struct('stage', {}, 'dvsv_rms', {}, 'chi2', {});
for stg = 0:nup
    [lin, cR0] = s1_lin_at(mdir, card_lin, opt.s1_work, opt.ref_qmod, periods, opt.maxF, zw, zmin, zmax);
    cfg.lin = lin; cfg.need = isfinite(cR0);
    if stg == 0                                     % bottom tied to the starting model at zmax
        [~, ib_] = max(lin.z); cfg.vbot = lin.vsv0(ib_); cfg.xbot = (lin.vsh0(ib_)/lin.vsv0(ib_))^2;
        cstart = cR0; idx = lin.idx; zn = lin.z;
    end
    R = cell(nch, 1); tic;
    for k = 1:nch
        R{k} = run_disp_inv(NaN(numel(periods), 1), cobs, sig, cfg, mc, periods, 1./periods, 1, zplot, ...
                            nit, burn, thin, 500 + 10*stg + k, [0 1]);
    end
    V = cell2mat(cellfun(@(x) x.vsv, R', 'UniformOutput', false));
    vmean = mean(V, 2);
    vnew = interp1(zplot, vmean, zn, 'linear', 'extrap');
    dv = vnew - lin.vsv0;
    CRs = cell2mat(cellfun(@(x) squeeze(x.CR), R', 'UniformOutput', false));
    chi2 = mean(((mean(CRs, 2) - cobs)./cstd).^2);
    hist_(end+1) = struct('stage', stg, 'dvsv_rms', sqrt(mean(dv.^2)), 'chi2', chi2); %#ok<SAGROW>
    fprintf('stage %d (%.0f s): posterior mean - linearization model, rms %.3f km/s; reduced chi^2 (linearized) %.2f; acceptance %s\n', ...
        stg, toc, hist_(end).dvsv_rms, chi2, mat2str(round(100*cellfun(@(x) x.acc, R'))/100));
    card_lin = s1_card_from_vsv(card_start, idx, vnew);
end

%% result: posterior statistics and the step-1 card
in = zplot >= zmin & zplot <= zmax;
q = sort(V, 2); ns = size(q, 2); pc = @(p) q(:, max(1, round(p*ns)));
rh = NaN(numel(zplot), 1);
for i = find(in)'
    X = cell2mat(cellfun(@(x) x.vsv(i, :)', R', 'UniformOutput', false)); rh(i) = gr_rhat(X);
end
qc = sort(CRs, 2); nsc = size(qc, 2); pcc = @(p) qc(:, max(1, round(p*nsc)));   % predicted dispersion (last stage, linearized)
lam = sort(cell2mat(cellfun(@(x) x.lam(:, 2), R, 'UniformOutput', false)));   % Rayleigh noise multiplier (1 if not sampled)
nml = struct('median', median(lam), 'p05', lam(max(1, round(0.05*numel(lam)))), 'p95', lam(round(0.95*numel(lam))), ...
    'sampled', mc.disp_noise_hier);
kvs = cell2mat(cellfun(@(x) x.kv(:), R, 'UniformOutput', false));   % interior Vsv control points
kvst = struct('median', median(kvs), 'p05', prctile_(kvs', 5), 'p95', prctile_(kvs', 95));
ssv = sort(cell2mat(cellfun(@(x) x.ss(:, 1), R, 'UniformOutput', false)));   % sampled Vsv curvature strength (NaN if fixed)
ssv = ssv(isfinite(ssv)); sms = struct('sampled', mc.smooth_hier, 'median', NaN, 'p05', NaN, 'p95', NaN);
if ~isempty(ssv), sms.median = median(ssv); sms.p05 = ssv(max(1, round(0.05*numel(ssv)))); sms.p95 = ssv(round(0.95*numel(ssv))); end
s1_unc = struct('z', zplot, 'vsv_mean', vmean, 'vsv_std', std(V, 0, 2), 'vsv_p05', pc(0.05), 'vsv_p16', pc(0.16), ...
    'vsv_p50', pc(0.5), 'vsv_p84', pc(0.84), 'vsv_p95', pc(0.95), 'vsv_samples', V(:, 1:max(1, floor(ns/400)):end), ...
    'c_p05', pcc(0.05), 'c_p16', pcc(0.16), 'c_p84', pcc(0.84), 'c_p95', pcc(0.95), ...
    'rhat', rh, 'stages', hist_, 'zmin', zmin, 'zmax', zmax, 'data', dlabel, 'noise_mult', nml, 'smooth_sig', sms, 'kv', kvst);
card_final = s1_card_from_vsv(card_start, idx, interp1(zplot, vmean, zn, 'linear', 'extrap'));
write_MINEOS_mod(card_final, fullfile(mdir, 'CARDS', [opt.s1_final '.card']));
copyfile(fullfile(mdir, 'CARDS', [opt.ref_qmod '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_final '.qmod']));
cd(mdir); cpre = run_mineos_phv_kernels(opt.s1_final, 'R', periods, 0, 0); cd(here);
discs = []; method = 'bayes_disp';
fprintf('Rayleigh fundamental mode, reduced chi^2 (MINEOS): start %.1f, posterior-mean model %.2f\n', ...
    mean(((cobs - cstart)./cstd).^2), mean(((cobs - cpre)./cstd).^2));
if mc.disp_noise_hier
    fprintf('data error multiplier (sampled): median %.2f [5-95%%: %.2f-%.2f]; reduced chi^2 with the scaled errors %.2f\n', ...
        nml.median, nml.p05, nml.p95, mean(((cobs - cpre)./(nml.median*cstd)).^2));
end
if mc.smooth_hier
    fprintf('smoothing strength (sampled): median %.2e [5-95%%: %.2e-%.2e] km/s/km^2 (hyperprior range %s)\n', ...
        sms.median, sms.p05, sms.p95, mat2str(mc.smooth_sig_range));
end
fprintf('interior Vsv control points: median %g [5-95%%: %g-%g]\n', kvst.median, kvst.p05, kvst.p95);
fprintf('Vsv posterior: R-hat median %.2f (max %.2f); std median %.3f km/s\n', median(rh(in & isfinite(rh))), max(rh(in & isfinite(rh))), median(s1_unc.vsv_std(in)));
save(opt.s1_file, 'card_prem', 'card_start', 'card_final', 'card_true', 'cobs', 'cstd', ...
     'cstart', 'cpre', 'periods', 'discs', 'opt', 'method', 's1_unc');

%% figures
cpt = @(c, f) card_at_depth(c, f, zplot + zw)/1000;
fig = figure('position', [50 50 2100 650]);
subplot(1, 3, 1); hold on; box on;
fill([s1_unc.vsv_p05(in); flipud(s1_unc.vsv_p95(in))], [zplot(in); flipud(zplot(in))], [0.80 0.88 1], 'edgecolor', 'none');
fill([s1_unc.vsv_p16(in); flipud(s1_unc.vsv_p84(in))], [zplot(in); flipud(zplot(in))], [0.55 0.70 0.95], 'edgecolor', 'none');
h1 = plot(cpt(card_true, 'vsv'), zplot, 'r-', 'linewidth', 2);
h2 = plot(cpt(card_start, 'vsv'), zplot, 'g-', 'linewidth', 1.5);
h3 = plot(vmean, zplot, 'b-', 'linewidth', 2);
set(gca, 'ydir', 'reverse'); ylim([0 400]); xlim([3.8 5]);
xlabel('V_{SV} (km/s)'); ylabel('Depth below seafloor (km)');
legend([h1 h2 h3], {'true (PA5)', 'start', 'posterior mean'}, 'location', 'southwest');
title('Step 1, Bayesian (dispersion): V_{SV} (shaded: 5-95%, 16-84%)');
subplot(1, 3, 2); hold on; box on;
fill([periods; flipud(periods)], [s1_unc.c_p05; flipud(s1_unc.c_p95)], [0.80 0.88 1], 'edgecolor', 'none');
fill([periods; flipud(periods)], [s1_unc.c_p16; flipud(s1_unc.c_p84)], [0.55 0.70 0.95], 'edgecolor', 'none');
hd = errorbar(periods, cobs, 2*cstd, 'ko');
hn = errorbar(periods*1.015, cobs, 2*nml.median*cstd, 'r.');            % posterior scaled errors (slightly offset)
set(hn, 'linewidth', 1.5);
hs = plot(periods, cstart, 'g-', 'linewidth', 1.5); hp = plot(periods, cpre, 'b-', 'linewidth', 1.5);
xlabel('Period (s)'); ylabel('Rayleigh phase velocity (km/s)');
legend([hd hn hs hp], {'data, assumed 2\sigma', ['posterior scaled 2\sigma (x' sprintf('%.2f', nml.median) ')'], 'start', ...
    'posterior mean model (MINEOS)'}, 'location', 'southeast');
title({sprintf('Data: %s', dlabel), 'shaded: posterior predicted 5-95%, 16-84% (linearized)'}, 'interpreter', 'none');
KVC = cellfun(@(x) x.kv, R, 'UniformOutput', false); KVMAX = mc.kvmax;
subplot(1, 3, 3); hold on; box on;                     % control points per chain (as D4's fig2)
cc_ = lines(max(numel(KVC), 7)); hk_ = zeros(1, numel(KVC)); lk_ = cell(1, numel(KVC));
for ic_ = 1:numel(KVC)
    cnt_ = accumarray(KVC{ic_}(:) + 1, 1, [KVMAX+1 1]) / numel(KVC{ic_});
    hk_(ic_) = plot(0:KVMAX, cnt_, '-o', 'color', cc_(ic_,:), 'markerfacecolor', cc_(ic_,:), 'markersize', 4);
    lk_{ic_} = sprintf('chain %d (median %g)', ic_, median(KVC{ic_}));
end
xlim([-0.5 KVMAX+0.5]); xlabel('interior V_{SV} control points'); ylabel('fraction of samples');
legend(hk_, lk_, 'location', 'northeast'); title('V_{SV} control points per chain');
print(fig, '-dpng', '-r110', fullfile(od, 'step1_bayes_disp.png')); close(fig);
fprintf('step-1 model CARDS/%s.card, %s, figures in %s/\n', opt.s1_final, opt.s1_file, od);
