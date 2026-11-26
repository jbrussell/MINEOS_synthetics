% D1_STEP1_BAYESIAN_WAVEFORM  Step 1 of the two-step inversion, Bayesian
% waveform version (alternative to D1_step1_rayleigh_isotropic_inv.m).
%
% Vsv-only transdimensional Bayesian inversion of the RAYLEIGH waveforms
% (ASWMS-Q event files of A1-A3, vertical component) with the step-2
% waveform framework (D3_step2_waveform_inv.m), set up for step 1:
%   - reference model = the step-1 STARTING model (opt.s1_start_model, as
%     in D1): its idagrn6 branch synthetics (D2's prep, opt.s1w_nm
%     branches, default 1 = fundamental-mode synthetics only) and MINEOS
%     linearization, re-linearized opt.s1w_lin_updates times about the
%     posterior mean during burn-in (opt.update_ref = 'lin')
%   - Rayleigh records only (opt.waves = 2)
%   - Vsv only: xi, Vpv/Vsv, Vph/Vsh, rho/Vsv and eta of the starting model
%     are kept exactly (opt.xi_fixed; s1_card_from_vsv.m), as in D1
%   - nuisance (B, tau, arrival angle per event), noise level and
%     correlation sampled as in step 2
% The posterior mean Vsv gives the step-1 card (the step-2 reference); the
% posterior spread its uncertainty.
%
% Output (same as D1, so D2 can follow): CARDS/<opt.s1_final>.card,
% opt.s1_file with the D1 fields (cobs: MINEOS branch-0 Rayleigh phase
% velocities of the true model at opt.s1_periods, for comparison only;
% cpre: of the step-1 card) plus method = 'bayes_waveform' and s1_unc (z,
% vsv_mean, vsv_std, vsv_p05/p50/p95, vsv_samples, rhat; as in
% D1_step1_Bayesian_disp.m); the D4 figures of the run in
% figs_step1_bayes_wf/ (opt.s1w_tag run) and step1_bayes_wf.png there.
%
% Options (opt.s1w_*): nm (1), sampler ('pt'), nens (4; or nchains for
% 'plain'), nit (1600), burn (800), lin_updates (2), smooth_sig (0 = off;
% fixed curvature strength), smooth_hier (false; true = sampled), no
% prior toward a reference Vsv (vsv_ref_sig = 0),
% tag ('_s1wf'), prepfile ('prep_s1wf.mat'), launch ('serial': all
% ensembles in this session - MATLAB parfor runs them in parallel with a
% pool; 'parallel': one background Octave process per ensemble, Octave
% only). Data: opt.aswms_dir, opt.aswms_components{2} (Rayleigh).
%
% Next: D2_step2_prep_mineos_idagrn.m

parameters_RL;
if fig_close_mode, close all; end
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
mdir = opt.mdir; here = pwd;
g = @(f, d) optget(opt, f, d);
od = 'figs_step1_bayes_wf'; if ~exist(od, 'dir'), mkdir(od); end

%% models (as in D1)
card_true = read_model_card(fullfile(mdir, 'CARDS', [opt.trueid '.card']));
card_start = pa5_iso_ref_card(card_true, opt.zmoho, opt.s1_zlin);
write_MINEOS_mod(card_start, fullfile(mdir, 'CARDS', [opt.s1_start '.card']));
copyfile(fullfile(mdir, 'CARDS', [opt.ref_qmod '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_start '.qmod']));
card_start = read_model_card(fullfile(mdir, 'CARDS', [opt.s1_start '.card']));

%% options of the step-1 waveform run (copied into D3's opt)
o1 = opt;
o1.tag = g('s1w_tag', '_s1wf'); o1.outdir = od;
o1.prepfile = g('s1w_prepfile', 'prep_s1wf.mat');
o1.refid = opt.s1_start; o1.ref_card = card_start;
o1.nm = g('s1w_nm', 1); o1.waves = 2; o1.xi_fixed = true;
o1.sampler = g('s1w_sampler', 'pt'); o1.update_ref = 'lin'; o1.lin_updates = g('s1w_lin_updates', 2);
o1.smooth_sig = g('s1w_smooth_sig', 0)*[1 1]; o1.smooth_hier = g('s1w_smooth_hier', false);   % (true: sampled strength; not recommended)
o1.vsv_ref_sig = 0; o1.xi_ref_sig = 0; o1.xi_ref_hier = false;
if strcmp(o1.sampler, 'pt')
    o1.nens = g('s1w_nens', 4); o1.pt_nit = g('s1w_nit', 1600); o1.pt_burn = g('s1w_burn', 800); nrun = o1.nens;
else
    o1.nchains = g('s1w_nens', 4); o1.nit = g('s1w_nit', 3000); o1.burn = g('s1w_burn', 1500); nrun = o1.nchains;
end

%% reference (starting-model) branch synthetics and linearization (D2's prep)
if exist(o1.prepfile, 'file') ~= 2
    prep_mineos_idagrn(o1);
end

%% the waveform inversion (D3), then its figures (D4)
if strcmp(o1.sampler, 'pt'), envv = 'ENSEMBLES'; else, envv = 'CHAINS'; end
opt0 = opt;
stf_ = ['s1wf_state' o1.tag '.mat'];         % D3/D4 are scripts in this workspace: keep our variables
save(stf_, 'o1', 'nrun', 'card_start', 'card_true', 'mdir', 'here', 'od', 'opt0', 'envv', 'stf_');
opt = o1;
setenv(envv, '[]'); D3_step2_waveform_inv;                        % setup only
if strcmp(g('s1w_launch', 'serial'), 'parallel') && exist('OCTAVE_VERSION', 'builtin')
    optf = ['s1wf_opt' o1.tag '.mat']; pth = path; save(optf, 'o1', 'pth');
    oc = fullfile(OCTAVE_HOME, 'bin', 'octave-cli'); if exist(oc, 'file') ~= 2, oc = 'octave-cli'; end
    for i = 1:nrun
        system(sprintf('cd "%s" && %s=%d "%s" -q --eval "load(''%s''); path(pth); opt = o1; D3_step2_waveform_inv" > run%s.%d.log 2>&1 &', ...
            here, envv, i, oc, optf, o1.tag, i));
    end
    chf = @(i) sprintf('chain%s_%02d.mat', o1.tag, i);
    while ~all(arrayfun(@(i) exist(chf(i), 'file') == 2, 1:nrun)), pause(60); end
else
    setenv(envv, ''); opt = o1; D3_step2_waveform_inv;              % all ensembles here
end
load(stf_); setenv(envv, '');
opt = o1; D4_step2_plot_results;
load(stf_); opt = opt0; g = @(f, d) optget(opt, f, d); delete(stf_);

%% posterior Vsv -> step-1 card and uncertainty
S = load(['setup' o1.tag '.mat'], 'cfg', 'zplot', 'Tdisp');
zplot = S.zplot; cfg = S.cfg; zw = cfg.zw; Td = S.Tdisp(:);
V = []; Vc = cell(nrun, 1); CRs = []; SS = []; KV = []; KVC = cell(nrun, 1);
for i = 1:nrun
    C = load(sprintf('chain%s_%02d.mat', o1.tag, i), 'sv'); Vc{i} = C.sv.vsv; V = [V, C.sv.vsv]; %#ok<AGROW>
    CRs = [CRs, squeeze(C.sv.CR(:, 1, :))]; %#ok<AGROW>   % branch-0 Rayleigh, posterior samples (linearized)
    if isfield(C.sv, 'ss'), SS = [SS; C.sv.ss(:, 1)]; end %#ok<AGROW>   % sampled Vsv curvature strength
    KV = [KV; C.sv.kv(:)]; %#ok<AGROW>                                    % interior Vsv control points
    KVC{i} = C.sv.kv(:);
end
in = zplot >= cfg.zmin & zplot <= cfg.zmax;
vmean = mean(V, 2); q = sort(V, 2); ns = size(q, 2); pc = @(p) q(:, max(1, round(p*ns)));
rh = NaN(numel(zplot), 1); nmin = min(cellfun(@(x) size(x, 2), Vc));
for j = find(in)'
    rh(j) = gr_rhat(cell2mat(cellfun(@(x) x(j, 1:nmin)', Vc', 'UniformOutput', false)));
end
qc = sort(CRs, 2); nsc = size(qc, 2); pcc = @(p) qc(:, max(1, round(p*nsc)));
SS = sort(SS(isfinite(SS))); sms = struct('sampled', o1.smooth_hier, 'median', NaN, 'p05', NaN, 'p95', NaN);
if ~isempty(SS), sms.median = median(SS); sms.p05 = SS(max(1, round(0.05*numel(SS)))); sms.p95 = SS(round(0.95*numel(SS))); end
s1_unc = struct('z', zplot, 'vsv_mean', vmean, 'vsv_std', std(V, 0, 2), 'vsv_p05', pc(0.05), 'vsv_p16', pc(0.16), ...
    'vsv_p50', pc(0.5), 'vsv_p84', pc(0.84), 'vsv_p95', pc(0.95), 'vsv_samples', V(:, 1:max(1, floor(ns/400)):end), ...
    'c_periods', Td, 'c_p05', pcc(0.05), 'c_p16', pcc(0.16), 'c_p84', pcc(0.84), 'c_p95', pcc(0.95), 'rhat', rh, ...
    'smooth_sig', sms, 'kv', struct('median', median(KV), 'p05', prctile_(KV', 5), 'p95', prctile_(KV', 95)), ...
    'zmin', cfg.zmin, 'zmax', cfg.zmax, 'data', sprintf('Rayleigh waveforms, %s', opt.aswms_dir));
% the card: the starting card's nodes between zmin and zmax (the inversion
% nodes of the setup, whose reference card is the starting card) get the
% posterior mean Vsv
card_final = s1_card_from_vsv(card_start, cfg.lin.idx, s1_vsv_at_nodes(zplot, vmean, cfg.lin.z, cfg.zmin, cfg.zmax));
write_MINEOS_mod(card_final, fullfile(mdir, 'CARDS', [opt.s1_final '.card']));
copyfile(fullfile(mdir, 'CARDS', [opt.ref_qmod '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_final '.qmod']));
periods = opt.s1_periods(:);
cd(mdir);
cobs = run_mineos_phv_kernels(opt.trueid, 'R', periods, 0, 0);
cstart = run_mineos_phv_kernels(opt.s1_start, 'R', periods, 0, 0);
cpre = run_mineos_phv_kernels(opt.s1_final, 'R', periods, 0, 0);
cd(here);
cstd = opt.s1_cstd * cobs; discs = []; method = 'bayes_waveform';
fprintf('branch-0 Rayleigh phase velocity vs the true model (MINEOS): rms start %.4f, step-1 card %.4f km/s\n', ...
    sqrt(mean((cstart - cobs).^2)), sqrt(mean((cpre - cobs).^2)));
fprintf('Vsv posterior: R-hat median %.2f (max %.2f); std median %.3f km/s\n', ...
    median(rh(in & isfinite(rh))), max(rh(in & isfinite(rh))), median(s1_unc.vsv_std(in)));
save(opt.s1_file, 'card_start', 'card_final', 'card_true', 'cobs', 'cstd', ...
     'cstart', 'cpre', 'periods', 'discs', 'opt', 'method', 's1_unc');

%% figure (as D1_step1_Bayesian_disp.m)
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
title('Step 1, Bayesian (Rayleigh waveforms): V_{SV} (shaded: 5-95%, 16-84%)');
subplot(1, 3, 2); hold on; box on;
ok_ = isfinite(s1_unc.c_p05);
fill([Td(ok_); flipud(Td(ok_))], [s1_unc.c_p05(ok_); flipud(s1_unc.c_p95(ok_))], [0.80 0.88 1], 'edgecolor', 'none');
fill([Td(ok_); flipud(Td(ok_))], [s1_unc.c_p16(ok_); flipud(s1_unc.c_p84(ok_))], [0.55 0.70 0.95], 'edgecolor', 'none');
hd = plot(periods, cobs, 'ko'); hs = plot(periods, cstart, 'g-', 'linewidth', 1.5); hp = plot(periods, cpre, 'b-', 'linewidth', 1.5);
xlabel('Period (s)'); ylabel('Rayleigh phase velocity, branch 0 (km/s)');
legend([hd hs hp], {'true model', 'start', 'step-1 card'}, 'location', 'southeast');
title({'Branch-0 Rayleigh dispersion (MINEOS)', 'shaded: posterior 5-95%, 16-84% (linearized)'});
KVMAX = optget(o1, 'kvmax', 15);
subplot(1, 3, 3); hold on; box on;                     % control points per chain (as D4's fig2)
cc_ = lines(max(numel(KVC), 7)); hk_ = zeros(1, numel(KVC)); lk_ = cell(1, numel(KVC));
for ic_ = 1:numel(KVC)
    cnt_ = accumarray(KVC{ic_}(:) + 1, 1, [KVMAX+1 1]) / numel(KVC{ic_});
    hk_(ic_) = plot(0:KVMAX, cnt_, '-o', 'color', cc_(ic_,:), 'markerfacecolor', cc_(ic_,:), 'markersize', 4);
    lk_{ic_} = sprintf('chain %d (median %g)', ic_, median(KVC{ic_}));
end
xlim([-0.5 KVMAX+0.5]); xlabel('interior V_{SV} control points'); ylabel('fraction of samples');
legend(hk_, lk_, 'location', 'northeast'); title('V_{SV} control points per chain');
print(fig, '-dpng', '-r110', fullfile(od, 'step1_bayes_wf.png')); if fig_close_mode, close(fig); end
fprintf('step-1 model CARDS/%s.card, %s, figures in %s/\n', opt.s1_final, opt.s1_file, od);
