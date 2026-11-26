% D5_STEP2_CHECK_LINEARIZATION  Step 2, part 4 (driver): how accurate is the
% linearized MINEOS dispersion that the waveform inversion uses?
%
% For a finished run (opt.tag; setup<tag>.mat, chain<tag>_XX.mat) and the
% reference model its post-burn-in samples used (the step-1 model, or the
% end-of-burn-in update refupd<tag>.mat, or opt.ref_prep), the phase
% velocities of branches 0..nm-1 (Love and Rayleigh, dispersion nodes
% Tdisp) are computed
%   (a) for the POSTERIOR MEAN model, by full MINEOS (a new card) and by
%       the linearized dispersion of the inversion;
%   (b) for the TRUE model projected on the inversion's parameterization
%       (PA5 Vsv and Vsh between the Moho and 350 km, eta too if it is
%       a model parameter; Vp and density scaled with the reference ratios;
%       the reference model elsewhere), by full MINEOS and linearized.
% Reported (rms over the periods where the branch exists, km/s):
%   E_lin(mean)  = linearized - MINEOS at the posterior mean: the
%                  linearization error at the solution
%   E_lin(true)  = linearized - MINEOS at the projected truth: the
%                  linearization error at the truth
%   E_param      = MINEOS(projected truth) - MINEOS(PA5): the part of the
%                  true dispersion the parameterization cannot represent
%                  (PA5 below 350 km, its Vp/Vs and rho/Vs ratios, eta)
%   E_mean       = MINEOS(posterior mean) - MINEOS(PA5)
% and E_param is split by adding the true model step by step to the
% projected truth (cumulative, MINEOS - PA5 after each step):
%   + eta: PA5 eta between the Moho and 350 km
%   + Vp, rho: PA5 Vpv, Vph and density there (instead of the scaling)
%   + below 350 km: PA5 below 350 km (what remains is the crust/water)
% Updating the reference (end of burn-in, D6) can reduce E_lin, not E_param.
% MINEOS labels branches by their order at each frequency; near an
% osculation (e.g. Love n = 0 and 1 in an oceanic lid) the labels of the
% linearized and the MINEOS dispersion can differ.
%
% Output in opt.outdir: d5_linearization.txt, fig10_linearization.png,
% d5<tag>.mat. Run after D3/D4. Takes ~5 min (10 MINEOS runs).
% Previous: D4_step2_plot_results.m. Next (optional):
% D6_step2_update_reference_and_rerun.m

parameters_RL;
tag = opt.tag; opt_ = opt;
load(['setup' tag '.mat']);
optS = opt; opt = opt_; clear opt_
if ~exist(outdir, 'dir'), mkdir(outdir); end
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
% the reference model the post-burn-in samples used, re-installed from its
% prep file (also gives older setup files the per-parameter kernels)
if exist(['refupd' tag '.mat'], 'file') == 2          % updated at the end of burn-in
    U = load(['refupd' tag '.mat']); U = U.U;
    Pr = load(U.prepfile); Pr = Pr.P;
    refname = sprintf('updated reference %s (end of burn-in)', U.refid);
    if isfield(U, 'reverted') && U.reverted, refname = 'reference (end-of-burn-in update rejected)'; end
elseif isfield(optS, 'ref_prep') && ~isempty(optS.ref_prep)   % D6 run
    Pr = load(optS.ref_prep); Pr = Pr.P;
    refname = sprintf('reference %s', optS.ref_prep);
else
    Pr = prep_mineos_idagrn(optS);
    refname = sprintf('reference %s', optS.refid);
end
[rec, cfg] = apply_reference(rec, cfg, Pr, Tdisp, fb);
card_true = Pr.card_true; clear Pr
eta_ = isfield(mc, 'eta_free') && mc.eta_free;

%% posterior mean
V = 0; Vh = 0; Ve = 0; ns = 0;
for ic = 1:mc.nchains
    D = load(sprintf('chain%s_%02d.mat', tag, ic));
    V = V + sum(D.sv.vsv, 2); Vh = Vh + sum(D.sv.vsh, 2); ns = ns + size(D.sv.vsv, 2);
    if eta_, Ve = Ve + sum(D.sv.eta, 2); end
end
pm.vsv = V/ns; pm.vsh = Vh/ns; pm.eta = []; if eta_, pm.eta = Ve/ns; end
tr.vsv = tru.vsv; tr.vsh = tru.vsh; tr.eta = []; if eta_, tr.eta = tru.eta; end

%% linearized and full-MINEOS dispersion
models = {pm, tr}; mname = {'posterior mean', 'projected truth'};
cardid = {['RLinv_d5_mean' tag], ['RLinv_d5_trueproj' tag]};
Clin = cell(2, 2); Cmin = cell(2, 2);           % {model, wave}
zn = cfg.lin.z;
for k = 1:2
    M = models{k};
    % profiles at the linearization nodes, interpolated within the model
    % domain only (not across the Moho; s1_vsv_at_nodes.m)
    vsv_n = s1_vsv_at_nodes(zplot, M.vsv, zn, cfg.zmin, cfg.zmax);
    vsh_n = s1_vsv_at_nodes(zplot, M.vsh, zn, cfg.zmin, cfg.zmax);
    eta_n = []; if ~isempty(M.eta), eta_n = s1_vsv_at_nodes(zplot, M.eta, zn, cfg.zmin, cfg.zmax); end
    [Clin{k,1}, Clin{k,2}] = mineos_lin_disp(cfg.lin, vsv_n, vsh_n, eta_n, 'both');
    card = ref_card_update(cfg.card_ref, cfg.zw, cfg.zmin, cfg.zmax, zplot, M.vsv, M.vsh, M.eta);
    write_MINEOS_mod(card, fullfile(opt.mdir, 'CARDS', [cardid{k} '.card']));
    copyfile(fullfile(opt.mdir, 'CARDS', [optget(opt, 'ref_qmod', opt.trueid) '.qmod']), fullfile(opt.mdir, 'CARDS', [cardid{k} '.qmod']));
    mineos_lock(opt.mdir, 'acquire');
    try
        for wv = 1:2
            TY = 'TS';
            o = mineos_modes(opt.mdir, cardid{k}, TY(wv), nm, opt.maxF, Tdisp, []);
            Cmin{k, wv} = o.phv(:, 1:nm);
        end
    catch err
        mineos_lock(opt.mdir, 'release'); rethrow(err);
    end
    mineos_lock(opt.mdir, 'release');
    fprintf('%s: MINEOS done\n', mname{k});
end
Ctr = {CtrueL, CtrueR};

%% split of E_param: add the true model step by step (cumulative)
zsf = cfg.card_ref.z - cfg.zw;
tr_eta = card_at_depth(card_true, 'eta', zplot + cfg.zw);
cA = ref_card_update(cfg.card_ref, cfg.zw, cfg.zmin, cfg.zmax, zplot, tru.vsv, tru.vsh, tr_eta);   % + eta
in = find(zsf >= cfg.zmin - 1e-6 & zsf <= cfg.zmax + 1e-6 & cfg.card_ref.vsv > 0);
atm = in(abs(zsf(in) - cfg.zmin) < 1e-6); if numel(atm) > 1, in = setdiff(in, atm(2:end)); end
% true values at reference-card nodes, keeping discontinuities: of a node
% pair at one depth (card order: deeper side first) the first takes the
% value below, the second the value above
zc_ = cfg.card_ref.z(:); nn_ = numel(zc_);
up_ = [false; zc_(2:end) == zc_(1:end-1)];             % upper node of a pair
tv_ = @(f, ii) arrayfun(@(i) card_at_depth(card_true, f, zc_(i), ifelse_side(up_(i))), ii(:));
cB = cA;                                                                   % + true Vp, rho
for f_ = {'vpv', 'vph', 'rho'}, cB.(f_{1})(in) = tv_(f_{1}, in); end
cC = cB; deep = find(zsf > cfg.zmax + 1e-6);                               % + true below 350 km
for f_ = {'vsv', 'vsh', 'vpv', 'vph', 'rho', 'eta'}, cC.(f_{1})(deep) = tv_(f_{1}, deep); end
dcards = {cA, cB, cC}; dname = {'+ eta', '+ Vp, rho', '+ below 350 km'};
Cdec = cell(3, 2);
for k = 1:3
    id_ = sprintf('RLinv_d5_dec%d%s', k, tag);
    write_MINEOS_mod(dcards{k}, fullfile(opt.mdir, 'CARDS', [id_ '.card']));
    copyfile(fullfile(opt.mdir, 'CARDS', [optget(opt, 'ref_qmod', opt.trueid) '.qmod']), fullfile(opt.mdir, 'CARDS', [id_ '.qmod']));
    mineos_lock(opt.mdir, 'acquire');
    try
        for wv = 1:2
            TY = 'TS'; o = mineos_modes(opt.mdir, id_, TY(wv), nm, opt.maxF, Tdisp, []);
            Cdec{k, wv} = o.phv(:, 1:nm);
        end
    catch err
        mineos_lock(opt.mdir, 'release'); rethrow(err);
    end
    mineos_lock(opt.mdir, 'release');
    fprintf('E_param split %s: MINEOS done\n', dname{k});
end

%% errors
rmsn = @(X) sqrt(mean(X.^2, 1, 'omitnan'));
if exist('OCTAVE_VERSION', 'builtin')
    rmsn = @(X) arrayfun(@(j) sqrt(mean(X(isfinite(X(:,j)), j).^2)), 1:size(X, 2));
end
wname = {'Love', 'Rayleigh'};
E = struct();
for wv = 1:2
    nanm = NaN(size(cfg.need)); nanm(cfg.need) = 1;  % 1 where the branch exists, NaN elsewhere
    msk = @(X) X .* nanm;
    E.lin_mean{wv} = rmsn(msk(Clin{1,wv} - Cmin{1,wv}));
    E.lin_true{wv} = rmsn(msk(Clin{2,wv} - Cmin{2,wv}));
    E.param{wv}    = rmsn(msk(Cmin{2,wv} - Ctr{wv}));
    E.mean{wv}     = rmsn(msk(Cmin{1,wv} - Ctr{wv}));
    for k = 1:3, E.dec{k, wv} = rmsn(msk(Cdec{k,wv} - Ctr{wv})); end
end
fid = fopen(fullfile(outdir, 'd5_linearization.txt'), 'w');
for fo = [1 fid]
    etastr = 'at the reference'; if eta_, etastr = 'free'; end
    fprintf(fo, 'Run %s, %s; posterior mean of %d samples; eta %s\n', tag, refname, ns, etastr);
    fprintf(fo, 'rms phase-velocity differences (km/s) per branch n = 0..%d\n', nm - 1);
    for wv = 1:2
        fprintf(fo, '%s\n', wname{wv});
        fprintf(fo, '  E_lin(mean)  linearized - MINEOS, posterior mean   %s\n', mat2str(round(E.lin_mean{wv}*1000)/1000));
        fprintf(fo, '  E_lin(true)  linearized - MINEOS, projected truth   %s\n', mat2str(round(E.lin_true{wv}*1000)/1000));
        fprintf(fo, '  E_param      MINEOS: projected truth - PA5          %s\n', mat2str(round(E.param{wv}*1000)/1000));
        fprintf(fo, '  E_mean       MINEOS: posterior mean - PA5           %s\n', mat2str(round(E.mean{wv}*1000)/1000));
        fprintf(fo, '  split of E_param (MINEOS - PA5, cumulative):\n');
        for k = 1:3
            fprintf(fo, '    %-16s                                 %s\n', dname{k}, mat2str(round(E.dec{k,wv}*1000)/1000));
        end
    end
end
fclose(fid);
save(fullfile(outdir, ['d5' tag '.mat']), 'Clin', 'Cmin', 'Ctr', 'Cdec', 'E', 'pm', 'tr', 'Tdisp', 'refname');

%% figure: differences vs period per branch
fig = figure('position', [50 50 1400 650]);
col = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.5 0.5 0.5];
for wv = 1:2
    for n = 1:nm
        subplot(2, nm, (wv-1)*nm + n); hold on; box on;
        ip = find(cfg.need(:, n)); T = Tdisp(ip);
        h = zeros(1, 4);
        h(1) = plot(T, Clin{1,wv}(ip,n) - Cmin{1,wv}(ip,n), '-o', 'color', col(1,:), 'markersize', 3);
        h(2) = plot(T, Clin{2,wv}(ip,n) - Cmin{2,wv}(ip,n), '-s', 'color', col(2,:), 'markersize', 3);
        h(3) = plot(T, Cmin{2,wv}(ip,n) - Ctr{wv}(ip,n), '-^', 'color', col(3,:), 'markersize', 3);
        h(4) = plot(T, Cmin{1,wv}(ip,n) - Ctr{wv}(ip,n), '--', 'color', col(4,:));
        plot(Tdisp([1 end]), [0 0], 'k:');
        set(gca, 'xscale', 'log'); xlim(Tdisp([1 end]));
        title(sprintf('%s n = %d', wname{wv}, n - 1));
        if n == 1, ylabel('\Deltac (km/s)'); end
        if wv == 2, xlabel('Period (s)'); end
        if wv == 1 && n == 1
            legend(h, {'lin - MINEOS (post. mean)', 'lin - MINEOS (proj. truth)', ...
                'MINEOS proj. truth - PA5', 'MINEOS post. mean - PA5'}, 'location', 'best', 'fontsize', 6);
        end
    end
end
print(fig, '-dpng', '-r130', fullfile(outdir, 'fig10_linearization.png'));
