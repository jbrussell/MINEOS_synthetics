% D3_STEP2_WAVEFORM_INV  Step 2, part 2 (driver): joint Love + Rayleigh
% multimode array waveform inversion (transdimensional Bayesian, rjMcMC)
% for Vsv, Vsh and xi, with the step-1 model as the reference model.
%
% Sampler (opt.sampler in parameters_RL.m):
%   'pt'    : (default) opt.nens parallel-tempering ensembles of opt.K
%             temperatures up to opt.Tmax (swaps between neighbouring
%             temperatures)
%   'plain' : opt.nchains independent rjMcMC chains from random prior draws
% Nuisance parameters, per record and branch: spectral amplitudes B_n(w),
% time shifts tau_n(w) and arrival angles dth_n(w) (opt.nuisance):
%   'marginal' : all integrated out (B exactly; tau, dth by linearization
%                about linearization points that are moved to their
%                posterior mean every opt.relin iterations, default 1 =
%                before every proposal; run_pt_LR.m, nuis_relin.m). Only Vsv,
%                xi and the noise parameters are sampled and tempered.
%   'sample'   : tau and dth sampled with the exact waveforms, B integrated
%                out (run_chain_LR.m; plain sampler only)
% The linearization points start from a Gauss-Newton fit of the nuisance
% parameters at the reference dispersion (spec_refine.m).
%
% Reference update (opt.update_ref):
%   'burnin' (default): the reference model (idagrn6_mask reference
%       spectra, MINEOS phase velocities and kernels) is recomputed ONCE, at
%       the end of burn-in, for the mean Vsv and Vsh of the last
%       opt.upd_frac of the burn-in samples of ALL chains/ensembles (T = 1);
%       every chain then continues from its own state with that common
%       reference, fixed for the rest of the run. With parallel processes,
%       the process that finds all burn-in files first does the update
%       (MINEOS + idagrn6, ~10 min) and the others wait for it.
%   'lin': as 'burnin', but only the MINEOS LINEARIZATION (cfg.lin: phase
%       velocities and kernels used to predict the dispersion of a model) is
%       moved to the mean model; the waveform reference (rec.Uref, rec.cref,
%       cfg.need) stays that of step 1, so the nuisance parameters keep their
%       meaning and every chain keeps its own linearization points. (A full
%       update changes the path delays tau_n(w) must absorb relative to the
%       reference; with a biased mean model these are larger and less smooth
%       than the 3-term tau_n(w) basis can fit.) Always kept; the
%       log-likelihood change of the end-of-burn-in states is only reported.
%       opt.lin_updates = n > 1: the burn-in is split into n equal stages
%       and the linearization is updated after each (lin_update.m), so the
%       kernels can follow the model as it moves away from the step-1
%       reference; the sampling uses the last one.
%   'none': the reference of opt.prepfile (step 1) throughout.
%   The update is kept only if the end-of-burn-in (T = 1) states of the
%   chains fit the data at least as well with it on average (refupd<tag>.mat
%   records the log-likelihood changes dL and whether it was reverted).
% opt.ref_prep: OPTIONAL prep file of another reference model to start from
%   (used by D6_step2_update_reference_and_rerun.m).
% The data, windows and normalizations are those of the step-1 reference
% in all cases, so likelihoods stay comparable.
%
%   >> D3_step2_waveform_inv                       % all ensembles here (parallel tempering)
%   $ ENSEMBLES="[]" octave-cli --eval D3_step2_waveform_inv     % setup only
%   $ ENSEMBLES="1 2" octave-cli --eval D3_step2_waveform_inv &  % ensembles in parallel
%   $ CHAINS="1" octave-cli --eval "opt.sampler='plain'; D3_step2_waveform_inv" &   % plain chains
%
% Parallel (MATLAB, Parallel Computing Toolbox): the chains/ensembles run in
% parfor loops (one per phase/stage; the MINEOS updates between them run in
% the client). opt.parallel = [] (default): MATLAB's own pool rule (an open
% pool, or an automatically created one); true: a pool of opt.par_workers
% (default: one per chain/ensemble) is started (D3_step2_waveform_inv_parallel.m);
% false: serial. In Octave the loops run serially (use separate processes,
% ENSEMBLES/CHAINS, as above).
%
% Checkpoints: every opt.ckpt_every (100) iterations each chain/ensemble
% saves ckpt<tag>_XX_<phase>.mat (run_pt_LR.m); rerunning the same command
% after an interruption resumes from them (and skips finished stages).
%
% Output: setup<tag>.mat, chain<tag>_XX.mat (tag = '_plain' or '_pt' by
% default), with the reference update also chain<tag>_XX_burn.mat (burn-in
% states), refupd<tag>.mat and prep<tag>_burnupd.mat (the new reference);
% when all chains/ensembles exist, the figures (D4).
% Previous: D2_step2_prep_mineos_idagrn.m. Next: D4_step2_plot_results.m

parameters_RL;
if fig_close_mode, close all; end
sub_setup_records;
tag = opt.tag;
setupf = ['setup' tag '.mat'];
pt_ = strcmp(opt.sampler, 'pt');
sample_ = strcmp(opt.nuisance, 'sample');
if pt_ && sample_, error('opt.nuisance = ''sample'' is only available with the plain sampler'); end
upd_ = any(strcmp(opt.update_ref, {'burnin', 'lin'}));
linonly_ = strcmp(opt.update_ref, 'lin');
nlu_ = max(1, round(optget(opt, 'lin_updates', 1)));   % linearization updates during the burn-in ('lin')
if upd_ && sample_
    fprintf('note: no reference update with opt.nuisance = ''sample''\n'); upd_ = false;
end
if ~isempty(opt.ref_prep)                         % start from another reference (D6)
    Pst = load(opt.ref_prep); Pst = Pst.P;
    [rec, cfg, Cref_f, CrefL, CrefR] = apply_reference(rec, cfg, Pst, Tdisp, fb);
    for r = 1:numel(rec)
        pd.true_delay{r} = rec(r).Del0 * (1./Ctrue_f{rec(r).wave} - 1./Cref_f{rec(r).wave});
    end
    fprintf('reference model: %s\n', opt.ref_prep); clear Pst
end
if pt_
    pt.T = logspace(0, log10(opt.Tmax), opt.K); mc.nchains = opt.nens;
    pt.nit = opt.pt_nit; pt.burn = opt.pt_burn;
else
    pt.T = 1; pt.nit = mc.nit; pt.burn = mc.burn;     % one temperature: ordinary chains
end
pt.thin = mc.thin; mc.nit = pt.nit; mc.burn = pt.burn;

if exist(setupf, 'file') == 2
    load(setupf, 'Ltrue', 'lin0');
else
    % initial linearization points: nuisance fit at the reference dispersion
    lin0 = spec_refine([], rec, Cref_f, sig_true, rho_emp, cfg, mc, nm, 30);
    % log-likelihood of the true structure, with nuisance parameters fitted
    % to it (idagrn data have no "true" nuisance parameters)
    if any(isfinite(Ctrue_f{1}(:))) && any(isfinite(Ctrue_f{2}(:)))   % (real data: no true model)
        Ltrue = fit_nuisance_LR(rec, Ctrue_f, sig_true, rho_emp, cfg, mc, nm);
        fprintf('logL of the true structure (fitted nuisance, sigma_true): %.1f\n', Ltrue);
    else
        Ltrue = NaN;
    end
    tmpf = tmp_name('setup');
    save(tmpf, 'tru', 'refp', 'cfg', 'mc', 'Tdisp', 'CtrueL', 'CtrueR', 'CrefL', 'CrefR', 'ClinL', 'ClinR', ...
         'zplot', 'sig_true', 'rho_emp', 'rec', 'sta', 'Ltrue', 'fb', 'nm', 'outdir', 'tag', 'pd', ...
         'opt', 'ev', 'lin0', 'pt');
    movefile(tmpf, setupf);
end

if pt_, envc = getenv('ENSEMBLES'); else, envc = getenv('CHAINS'); end
runs = 1:mc.nchains;
if ~isempty(envc), runs = str2num(envc); end %#ok<ST2NM>
parw_ = par_workers(optget(opt, 'parallel', []), optget(opt, 'par_workers', []), numel(runs));   % (MATLAB pool)
chainf = @(i) sprintf('chain%s_%02d.mat', tag, i);
burnf  = @(i) sprintf('chain%s_%02d_burn.mat', tag, i);
updf   = ['refupd' tag '.mat'];
ckf    = @(i, ph) sprintf('ckpt%s_%02d_%s.mat', tag, i, ph);   % checkpoints (run_pt_LR; deleted when done)
pt.ckpt_every = optget(opt, 'ckpt_every', 100);
rmck   = @rm_if_exists;
if ~upd_
    parfor (q = 1:numel(runs), parw_)
        ic = runs(q);
        if sample_
            [sv, st] = run_chain_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 1000 + ic, ic, lin0);
            save_chain(chainf(ic), sv, st, ic);
        else
            [sv, st, info] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ic, lin0, setfield(pt, 'ckpt', ckf(ic, 'run')), ic);
            save_pt(chainf(ic), sv, st, ic, info); rmck(ckf(ic, 'run'));
        end
    end
elseif linonly_ && nlu_ > 1 && ~isempty(runs)
    % ---------- burn-in in opt.lin_updates stages, each followed by a ----------
    % ---------- linearization update (lin_update.m), then sampling    ----------
    % stage k: chain<tag>_XX_burn<k>.mat (its last opt.upd_frac saved, T = 1),
    % update k: refupd<tag>_<k>.mat, MINEOS of prep<tag>_lin<k>.mat; every
    % chain continues from its own states with the new linearization
    stgf = @(i, k) sprintf('chain%s_%02d_burn%d.mat', tag, i, k);
    nk = diff(round(linspace(0, pt.burn, nlu_ + 1)));    % stage lengths
    cfg_k = cfg; Uall = cell(1, nlu_);
    for k = 1:nlu_
        ptk = pt; ptk.nit = nk(k); ptk.jump_last = nk(k);   % cycle-skip checks over the whole burn-in
        ptk.burn = nk(k) - max(pt.thin, pt.thin*round(opt.upd_frac*nk(k)/pt.thin));
        parfor (q = 1:numel(runs), parw_)
            ic = runs(q);
            if exist(chainf(ic), 'file') == 2 || exist(stgf(ic, k), 'file') == 2, continue; end
            if k == 1, S0 = {}; else, B = load(stgf(ic, k - 1)); S0 = B.Sf; end
            [sv1, ~, info1, Sf] = run_pt_LR(rec, cfg_k, mc, Tdisp, fb, nm, zplot, 2000 + 100*(k-1) + ic, lin0, ...
                setfield(ptk, 'ckpt', ckf(ic, sprintf('stage%d', k))), ic, S0);
            save_burn(stgf(ic, k), sv1, info1, Sf); rmck(ckf(ic, sprintf('stage%d', k)));
        end
        Uall{k} = lin_update(arrayfun(@(i) stgf(i, k), 1:mc.nchains, 'UniformOutput', false), ...
            sprintf('refupd%s_%d.mat', tag, k), sprintf('%s_lin%s_%d', opt.refid, tag, k), ...
            sprintf('prep%s_lin%d.mat', tag, k), cfg_k, cfg, rec, opt, mc, zplot, Tdisp, fb, nm, lin0, Ltrue);
        Pu = load(Uall{k}.prepfile); cfg_k = apply_linearization(cfg, Pu.P); clear Pu
    end
    if exist(updf, 'file') ~= 2                    % final update + history (for D4)
        U = Uall{end}; U.stages = cellfun(@(x) rmfield(x, 'lin0'), Uall, 'UniformOutput', false);
        tmpf = tmp_name('refupd'); save(tmpf, 'U'); movefile(tmpf, updf);
    end
    pt2 = pt; pt2.nit = pt.nit - pt.burn; pt2.burn = 0; pt2.jump_last = 0;
    parfor (q = 1:numel(runs), parw_)
        ic = runs(q);
        if exist(chainf(ic), 'file') == 2, continue; end
        B = load(stgf(ic, nlu_)); S0 = B.Sf;
        [sv, st, info] = run_pt_LR(rec, cfg_k, mc, Tdisp, fb, nm, zplot, 3000 + ic, lin0, setfield(pt2, 'ckpt', ckf(ic, 'samp')), ic, S0);
        Lb = []; Lt = []; sw = 0;
        for k = 1:nlu_                                % full trace: all burn-in stages
            B = load(stgf(ic, k));
            Lb = [Lb; B.sv1.L]; Lt = [Lt; B.info1.Ltrace]; sw = sw + B.info1.swap; %#ok<AGROW>
        end
        sv.L = [Lb; sv.L]; info.Ltrace = [Lt; info.Ltrace]; info.swap = info.swap + sw;
        info.ref_update = updf;
        save_pt(chainf(ic), sv, st, ic, info); rmck(ckf(ic, 'samp'));
    end
elseif ~isempty(runs)
    % ---------- phase 1: burn-in with the current reference ----------
    % the last opt.upd_frac of the burn-in is saved (T = 1) for the update
    pt1 = pt; pt1.nit = pt.burn; pt1.jump_last = pt.burn;   % cycle-skip checks over the whole burn-in
    pt1.burn = pt.burn - max(pt.thin, pt.thin*round(opt.upd_frac*pt.burn/pt.thin));
    parfor (q = 1:numel(runs), parw_)
        ic = runs(q);
        if exist(chainf(ic), 'file') == 2 || exist(burnf(ic), 'file') == 2, continue; end
        [sv1, ~, info1, Sf] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ic, lin0, setfield(pt1, 'ckpt', ckf(ic, 'burn')), ic);
        save_burn(burnf(ic), sv1, info1, Sf); rmck(ckf(ic, 'burn'));
    end
    % ---------- reference update: once, by one process ----------
    lockd = [updf '.lock'];
    while exist(updf, 'file') ~= 2
        allb = all(arrayfun(@(i) exist(burnf(i), 'file') == 2, 1:mc.nchains));
        if allb && ~(exist(lockd, 'dir') == 7) && system(sprintf('mkdir "%s" 2>/dev/null', lockd)) == 0
            V = 0; Vh = 0; Ve = 0; ns = 0;
            for i = 1:mc.nchains
                B = load(burnf(i));
                V = V + sum(B.sv1.vsv, 2); Vh = Vh + sum(B.sv1.vsh, 2); Ve = Ve + sum(B.sv1.eta, 2);
                ns = ns + size(B.sv1.vsv, 2);
            end
            U = struct();                         % (U was the group-velocity table of the setup)
            U.vsv = V/ns; U.vsh = Vh/ns; U.eta = []; U.nsamples = ns;
            if mc.eta_free, U.eta = Ve/ns; end
            U.refid = sprintf('%s_upd%s', opt.refid, tag);
            U.prepfile = sprintf('prep%s_burnupd.mat', tag);
            fprintf('reference update: mean of %d late burn-in samples -> %s (MINEOS + idagrn6)\n', ns, U.refid);
            ou_ = opt; if linonly_, ou_.lin_only = true; end   % ('lin': MINEOS only, no idagrn6)
            Pu = update_reference(U.vsv, U.vsh, zplot, cfg, ou_, U.refid, U.prepfile, U.eta);
            U.mode = opt.update_ref;
            if linonly_
                % linearization only: same records and nuisance meaning; report
                % how the states' dispersion and log-likelihood change (each
                % state with its own linearization points)
                cfg_u = apply_linearization(cfg, Pu);
                dL = zeros(mc.nchains, 1); dCmax = dL; dCrms = dL;
                for i = 1:mc.nchains
                    B = load(burnf(i)); s1 = B.Sf{1};
                    [sk, okk] = lr_common('disp', struct(), s1.m, cfg_u, Tdisp, fb, nm, 'both');
                    Ln = -Inf;
                    if okk
                        Ln = 0;
                        for r = 1:numel(rec)
                            [~, Lr] = nuis_relin(rec(r), sk.Cf{rec(r).wave}, s1.lin{r}, cfg_u, s1.sig, s1.rho, mc, nm, false);
                            Ln = Ln + Lr;
                        end
                        dC = 100*[(sk.CL - s1.CL)./s1.CL, (sk.CR - s1.CR)./s1.CR];
                        dC = dC(isfinite(dC)); dCmax(i) = max(abs(dC)); dCrms(i) = sqrt(mean(dC.^2));
                    end
                    dL(i) = Ln - B.info1.Ltrace(end, 1);
                end
                U.dL = dL; U.dCmax = dCmax; U.dCrms = dCrms; U.reverted = false;
                U.lin0 = lin0; U.Ltrue = Ltrue;   % waveform reference unchanged: log L(true) unchanged
                fprintf('linearization update: dispersion change of the end-of-burn-in states (%%) rms %s, max %s\n', ...
                    mat2str(round(100*dCrms')/100), mat2str(round(100*dCmax')/100));
                fprintf('linearization update: logL change of the end-of-burn-in states %s (mean %.1f; reported, not tested)\n', ...
                    mat2str(round(dL')), mean(dL));
                tmpf = tmp_name('refupd'); save(tmpf, 'U'); movefile(tmpf, updf);
                rmdir(lockd); clear Pu cfg_u B
                continue
            end
            [rec_u, cfg_u, Cu_f] = apply_reference(rec, cfg, Pu, Tdisp, fb);
            U.lin0 = spec_refine([], rec_u, Cu_f, sig_true, rho_emp, cfg_u, mc, nm, 30);
            % safeguard: the T = 1 end-of-burn-in states of all chains must fit
            % the data at least as well with the new reference (after moving
            % their nuisance linearization points to it); otherwise keep the
            % old reference
            dL = zeros(mc.nchains, 1);
            for i = 1:mc.nchains
                B = load(burnf(i)); s1 = B.Sf{1};
                [sk, okk] = lr_common('disp', struct(), s1.m, cfg_u, Tdisp, fb, nm, 'both');
                Ln = -Inf;
                if okk
                    lk = spec_refine(U.lin0, rec_u, sk.Cf, s1.sig, s1.rho, cfg_u, mc, nm, opt.upd_ngn);
                    Ln = 0;
                    for r = 1:numel(rec_u)
                        [~, Lr] = nuis_relin(rec_u(r), sk.Cf{rec_u(r).wave}, lk{r}, cfg_u, s1.sig, s1.rho, mc, nm, false);
                        Ln = Ln + Lr;
                    end
                end
                dL(i) = Ln - B.info1.Ltrace(end, 1);
            end
            U.dL = dL; U.reverted = mean(dL) < 0;
            fprintf('reference update: logL change of the end-of-burn-in states %s (mean %.1f)\n', mat2str(round(dL')), mean(dL));
            if U.reverted
                fprintf('reference update REJECTED (the states fit worse): keeping the reference of %s\n', ...
                    ifempty_(opt.ref_prep, opt.prepfile));
                U.prepfile = ifempty_(opt.ref_prep, opt.prepfile); U.refid = 'not updated';
                U.lin0 = lin0; U.Ltrue = Ltrue;
            else
                U.Ltrue = NaN;
                if any(isfinite(Ctrue_f{1}(:))), U.Ltrue = fit_nuisance_LR(rec_u, Ctrue_f, sig_true, rho_emp, cfg_u, mc, nm); end
                fprintf('logL of the true structure with the updated reference: %.1f (was %.1f)\n', U.Ltrue, Ltrue);
            end
            tmpf = tmp_name('refupd'); save(tmpf, 'U'); movefile(tmpf, updf);
            rmdir(lockd); clear Pu rec_u cfg_u Cu_f B
        else
            fprintf('waiting for the burn-in of all chains/ensembles and the reference update (%s)\n', updf);
            pause(60);
        end
    end
    % ---------- phase 2: sampling with the updated reference (fixed) ----------
    U = load(updf); U = U.U;
    Pu = load(U.prepfile); Pu = Pu.P;
    lino_ = isfield(U, 'mode') && strcmp(U.mode, 'lin');
    if lino_                                      % linearization only
        rec_u = rec; cfg_u = apply_linearization(cfg, Pu);
    else
        [rec_u, cfg_u] = apply_reference(rec, cfg, Pu, Tdisp, fb);
    end
    clear Pu
    pt2 = pt; pt2.nit = pt.nit - pt.burn; pt2.burn = 0; pt2.jump_last = 0;
    parfor (q = 1:numel(runs), parw_)
        ic = runs(q);
        if exist(chainf(ic), 'file') == 2, continue; end
        B = load(burnf(ic)); S0 = B.Sf;
        for k = 1:numel(S0)                       % linearization points at the new reference
            [sk, okk] = lr_common('disp', struct(), S0{k}.m, cfg_u, Tdisp, fb, nm, 'both');
            if ~okk, error('state %d of chain %d has no valid dispersion with the new reference', k, ic); end
            if ~lino_                             % (linearization only: own points stay valid)
                S0{k}.lin = spec_refine(U.lin0, rec_u, sk.Cf, S0{k}.sig, S0{k}.rho, cfg_u, mc, nm, opt.upd_ngn);
            end
        end
        [sv, st, info] = run_pt_LR(rec_u, cfg_u, mc, Tdisp, fb, nm, zplot, 3000 + ic, U.lin0, setfield(pt2, 'ckpt', ckf(ic, 'samp')), ic, S0);
        sv.L = [B.sv1.L; sv.L];                   % full trace (burn-in: old reference)
        info.Ltrace = [B.info1.Ltrace; info.Ltrace];
        info.swap = info.swap + B.info1.swap;
        info.ref_update = updf;
        save_pt(chainf(ic), sv, st, ic, info); rmck(ckf(ic, 'samp'));
    end
end
done = arrayfun(@(i) exist(sprintf('chain%s_%02d.mat', tag, i), 'file') == 2, 1:mc.nchains);
if all(done)
    D4_step2_plot_results
else
    fprintf('Chains/ensembles %s still missing; then run D4_step2_plot_results.\n', mat2str(find(~done)));
end
