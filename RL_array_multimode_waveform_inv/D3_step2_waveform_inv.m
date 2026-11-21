% D3_STEP2_WAVEFORM_INV  Step 2, part 2 (driver): joint Love + Rayleigh
% multimode array waveform inversion (transdimensional Bayesian, rjMcMC)
% for Vsv, Vsh and xi, with the step-1 model as the reference model.
%
% Sampler (opt.sampler in parameters_RL.m):
%   'plain' : opt.nchains independent rjMcMC chains from random prior draws
%   'pt'    : opt.nens parallel-tempering ensembles of opt.K temperatures up
%             to opt.Tmax (swaps between neighbouring temperatures)
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
%   'none': the reference of opt.prepfile (step 1) throughout.
% opt.ref_prep: OPTIONAL prep file of another reference model to start from
%   (used by D6_step2_update_reference_and_rerun.m).
% The data, windows and normalizations are those of the step-1 reference
% in all cases, so likelihoods stay comparable.
%
%   >> D3_step2_waveform_inv                       % all chains/ensembles here
%   $ CHAINS="[]" octave-cli --eval D3_step2_waveform_inv     % setup only
%   $ CHAINS="1 2" octave-cli --eval D3_step2_waveform_inv &  % chains in parallel
%   $ ENSEMBLES="1" octave-cli --eval "opt.sampler='pt'; D3_step2_waveform_inv" &
%
% Output: setup<tag>.mat, chain<tag>_XX.mat (tag = '_plain' or '_pt' by
% default), with the reference update also chain<tag>_XX_burn.mat (burn-in
% states), refupd<tag>.mat and prep<tag>_burnupd.mat (the new reference);
% when all chains/ensembles exist, the figures (D4).
% Previous: D2_step2_prep_mineos_idagrn.m. Next: D4_step2_plot_results.m

close all;
parameters_RL;
sub_setup_records;
tag = opt.tag;
setupf = ['setup' tag '.mat'];
pt_ = strcmp(opt.sampler, 'pt');
sample_ = strcmp(opt.nuisance, 'sample');
if pt_ && sample_, error('opt.nuisance = ''sample'' is only available with the plain sampler'); end
upd_ = strcmp(opt.update_ref, 'burnin');
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
    Ltrue = fit_nuisance_LR(rec, Ctrue_f, sig_true, rho_emp, cfg, mc, nm);
    fprintf('logL of the true structure (fitted nuisance, sigma_true): %.1f\n', Ltrue);
    tmpf = sprintf('setup_%d.mat', round(1e6*rand));
    save(tmpf, 'tru', 'refp', 'cfg', 'mc', 'Tdisp', 'CtrueL', 'CtrueR', 'CrefL', 'CrefR', 'ClinL', 'ClinR', ...
         'zplot', 'sig_true', 'rho_emp', 'rec', 'sta', 'Ltrue', 'fb', 'nm', 'outdir', 'tag', 'pd', ...
         'opt', 'ev', 'lin0', 'pt');
    movefile(tmpf, setupf);
end

if pt_, envc = getenv('ENSEMBLES'); else, envc = getenv('CHAINS'); end
runs = 1:mc.nchains;
if ~isempty(envc), runs = str2num(envc); end %#ok<ST2NM>
chainf = @(i) sprintf('chain%s_%02d.mat', tag, i);
burnf  = @(i) sprintf('chain%s_%02d_burn.mat', tag, i);
updf   = ['refupd' tag '.mat'];
if ~upd_
    parfor q = 1:numel(runs)
        ic = runs(q);
        if sample_
            [sv, st] = run_chain_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 1000 + ic, ic, lin0);
            save_chain(chainf(ic), sv, st, ic);
        else
            [sv, st, info] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ic, lin0, pt, ic);
            save_pt(chainf(ic), sv, st, ic, info);
        end
    end
elseif ~isempty(runs)
    % ---------- phase 1: burn-in with the current reference ----------
    % the last opt.upd_frac of the burn-in is saved (T = 1) for the update
    pt1 = pt; pt1.nit = pt.burn;
    pt1.burn = pt.burn - max(pt.thin, pt.thin*round(opt.upd_frac*pt.burn/pt.thin));
    parfor q = 1:numel(runs)
        ic = runs(q);
        if exist(chainf(ic), 'file') == 2 || exist(burnf(ic), 'file') == 2, continue; end
        [sv1, ~, info1, Sf] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ic, lin0, pt1, ic);
        save_burn(burnf(ic), sv1, info1, Sf);
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
            Pu = update_reference(U.vsv, U.vsh, zplot, cfg, opt, U.refid, U.prepfile, U.eta);
            [rec_u, cfg_u, Cu_f] = apply_reference(rec, cfg, Pu, Tdisp, fb);
            U.lin0 = spec_refine([], rec_u, Cu_f, sig_true, rho_emp, cfg_u, mc, nm, 30);
            U.Ltrue = fit_nuisance_LR(rec_u, Ctrue_f, sig_true, rho_emp, cfg_u, mc, nm);
            fprintf('logL of the true structure with the updated reference: %.1f (was %.1f)\n', U.Ltrue, Ltrue);
            tmpf = sprintf('refupd_%d.mat', round(1e6*rand)); save(tmpf, 'U'); movefile(tmpf, updf);
            rmdir(lockd); clear Pu rec_u cfg_u Cu_f B
        else
            fprintf('waiting for the burn-in of all chains/ensembles and the reference update (%s)\n', updf);
            pause(60);
        end
    end
    % ---------- phase 2: sampling with the updated reference (fixed) ----------
    U = load(updf); U = U.U;
    Pu = load(U.prepfile); Pu = Pu.P;
    [rec_u, cfg_u] = apply_reference(rec, cfg, Pu, Tdisp, fb); clear Pu
    pt2 = pt; pt2.nit = pt.nit - pt.burn; pt2.burn = 0;
    parfor q = 1:numel(runs)
        ic = runs(q);
        if exist(chainf(ic), 'file') == 2, continue; end
        B = load(burnf(ic)); S0 = B.Sf;
        for k = 1:numel(S0)                       % linearization points at the new reference
            [sk, okk] = lr_common('disp', struct(), S0{k}.m, cfg_u, Tdisp, fb, nm, 'both');
            if ~okk, error('state %d of chain %d has no valid dispersion with the new reference', k, ic); end
            S0{k}.lin = spec_refine(U.lin0, rec_u, sk.Cf, S0{k}.sig, S0{k}.rho, cfg_u, mc, nm, opt.upd_ngn);
        end
        [sv, st, info] = run_pt_LR(rec_u, cfg_u, mc, Tdisp, fb, nm, zplot, 3000 + ic, U.lin0, pt2, ic, S0);
        sv.L = [B.sv1.L; sv.L];                   % full trace (burn-in: old reference)
        info.Ltrace = [B.info1.Ltrace; info.Ltrace];
        info.swap = info.swap + B.info1.swap;
        info.ref_update = updf;
        save_pt(chainf(ic), sv, st, ic, info);
    end
end
done = arrayfun(@(i) exist(sprintf('chain%s_%02d.mat', tag, i), 'file') == 2, 1:mc.nchains);
if all(done)
    D4_step2_plot_results
else
    fprintf('Chains/ensembles %s still missing; then run D4_step2_plot_results.\n', mat2str(find(~done)));
end
