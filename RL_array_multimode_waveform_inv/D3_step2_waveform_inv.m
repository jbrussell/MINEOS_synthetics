% D3_STEP2_WAVEFORM_INV  Step 2, part 2 (driver): joint Love + Rayleigh
% multimode array waveform inversion (transdimensional Bayesian, rjMcMC)
% for Vsv, Vsh and xi, with the step-1 model as the reference model.
%
% Sampler (opt.sampler in parameters_RL.m):
%   'plain' : opt.nchains independent rjMcMC chains from random prior draws
%             (run_chain_LR.m); time shifts and arrival angles are sampled,
%             the spectral amplitudes B_n(w) are marginalized
%   'pt'    : opt.nens parallel-tempering ensembles of opt.K temperatures up
%             to opt.Tmax (run_pt_LR.m); B_n(w), tau_n and dth_n(w) are all
%             marginalized, by linearizing tau_n and dth_n(w) about a pilot
%             run (pt_pilot_LR.m + spec_refine.m); only Vsv, xi and the
%             noise parameters are tempered
%
%   >> D3_step2_waveform_inv                       % all chains/ensembles here
%   $ CHAINS="1 2" octave-cli --eval D3_step2_waveform_inv &      % plain, in parallel
%   $ ENSEMBLES="1" octave-cli --eval "opt.sampler='pt'; D3_step2_waveform_inv" &
%
% Output: setup<tag>.mat, chain<tag>_XX.mat (tag = '_plain' or '_pt' by
% default); when all chains/ensembles exist, the figures (D4).
% Previous: D2_step2_prep_mineos_idagrn.m. Next: D4_step2_plot_results.m

close all;
parameters_RL;
sub_setup_records;
tag = opt.tag;
setupf = ['setup' tag '.mat'];
pt_ = strcmp(opt.sampler, 'pt');
if pt_
    pt.K = opt.K; pt.Tmax = opt.Tmax; pt.T = logspace(0, log10(pt.Tmax), pt.K);
    pt.nit = opt.pt_nit; pt.burn = opt.pt_burn; pt.thin = 5;
    mc.nchains = opt.nens; mc.nit = pt.nit; mc.burn = pt.burn;
else
    pt = [];
end
lin = [];

if exist(setupf, 'file') == 2
    load(setupf, 'Ltrue', 'lin');
else
    if pt_
        % linearization point of tau_n and dth_n(w): pilot chain + Gauss-Newton
        lin = pt_pilot_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 999, opt.pilot_nit);
        Cp = {disp_to_band(lin.pilot_CL, Tdisp, fb), disp_to_band(lin.pilot_CR, Tdisp, fb)};
        lin = spec_refine(lin, rec, Cp, lin.sig, lin.rho, cfg, mc, mc.spec_K, mc.spec_Ka, mc.dthp);
    else
        lin = [];
    end
    % log-likelihood of the true structure, with nuisance parameters fitted
    % to it (idagrn data have no "true" nuisance parameters)
    Ltrue = fit_nuisance_LR(rec, Ctrue_f, sig_true, rho_emp, cfg, mc, nm);
    fprintf('logL of the true structure (fitted nuisance, sigma_true): %.1f\n', Ltrue);
    tmpf = sprintf('setup_%d.mat', round(1e6*rand));
    save(tmpf, 'tru', 'refp', 'cfg', 'mc', 'Tdisp', 'CtrueL', 'CtrueR', 'CrefL', 'CrefR', 'ClinL', 'ClinR', ...
         'zplot', 'sig_true', 'rho_emp', 'rec', 'sta', 'Ltrue', 'fb', 'nm', 'outdir', 'tag', 'pd', ...
         'opt', 'ev', 'lin');
    movefile(tmpf, setupf);
end

if pt_
    runs = 1:mc.nchains; envc = getenv('ENSEMBLES');
else
    runs = 1:mc.nchains; envc = getenv('CHAINS');
end
if ~isempty(envc), runs = str2num(envc); end %#ok<ST2NM>
parfor q = 1:numel(runs)
    ic = runs(q);
    if pt_
        [sv, st, info] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ic, lin, pt, ic);
        save_pt(sprintf('chain%s_%02d.mat', tag, ic), sv, st, ic, info);
    else
        [sv, st] = run_chain_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 1000 + ic, ic);
        save_chain(sprintf('chain%s_%02d.mat', tag, ic), sv, st, ic);
    end
end
done = arrayfun(@(i) exist(sprintf('chain%s_%02d.mat', tag, i), 'file') == 2, 1:mc.nchains);
if all(done)
    D4_step2_plot_results
else
    fprintf('Chains/ensembles %s still missing; then run D4_step2_plot_results.\n', mat2str(find(~done)));
end
