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
% parameters at the reference dispersion (fit_nuisance_LR.m).
%
%   >> D3_step2_waveform_inv                       % all chains/ensembles here
%   $ CHAINS="[]" octave-cli --eval D3_step2_waveform_inv     % setup only
%   $ CHAINS="1 2" octave-cli --eval D3_step2_waveform_inv &  % chains in parallel
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
sample_ = strcmp(opt.nuisance, 'sample');
if pt_ && sample_, error('opt.nuisance = ''sample'' is only available with the plain sampler'); end
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
parfor q = 1:numel(runs)
    ic = runs(q);
    if sample_
        [sv, st] = run_chain_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 1000 + ic, ic, lin0);
        save_chain(sprintf('chain%s_%02d.mat', tag, ic), sv, st, ic);
    else
        [sv, st, info] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ic, lin0, pt, ic);
        save_pt(sprintf('chain%s_%02d.mat', tag, ic), sv, st, ic, info);
    end
end
done = arrayfun(@(i) exist(sprintf('chain%s_%02d.mat', tag, i), 'file') == 2, 1:mc.nchains);
if all(done)
    D4_step2_plot_results
else
    fprintf('Chains/ensembles %s still missing; then run D4_step2_plot_results.\n', mat2str(find(~done)));
end
