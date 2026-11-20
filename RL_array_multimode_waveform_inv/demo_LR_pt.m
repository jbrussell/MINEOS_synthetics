% DEMO_LR_PT  Parallel-tempering version of demo_LR (same data).
% Per-record spectral amplitudes B_n(w), time shifts tau_n and
% frequency-dependent arrival angles dth_n(w) are all marginalized
% analytically: B_n exactly, tau_n and dth_n(w) by linearizing about a
% pilot run refined with spec_refine (event_G_spec). Only Vsv, xi and the
% noise parameters are tempered. Several independent ensembles give R-hat.
%
%   >> demo_LR_pt
%   $ ENSEMBLES="1" octave-cli --eval demo_LR_pt &
%   $ ENSEMBLES="2" octave-cli --eval demo_LR_pt &
%
% Options as in demo_LR, plus K (temperatures), Tmax, pt_nit, pt_burn, nens,
% pilot_nit.

close all;
if ~exist('opt', 'var'), opt = struct(); end
if ~isfield(opt, 'tag'), opt.tag = '_pt'; end
if ~isfield(opt, 'outdir'), opt.outdir = 'figs_pt'; end
setup_LR;
tag = opt.tag;
pt.K = optget(opt, 'K', 8);   pt.Tmax = optget(opt, 'Tmax', 20);
pt.T = logspace(0, log10(pt.Tmax), pt.K);
pt.nit = optget(opt, 'pt_nit', 2500); pt.burn = optget(opt, 'pt_burn', 1200); pt.thin = 5;
mc.nchains = optget(opt, 'nens', 2); mc.nit = pt.nit; mc.burn = pt.burn;

% ---- linearization point: pilot chain + Gauss-Newton refinement (cached) ----
linf = ['lin_pilot' tag '.mat'];
if exist(linf, 'file') == 2, load(linf, 'lin');
else
    lin = pt_pilot_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 999, optget(opt, 'pilot_nit', 1500));
    Cp = {disp_to_band(lin.pilot_CL, Tdisp, fb), disp_to_band(lin.pilot_CR, Tdisp, fb)};
    lin = spec_refine(lin, rec, Cp, lin.sig, lin.rho, cfg, mc, mc.spec_K, mc.spec_Ka, mc.dthp);
    tmpf = sprintf('lin_%d.mat', round(1e6*rand)); save(tmpf, 'lin'); movefile(tmpf, linf);
end

% ---- log L of the true structure on the same (linearized-marginal) scale ----
% (idagrn data have no "true" nuisance parameters: they are fitted)
Ltrue = fit_nuisance_LR(rec, Ctrue_f, sig_true, rho_emp, cfg, mc, nm);
fprintf('linearized-marginal logL of the true structure (fitted nuisance): %.1f\n', Ltrue);
tmpf = sprintf('setup_%d.mat', round(1e6*rand));
save(tmpf, 'tru', 'refp', 'cfg', 'mc', 'Tdisp', 'CtrueL', 'CtrueR', 'CrefL', 'CrefR', 'ClinL', 'ClinR', ...
     'zplot', 'sig_true', 'rho_emp', 'rec', 'sta', 'Ltrue', 'fb', 'nm', 'outdir', 'tag', 'pt', 'lin', 'pc', 'opt', 'ev');
movefile(tmpf, ['setup' tag '.mat']);

ens = 1:mc.nchains;
envc = getenv('ENSEMBLES');
if ~isempty(envc), ens = str2num(envc); end %#ok<ST2NM>
parfor q = 1:numel(ens)
    ie = ens(q);
    [sv, st, info] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 2000 + ie, lin, pt, ie);
    save_pt(sprintf('chain%s_%02d.mat', tag, ie), sv, st, ie, info);
end
done = arrayfun(@(i) exist(sprintf('chain%s_%02d.mat', tag, i), 'file') == 2, 1:mc.nchains);
if all(done), plot_LR; else
    fprintf('Ensembles %s still missing; then run plot_LR with tag = ''%s''.\n', mat2str(find(~done)), tag);
end
