% DEMO_LR  Joint Love + Rayleigh, radially anisotropic, multimode array
% waveform inversion (NoMelt geometry, 15 real 2012 events, idagrn6 data
% for PA5, MINEOS dispersion, PREM reference) with ordinary rjMcMC.
%
% Vsv(z) and xi(z) = (Vsh/Vsv)^2 are independent transdimensional B-spline
% profiles; Vsh = Vsv*sqrt(xi). Several independent chains start from
% random prior draws.
%
%   >> demo_LR
%   $ CHAINS="1 2" octave-cli --eval demo_LR &     % chains in parallel
%   $ CHAINS="3 4" octave-cli --eval demo_LR &
%
% Options (struct opt before calling): tag, outdir, prepfile (MINEOS and
% idagrn6 cache; the reference model), nit, burn, nchains, ... (setup_LR).
% Iterated reference (stage 2): after a run,
%   >> update_reference('_plain', 'prep_s2.mat', 'RLinv_ref_s2');
%   >> opt = struct('tag', '_s2', 'outdir', 'figs_s2', 'prepfile', 'prep_s2.mat'); demo_LR
%
% Outputs setup<tag>.mat, chain<tag>_XX.mat and, when all chains exist,
% figures and summary.txt in outdir.

close all;
if ~exist('opt', 'var'), opt = struct(); end
if ~isfield(opt, 'tag'), opt.tag = '_plain'; end
if ~isfield(opt, 'outdir'), opt.outdir = 'figs_plain'; end
setup_LR;
tag = opt.tag;
setupf = ['setup' tag '.mat'];
if exist(setupf, 'file') == 2
    load(setupf, 'Ltrue');
else
    % log-likelihood of the true structure, with nuisance parameters fitted
    % to it (the idagrn data have no "true" nuisance parameters)
    Ltrue = fit_nuisance_LR(rec, Ctrue_f, sig_true, rho_emp, cfg, mc, nm);
    fprintf('logL of the true structure (fitted nuisance, sigma_true): %.1f\n', Ltrue);
    tmpf = sprintf('setup_%d.mat', round(1e6*rand));
    save(tmpf, 'tru', 'refp', 'cfg', 'mc', 'Tdisp', 'CtrueL', 'CtrueR', 'CrefL', 'CrefR', 'ClinL', 'ClinR', ...
         'zplot', 'sig_true', 'rho_emp', 'rec', 'sta', 'Ltrue', 'fb', 'nm', 'outdir', 'tag', 'pc', 'opt', 'ev');
    movefile(tmpf, setupf);
end

chains = 1:mc.nchains;
envc = getenv('CHAINS');
if ~isempty(envc), chains = str2num(envc); end %#ok<ST2NM>
parfor q = 1:numel(chains)
    ic = chains(q);
    [sv, st] = run_chain_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, 1000 + ic, ic);
    save_chain(sprintf('chain%s_%02d.mat', tag, ic), sv, st, ic);
end
done = arrayfun(@(i) exist(sprintf('chain%s_%02d.mat', tag, i), 'file') == 2, 1:mc.nchains);
if all(done), plot_LR; else
    fprintf('Chains %s still missing; then run plot_LR with tag = ''%s''.\n', mat2str(find(~done)), tag);
end
