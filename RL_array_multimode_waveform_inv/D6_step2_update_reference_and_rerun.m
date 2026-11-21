% D6_STEP2_UPDATE_REFERENCE_AND_RERUN  Step 2, optional outer iteration
% (driver): new reference model from the posterior mean of a finished run,
% then the waveform inversion again about it.
%
%   1. posterior mean Vsv, Vsh (and eta if it is a model parameter) of the
%      run opt.d6_from (default opt.tag);
%   2. new reference card: the run's reference (after its end-of-burn-in
%      update, if any) with these profiles between the Moho and 350 km,
%      Vp and density scaled with the reference ratios
%      (functions/update_reference.m): MINEOS mode tables, phase/group
%      velocities and kernels, and idagrn6_mask reference spectra are
%      recomputed (~10 min); the data are reused;
%   3. D3_step2_waveform_inv with this reference (opt.ref_prep) and the tag
%      opt.d6_tag (default [d6_from '_it2']), then D4 (figures in
%      figs<d6_tag>). The data, windows and normalizations are unchanged,
%      so the log-likelihoods are comparable with the first run.
% The second run uses the sampler settings and parameterization of the
% first (sampler, chains, iterations, eta_free, ...; stored in its setup
% file) and by default does not update its reference again at the end of
% burn-in (opt.d6_update_ref = 'none'; 'burnin' to do so).
% Repeat with opt.d6_from = the new tag for further iterations; stop when the
% posterior mean no longer changes (and D5 shows a small linearization
% error).
%
%   >> D6_step2_update_reference_and_rerun                    % from opt.tag
%   >> opt.d6_from = '_pt'; D6_step2_update_reference_and_rerun
%
% Previous: D5_step2_check_linearization.m

parameters_RL;
if ~isfield(opt, 'd6_from'), opt.d6_from = opt.tag; end
if ~isfield(opt, 'd6_tag'), opt.d6_tag = [opt.d6_from '_it2']; end
if ~isfield(opt, 'd6_update_ref'), opt.d6_update_ref = 'none'; end
from = opt.d6_from; opt_ = opt;
S0 = load(['setup' from '.mat'], 'rec', 'cfg', 'mc', 'zplot', 'Tdisp', 'fb', 'opt');
% the reference the post-burn-in samples of that run used
if exist(['refupd' from '.mat'], 'file') == 2            % updated at the end of burn-in
    U = load(['refupd' from '.mat']); U = U.U; Pu = load(U.prepfile); Pu = Pu.P;
elseif isfield(S0.opt, 'ref_prep') && ~isempty(S0.opt.ref_prep)
    Pu = load(S0.opt.ref_prep); Pu = Pu.P;
else
    Pu = prep_mineos_idagrn(S0.opt);
end
[~, cfg0] = apply_reference(S0.rec, S0.cfg, Pu, S0.Tdisp, S0.fb); clear Pu
eta_ = isfield(S0.mc, 'eta_free') && S0.mc.eta_free;
V = 0; Vh = 0; Ve = 0; ns = 0;
for ic = 1:S0.mc.nchains
    D = load(sprintf('chain%s_%02d.mat', from, ic));
    V = V + sum(D.sv.vsv, 2); Vh = Vh + sum(D.sv.vsh, 2); ns = ns + size(D.sv.vsv, 2);
    if eta_, Ve = Ve + sum(D.sv.eta, 2); end
end
eta_mean = []; if eta_, eta_mean = Ve/ns; end
refid_new = sprintf('%s_ref%s', opt.refid, opt.d6_tag);
prep_new = sprintf('prep%s.mat', opt.d6_tag);
fprintf('D6: reference from the posterior mean of %s (%d samples) -> %s\n', from, ns, refid_new);
update_reference(V/ns, Vh/ns, S0.zplot, cfg0, opt, refid_new, prep_new, eta_mean);

% rerun the inversion about the new reference (same data and windows)
opt = opt_;
% same sampler settings and model parameterization as the first run
f_ = {'sampler', 'nuisance', 'relin', 'nchains', 'nit', 'burn', 'nens', 'K', 'Tmax', ...
      'pt_nit', 'pt_burn', 'eta_free', 'spec_Kt', 'tau_prior', 'dtau_prior', 'N', 'noise'};
for i_ = 1:numel(f_), if isfield(S0.opt, f_{i_}), opt.(f_{i_}) = S0.opt.(f_{i_}); end, end
opt.ref_prep = prep_new; opt.tag = opt.d6_tag; opt.outdir = ['figs' opt.d6_tag];
opt.update_ref = opt.d6_update_ref;
if exist('OCTAVE_VERSION', 'builtin'), clear -x opt; else, clearvars -except opt; end
D3_step2_waveform_inv
