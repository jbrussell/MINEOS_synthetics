% RUN_RECOVERY_TESTS  The six recovery tests of README.md (results/).
% Each run: plain rjMcMC, 4 chains (here run one after the other; with
% Octave, run chains in parallel processes with CHAINS="k", see demo_LR.m).
% The data are always the idagrn6 seismograms of PA5 with branches 0-7,
% except run 4 (branches 0-2 only).
%
%   1  PREM reference, branches 0-2 modeled              (the requested test)
%   2  reference updated to the posterior mean of run 1, branches 0-2
%   3  PA5 reference (= the truth), branches 0-2
%   4  PA5 reference, branches 0-2, DATA with branches 0-2 only
%   5  PA5 reference, branches 0-7 modeled
%   6  PREM reference, branches 0-7 modeled
%
% Run time with 4 Octave processes: ~30 min per run with 3 branches, ~1 h
% with 8; MINEOS + idagrn6 for each new reference: 1-2 min.

P0 = prep_mineos_idagrn(struct());                       % PREM reference, PA5 data (8 branches)
pa5 = P0.card_true;
base = struct('nchains', 4, 'nit', 4000, 'burn', 2000);
runs = {
  '_s1',   struct('prepfile', 'prep_mineos.mat')
  '_s2',   struct('prepfile', 'prep_s2.mat')
  '_ctl',  struct('prepfile', 'prep_pa5ref.mat', 'refid', 'RLinv_pa5ref', 'ref_card', pa5, 'P_true', P0)
  '_ctl3', struct('prepfile', 'prep_ctl3.mat', 'refid', 'RLinv_pa5ref', 'ref_card', pa5, 'nbranch_data', 3)
  '_ctl8', struct('prepfile', 'prep_ctl8.mat', 'refid', 'RLinv_pa5ref', 'ref_card', pa5, 'nm', 8, 'P_true', P0, 'nit', 3000, 'burn', 1500)
  '_s1n8', struct('prepfile', 'prep_s1n8.mat', 'nm', 8, 'P_true', P0, 'nit', 3000, 'burn', 1500)
};
for irun = 1:size(runs, 1)
    tag = runs{irun, 1}; o = runs{irun, 2};
    if strcmp(tag, '_s2')                                 % iterated reference from run 1
        update_reference('_s1', 'prep_s2.mat', 'RLinv_ref_s2');
    end
    prep_mineos_idagrn(o);                                % cached after the first call
    f = fieldnames(base);
    for i = 1:numel(f), if ~isfield(o, f{i}), o.(f{i}) = base.(f{i}); end, end
    o.tag = tag; o.outdir = ['figs' tag];
    opt = o;
    demo_LR
end
compare_runs(runs(:,1)', {'1: PREM ref., n=0-2', '2: iterated ref., n=0-2', '3: PA5 ref., n=0-2', ...
    '4: PA5 ref., data n=0-2', '5: PA5 ref., n=0-7', '6: PREM ref., n=0-7'}, 'results/compare_runs.png');
