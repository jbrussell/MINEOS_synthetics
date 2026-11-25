% D4_STEP2_PLOT_RESULTS  Step 2, part 3 (driver): figures and summary of a
% finished inversion (setup<tag>.mat, chain<tag>_XX.mat), in opt.outdir
% (functions/sub_plot_results.m), and, when both samplers have been run,
% a comparison of the plain and parallel-tempering posteriors
% (figs_compare.png). Called automatically by D3 when all chains exist.
%
%   >> D4_step2_plot_results                            % parallel tempering (default)
%   >> opt.sampler = 'plain'; D4_step2_plot_results     % plain chains
%
% Previous: D3_step2_waveform_inv.m. Next: D5_step2_check_linearization.m

parameters_RL;
tag = opt.tag;
opt_ = opt;
sub_plot_results;
opt = opt_; clear opt_
tags = {'_plain', '_pt'};
have = cellfun(@(t) exist(['setup' t '.mat'], 'file') == 2 && exist(sprintf('chain%s_01.mat', t), 'file') == 2, tags);
if all(have)
    compare_runs(tags, {'plain rjMcMC', 'parallel tempering'}, 'figs_compare.png');
end
