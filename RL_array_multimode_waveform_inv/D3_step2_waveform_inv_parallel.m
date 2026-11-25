% D3_STEP2_WAVEFORM_INV_PARALLEL  D3_step2_waveform_inv.m with its chains /
% parallel-tempering ensembles run in parallel in one MATLAB session
% (Parallel Computing Toolbox): a parallel pool with one worker per
% chain/ensemble (at most opt.par_workers) is started if none is open, and
% every phase (burn-in stages, sampling) runs the chains/ensembles at the
% same time (parfor); the setup, the MINEOS linearization/reference updates
% between the stages and the figures (D4) run in the client. Same results
% as the serial runs (each chain/ensemble has its own random seed);
% checkpoints and resuming as in D3.
%
% Each worker runs one ensemble with all its opt.K temperatures, so the
% useful number of workers is opt.nens (default 4). Memory: every worker
% holds a copy of the records and the linearization (of order 100 MB for
% the 4-branch setup).
%
%   >> D3_step2_waveform_inv_parallel
%   >> opt.par_workers = 2; D3_step2_waveform_inv_parallel   % 2 at a time
%
% Without the toolbox, or in Octave, it runs serially (for Octave use
% separate processes: ENSEMBLES="1" octave-cli --eval D3_step2_waveform_inv &).
% Next: D4_step2_plot_results.m (called automatically when all are done).
if ~exist('opt', 'var'), opt = struct(); end
opt.parallel = true;
D3_step2_waveform_inv
