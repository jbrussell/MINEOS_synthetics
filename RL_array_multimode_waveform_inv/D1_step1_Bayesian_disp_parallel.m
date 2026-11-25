% D1_STEP1_BAYESIAN_DISP_PARALLEL  D1_step1_Bayesian_disp.m with its rjMcMC
% chains run in parallel (MATLAB, Parallel Computing Toolbox): in each
% re-linearization stage the opt.s1b_chains chains run at the same time on
% a parallel pool (parfor, s1b_run_chains.m); the MINEOS runs between the
% stages stay serial. Same results as the serial script (same seeds).
%
% Options (in OPT before running, or in parameters_RL.m):
%   opt.s1b_workers : maximum number of workers (default opt.s1b_chains,
%                     4). A pool is started with that many workers if none
%                     is open (parpool); an open pool is reused.
% Without the toolbox, or in Octave, the chains run serially (warning).
%
%   >> D1_step1_Bayesian_disp_parallel
%   >> opt.s1b_workers = 8; opt.s1b_chains = 8; D1_step1_Bayesian_disp_parallel
%
% Next: D2_step2_prep_mineos_idagrn.m
if ~exist('opt', 'var'), opt = struct(); end
opt.s1b_parallel = true;
D1_step1_Bayesian_disp
