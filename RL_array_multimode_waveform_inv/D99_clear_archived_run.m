% D99_CLEAR_ARCHIVED_RUN  Clear the working directory of a run that has been
% archived: remove the files and folders of the run (project opt.project_tag)
% that save_project_results.m copied to results/<project_tag>/ (D4 does this
% at its end), so that the next run starts from a clean folder. Only items
% with an identical-size copy in the archive are removed; the list is
% printed and removal asks for confirmation (functions/clear_archived_run.m).
%
% By default only what belongs to this run is removed: setup<tag>.mat,
% chain<tag>_*.mat, refupd<tag>*.mat, prep<tag>_*.mat, nuis_true<tag>.mat,
% the logs and summary files of the tag and the figure folder (opt.outdir).
% What other runs may use stays: the step-1 result and figures, the D2
% record sections (figs_step2), the waveform
% data (opt.aswms_dir) and the reference prep file (opt.prepfile);
% opt.clear_shared = true removes those too. Model cards in
% run_MINEOS/CARDS are never removed.
% Checkpoints (ckpt<tag>_*) are not archived and are not removed.
%
%   >> opt.project_tag = 'aa1k12'; D99_clear_archived_run
%   >> opt = load_project('aa1k12'); D99_clear_archived_run      % options of the archived run
%   >> opt.clear_shared = true; D99_clear_archived_run            % also the shared inputs
%   >> opt.clear_force = true; D99_clear_archived_run             % no confirmation
%
% Previous: D4_step2_plot_results.m (archives the run).
parameters_RL;
clear_archived_run(opt, optget(opt, 'clear_shared', false), optget(opt, 'clear_force', false));
