function rd = save_project_results(opt)
% rd = save_project_results(opt)
%
% Archive a finished run in <opt.results_root>/<opt.project_tag>/ (default
% results/<project_tag>/) so that runs with different settings can be found,
% compared and rerun. Called at the end of D4_step2_plot_results.m when
% opt.save_results is true; can also be called by hand after a run
% (parameters_RL; save_project_results(opt)). Existing files of the same
% project are overwritten. Run from the directory of the run (the working
% directory of D1-D4).
%
% Contents of results/<project_tag>/:
%   README.txt          what the run is, where everything is, how to rerun
%   settings.txt        every option of the run (opt), one per line
%   parameters/         parameters_RL.m as it was at the time of the run,
%                       opt_used.mat (the complete resolved opt struct,
%                       including overrides given as a struct before the
%                       drivers: the exact settings of the run)
%   figures/            the figure folder of the run (opt.outdir): all
%                       figures, summary.txt, posterior_summary<tag>.mat,
%                       B_posterior<tag>.mat, dth_posterior<tag>.mat
%   step2/              setup<tag>.mat, chain<tag>_XX*.mat (samples; also
%                       the burn-in files), refupd<tag>*.mat, prep<tag>_*.mat
%                       (linearization updates), nuis_true<tag>.mat, logs
%                       *<tag>*.log, the D2 record sections (figs_step2/),
%                       and (opt.results_data) the reference prep file
%                       opt.prepfile
%   step1/              opt.s1_file (step1_result.mat) and the step-1 figure
%                       folders (figs_step1*) as they are at archiving time
%   data/               (opt.results_data, default true) the waveform data:
%                       a copy of opt.aswms_dir (sacdata/, eventmat/,
%                       winpara/, branchmat/, synthetic truth file) and the
%                       MINEOS cards of the true and step-1 models
%                       (data/CARDS/)
% Rerun: opt = load_project('<project_tag>'); then the drivers (README.txt).
pt = opt.project_tag;
rd = fullfile(optget(opt, 'results_root', 'results'), pt);
tag = opt.tag;
dirs = {'parameters', 'figures', 'step1', 'step2', 'data'};
for k = 1:numel(dirs), mk_(fullfile(rd, dirs{k})); end
log_ = {};

% parameters: the file and the resolved options
pf = which('parameters_RL');
if ~isempty(pf), cp_(pf, fullfile(rd, 'parameters')); end
save(fullfile(rd, 'parameters', 'opt_used.mat'), 'opt', '-v7');
fid = fopen(fullfile(rd, 'settings.txt'), 'w');
fprintf(fid, '%% options of project %s (run tag %s), %s\n', pt, tag, datestr(now, 'yyyy-mm-dd HH:MM'));
write_opt_(fid, opt, 'opt');
% the sampler settings as the run used them (mc of setup<tag>.mat: also the
% defaults set in sub_setup_records.m, e.g. spec_K, sb, tau_prior)
mc = [];
if exist(['setup' tag '.mat'], 'file') == 2
    M_ = load(['setup' tag '.mat'], 'mc'); if isfield(M_, 'mc'), mc = M_.mc; end
end
if ~isempty(mc)
    fprintf(fid, '\n%% sampler settings of the run (mc in setup%s.mat)\n', tag);
    write_opt_(fid, mc, 'mc');
end
fclose(fid);

% figures, step-2 results, step 1, waveform data and model cards: the list of
% project_file_list.m (also used by D99_clear_archived_run.m)
data_ = optget(opt, 'results_data', true);
L = project_file_list(opt);
for i = 1:numel(L)
    dst = fullfile(rd, L(i).dst); mk_(dst);
    if L(i).sub
        d1 = dir(L(i).pat);
        for k = 1:numel(d1)
            if d1(k).isdir && ~any(strcmp(d1(k).name, {'.', '..'}))
                mk_(fullfile(dst, d1(k).name));
                log_{end+1} = cp_(fullfile(fileparts(L(i).pat), d1(k).name, '*'), fullfile(dst, d1(k).name)); %#ok<AGROW>
            end
        end
    else
        log_{end+1} = cp_(L(i).pat, dst); %#ok<AGROW>
    end
end

% README
git_ = '';
try
    [st_, out_] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', fileparts(pf)));
    if st_ == 0, git_ = strtrim(out_); end
catch
end
fid = fopen(fullfile(rd, 'README.txt'), 'w');
fprintf(fid, 'Project %s\n', pt);
fprintf(fid, '  run tag %s, sampler %s, archived %s from %s\n', tag, optget(opt, 'sampler', ''), ...
        datestr(now, 'yyyy-mm-dd HH:MM'), pwd);
if ~isempty(git_), fprintf(fid, '  code version (git commit): %s\n', git_); end
fprintf(fid, '\nKey settings (all of them in settings.txt and parameters/opt_used.mat):\n');
keys = {'trueid', 'aswms_dir', 'nbranch_data', 'nm', 'sampler', 'nens', 'K', 'Tmax', 'pt_nit', 'pt_burn', ...
        'nchains', 'nit', 'burn', 'noise_model', 'noise_levels', 'windows', 'nuisance', 'spec_K', 'spec_Ka', ...
        'spec_Kt', 'update_ref', 'lin_updates', 'smooth_sig', 'smooth_hier', 'vsv_ref_sig', 'xi_ref_hier', ...
        'synth_arrival_shift', 'noise', 's1_cstd', 'prepfile'};
for k = 1:numel(keys)
    if isfield(opt, keys{k}), fprintf(fid, '  %-20s = %s\n', keys{k}, val2str_(opt.(keys{k})));
    elseif ~isempty(mc) && isfield(mc, keys{k}), fprintf(fid, '  %-20s = %s  (mc)\n', keys{k}, val2str_(mc.(keys{k}))); end
end
% synthetic data: how they were made (A0 settings in the truth file; these
% are not in the step-2 options)
tf_ = fullfile(optget(opt, 'aswms_dir', ''), optget(opt, 'truthfile', 'synth_truth.mat'));
if exist(tf_, 'file') == 2
    T_ = load(tf_);
    if isfield(T_, 'truth')
        t_ = T_.truth;
        fprintf(fid, '\nSynthetic data (A0, %s): true model %s, data branches %d, noise %g x peak, seed %g', ...
                tf_, val2str_(t_.trueid), t_.nbranch_data, t_.noise, t_.seed);
        if isfield(t_, 'noise_spread'), fprintf(fid, ', station noise spread %g', t_.noise_spread); end
        if isfield(t_, 'dth_true') && any(t_.dth_true ~= 0)
            fprintf(fid, ',\n  arrival-angle anomalies (deg, by event) %s', mat2str(round(t_.dth_true*100)/100));
        end
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\nResults: figures/summary.txt (statistics), figures/posterior_summary%s.mat\n', tag);
fprintf(fid, '  (pooled posterior and its statistics), figures/*.png, step2/chain%s_XX.mat (samples).\n', tag);
fprintf(fid, '\nFolders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/%s\n', ...
        repmat(', data/ (waveform data and model cards)', 1, data_));
fprintf(fid, '\nRerun (from RL_array_multimode_waveform_inv/):\n');
fprintf(fid, '  >> opt = load_project(''%s'');          %% the exact options of this run\n', pt);
fprintf(fid, '  >> opt = load_project(''%s'', true);    %% the same, with the archived waveform data (data/)\n', pt);
fprintf(fid, '  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,\n');
fprintf(fid, '  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data\n');
fprintf(fid, '  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.\n');
fprintf(fid, '  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.\n');
fclose(fid);
ncp = sum(cellfun(@(x) numel(x), log_));
fprintf('run archived in %s (%d files)\n', rd, ncp);
end

function mk_(d)
if exist(d, 'dir') ~= 7, mkdir(d); end
end

function n = cp_(pat, dst)
% copy the files (and folders) matching pat to dst; returns the copied names
n = {};
d = dir(pat);
src = fileparts(pat);
for k = 1:numel(d)
    if any(strcmp(d(k).name, {'.', '..'})), continue; end
    s = fullfile(src, d(k).name);
    if d(k).isdir
        copyfile(s, fullfile(dst, d(k).name));
    else
        copyfile(s, dst);
    end
    n{end+1} = d(k).name; %#ok<AGROW>
end
end

function write_opt_(fid, s, pre)
f = fieldnames(s);
for k = 1:numel(f)
    v = s.(f{k});
    if isstruct(v) && numel(v) == 1 && ~isempty(fieldnames(v))
        write_opt_(fid, v, [pre '.' f{k}]);
    else
        fprintf(fid, '%s.%s = %s;\n', pre, f{k}, val2str_(v));
    end
end
end

function s = val2str_(v)
if ischar(v)
    s = ['''' v ''''];
elseif (isnumeric(v) || islogical(v)) && numel(v) <= 64
    if islogical(v) && isscalar(v), if v, s = 'true'; else, s = 'false'; end
    else, s = mat2str(v, 6); end
elseif isnumeric(v) || islogical(v)
    s = sprintf('[%s array, %g to %g]', mat2str(size(v)), min(v(:)), max(v(:)));
elseif iscellstr(v)
    s = ['{' strjoin(cellfun(@(x) ['''' x ''''], v, 'uniformoutput', false), ', ') '}'];
elseif isstruct(v)
    s = sprintf('[struct %s]', mat2str(size(v)));
else
    s = sprintf('[%s]', class(v));
end
end
