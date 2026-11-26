function opt = load_project(project_tag, use_archived_data, results_root)
% opt = load_project(project_tag, use_archived_data, results_root)
%
% The exact options of an archived run (results/<project_tag>/parameters/
% opt_used.mat, written by save_project_results.m), to rerun it or start a
% variant: set opt, then run the drivers (they call parameters_RL, which
% keeps every field already in opt).
%
%   >> opt = load_project('aa1k12');                 % same settings
%   >> D3_step2_waveform_inv
%   >> opt = load_project('aa1k12'); opt.spec_K = 8; opt.project_tag = 'aa1k8';   % a variant
%
% use_archived_data (default false): read the waveform data from the
%   archive (results/<project_tag>/data/<aswms dir>/) instead of
%   opt.aswms_dir, and the reference prep file from results/<project_tag>/
%   step2/ when it was archived.
% Paths of this installation replace those of the archived run (code
% folder, MINEOS folder, event and station files) when the archived ones
% do not exist here. The automatic-value markers are dropped, so the run
% tag, figure folder and project tag of the archived run are kept as set
% values.
if nargin < 2 || isempty(use_archived_data), use_archived_data = false; end
if nargin < 3 || isempty(results_root), results_root = 'results'; end
rd = fullfile(results_root, project_tag);
f = fullfile(rd, 'parameters', 'opt_used.mat');
if exist(f, 'file') ~= 2, error('load_project: %s not found', f); end
S = load(f); opt = S.opt;
for fl = {'tag_auto_', 'outdir_auto_', 'project_tag_auto_'}
    if isfield(opt, fl{1}), opt = rmfield(opt, fl{1}); end
end
% installation paths: let parameters_RL set them again if they do not exist here
for fl = {'mdir', 'eventfile', 'stationfile'}
    if isfield(opt, fl{1}) && ~(exist(opt.(fl{1}), 'file') || exist(opt.(fl{1}), 'dir'))
        opt = rmfield(opt, fl{1});
    end
end
if use_archived_data
    ad = regexprep(opt.aswms_dir, '[/\\]+$', ''); [~, an] = fileparts(ad);
    da = fullfile(rd, 'data', an);
    if exist(da, 'dir') ~= 7, error('load_project: no archived waveform data in %s', da); end
    opt.aswms_dir = [da filesep];
    [~, pn, pe] = fileparts(opt.prepfile);
    pa = fullfile(rd, 'step2', [pn pe]);
    if exist(pa, 'file') == 2, opt.prepfile = pa; end
end
fprintf('options of project %s loaded (run tag %s)%s\n', project_tag, opt.tag, ...
        repmat(sprintf('; waveform data from %s', opt.aswms_dir), 1, use_archived_data));
end
