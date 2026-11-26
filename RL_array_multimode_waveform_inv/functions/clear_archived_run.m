function [removed, kept] = clear_archived_run(opt, include_shared, force)
% [removed, kept] = clear_archived_run(opt, include_shared, force)
%
% Remove from the run directory the files and folders of a run that are in
% its project archive <opt.results_root>/<opt.project_tag>/ (written by
% save_project_results.m; the list is project_file_list.m). Driver:
% D99_clear_archived_run.m. Run from the directory of the run.
%
% Model cards in run_MINEOS/CARDS are never removed (they are archived in
% data/CARDS/ but stay in place). A file is removed only if its copy in the archive exists and has the same
% size; a folder only if every file in it (recursively) is there. Anything
% else is kept and listed. Files the archive does not hold (checkpoints
% ckpt<tag>_*, MINEOS scratch in run_MINEOS) are not touched.
%
% include_shared (default false): also remove what other runs may use - the
%   step-1 result (opt.s1_file) and figure folders figs_step1*, the D2
%   record sections (figs_step2), the waveform
%   data (opt.aswms_dir), the reference prep file (opt.prepfile) and
%   figs_compare.png. Off by default, because a
%   later run of another project (or a rerun) reads them from here.
% force (default false): remove without asking; otherwise the list is
%   printed and removal needs a 'y' answer.
%
% removed, kept : cell arrays of the paths removed and of those kept because
%   their archived copy is missing or differs.
if nargin < 2 || isempty(include_shared), include_shared = false; end
if nargin < 3 || isempty(force), force = false; end
rd = fullfile(optget(opt, 'results_root', 'results'), opt.project_tag);
if exist(rd, 'dir') ~= 7, error('clear_archived_run: no archive %s (run D4 or save_project_results first)', rd); end
L = project_file_list(opt);
% never remove model cards from the MINEOS folder (run_MINEOS/CARDS): other
% runs, D5/D6 and MINEOS itself read them
cd_ = fullfile(opt.mdir, 'CARDS');
L = L(~strncmp({L.pat}, cd_, numel(cd_)));
cand = {}; kept = {}; why = {}; skipped = {}; folders = {};
for i = 1:numel(L)
    d = dir(L(i).pat);
    src = fileparts(L(i).pat);
    for k = 1:numel(d)
        if any(strcmp(d(k).name, {'.', '..'})), continue; end
        s = fullfile(src, d(k).name);
        if L(i).sub && ~d(k).isdir, continue; end
        if ~L(i).shared || include_shared
            if same_(s, fullfile(rd, L(i).dst, d(k).name)), cand{end+1} = s; %#ok<AGROW>
            else, kept{end+1} = s; why{end+1} = why_(s, fullfile(rd, L(i).dst, d(k).name)); end %#ok<AGROW>
        else
            skipped{end+1} = s; %#ok<AGROW>
        end
    end
    % a folder whose contents are archived (the figure folder, the waveform
    % data folder): removed too once it is empty
    if ~L(i).sub && ~isempty(regexp(L(i).pat, '[/\\]\*$', 'once')) && (~L(i).shared || include_shared)
        folders{end+1} = src; %#ok<AGROW>
    end
end
removed = {};
fprintf('project %s (run tag %s), archive %s\n', opt.project_tag, opt.tag, rd);
if ~isempty(skipped)
    fprintf('kept, may be used by other runs (include_shared = true removes them): %d item(s)\n', numel(skipped));
    fprintf('    %s\n', skipped{:});
end
if ~isempty(kept)
    fprintf('kept, the archive has no identical copy: %d item(s)\n', numel(kept));
    for k = 1:numel(kept), fprintf('    %s (%s)\n', kept{k}, why{k}); end
    fprintf(['  (changed after the run was archived, e.g. D4 rerun with opt.save_results = false: to archive\n' ...
             '   the current version and then remove it, run  parameters_RL; save_project_results(opt)  and D99 again)\n']);
end
if isempty(cand), fprintf('nothing to remove\n'); return; end
fprintf('to remove (identical copy in the archive): %d item(s)\n', numel(cand));
fprintf('    %s\n', cand{:});
if ~force
    a = input('remove them? [y/N] ', 's');
    if ~any(strcmpi(strtrim(a), {'y', 'yes'})), fprintf('nothing removed\n'); return; end
end
for k = 1:numel(cand)
    if exist(cand{k}, 'dir') == 7, rmdir(cand{k}, 's'); else, delete(cand{k}); end
    removed{end+1} = cand{k}; %#ok<AGROW>
end
for k = 1:numel(folders)
    f = folders{k};
    if exist(f, 'dir') == 7
        e = dir(f); e = e(~ismember({e.name}, {'.', '..'}));
        if isempty(e), rmdir(f); removed{end+1} = f; end %#ok<AGROW>
    end
end
fprintf('removed %d item(s)\n', numel(removed));
end

function tf = same_(s, a)
% s (file or folder) has an identical-size copy at a (recursively)
if exist(s, 'dir') == 7
    tf = exist(a, 'dir') == 7;
    if ~tf, return; end
    e = dir(s);
    for k = 1:numel(e)
        if any(strcmp(e(k).name, {'.', '..'})), continue; end
        if ~same_(fullfile(s, e(k).name), fullfile(a, e(k).name)), tf = false; return; end
    end
else
    ds = dir(s); da = dir(a);
    tf = numel(ds) == 1 && numel(da) == 1 && ~da.isdir && ds.bytes == da.bytes;
end
end

function w = why_(s, a)
% why the copy of s at a does not count as identical
if exist(a, 'file') == 0, w = 'not in the archive'; return; end
if exist(s, 'dir') == 7, w = 'archived folder differs'; return; end
ds = dir(s); da = dir(a);
if ds.datenum > da.datenum, t = 'newer here'; else, t = 'older here'; end
w = sprintf('%d bytes here, %d in the archive; %s', ds.bytes, da.bytes, t);
end
