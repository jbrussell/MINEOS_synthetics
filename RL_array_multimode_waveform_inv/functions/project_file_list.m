function L = project_file_list(opt)
% L = project_file_list(opt)
%
% The files and folders of a run that save_project_results.m copies to the
% project archive <opt.results_root>/<opt.project_tag>/, and that
% D99_clear_archived_run.m removes from the working directory once they are
% there. One entry per source pattern (relative to the run directory):
%   L(i).pat    : dir() pattern of the sources
%   L(i).dst    : destination folder in the archive (relative to it)
%   L(i).sub    : true = each match is a folder copied as <dst>/<name>
%                 (figs_step1*: step1/figs_step1...), false = copied into dst
%   L(i).shared : true = may also be used by other runs (step-1 result and
%                 figures, waveform data, reference prep file, model cards,
%                 figs_compare.png); false = belongs to this run only (its
%                 tag: setup, chains, figures folder, logs, linearization
%                 cards)
% Step-2 patterns end at the tag or continue with '_' or '.', so a project
% whose tag extends this one (e.g. _aa1k12 for _aa1) is not picked up.
tag = opt.tag;
L = struct('pat', {}, 'dst', {}, 'sub', {}, 'shared', {});
% figures and summaries
if isfield(opt, 'outdir') && ~isempty(opt.outdir), L = add_(L, fullfile(opt.outdir, '*'), 'figures', false, false); end
L = add_(L, '*figs_compare.png', 'figures', false, true);
% step 2 results
pats = {['setup' tag '.mat'], ['chain' tag '_*.mat'], ['refupd' tag '.mat'], ['refupd' tag '_*.mat'], ...
        ['prep' tag '_*.mat'], ['nuis_true' tag '.mat'], ['*' tag '.*.log'], ['*' tag '.txt']};
for k = 1:numel(pats), L = add_(L, pats{k}, 'step2', false, false); end
data_ = optget(opt, 'results_data', true);
if data_ && isfield(opt, 'prepfile') && ~isempty(opt.prepfile), L = add_(L, opt.prepfile, 'step2', false, true); end
% step 1
if isfield(opt, 's1_file') && ~isempty(opt.s1_file), L = add_(L, opt.s1_file, 'step1', false, true); end
L = add_(L, 'figs_step1*', 'step1', true, true);
% waveform data and model cards
if data_
    ad = optget(opt, 'aswms_dir', '');
    if ~isempty(ad)
        ad = regexprep(ad, '[/\\]+$', '');
        [~, an] = fileparts(ad);
        L = add_(L, fullfile(ad, '*'), fullfile('data', an), false, true);
    end
    cd_ = fullfile(opt.mdir, 'CARDS');
    L = add_(L, fullfile(cd_, [opt.trueid '.card']), fullfile('data', 'CARDS'), false, true);
    L = add_(L, fullfile(cd_, [opt.s1_final '.card']), fullfile('data', 'CARDS'), false, true);
    L = add_(L, fullfile(cd_, ['RLinv_step1_lin' tag '_*.card']), fullfile('data', 'CARDS'), false, false);
end
end

function L = add_(L, pat, dst, sub, shared)
L(end+1) = struct('pat', pat, 'dst', dst, 'sub', sub, 'shared', shared);
end
