% A3B_RUN_EXAM_WINPARA  Data preparation, optional step (MATLAB,
% interactive): review and manually fix the automatic windows.
%
% Port of a3_b_run_exam_winpara.m of ASWMS-Q (github.com/jbrussell/ASWMS-Q),
% run for each component of opt.aswms_components, with exam_winpara_phase.m
% (functions/aswms): record section of each event with its window; keys:
% p pick a window, g toggle event.isgood, q next event, c cheat sheet.
% Events with event.isgood = 0 are not used by the inversion. Needs a
% figure window and ginput (MATLAB).
% Input/output: <opt.aswms_dir>eventmat/<event id>_<component>.mat

parameters_RL;
for icomp = 1:numel(opt.aswms_components)
parameters = aswms_setup(opt, opt.aswms_components{icomp});
workingdir = parameters.workingdir;
eventmat_files = dir([workingdir,'eventmat/*_',parameters.component,'.mat']);

for ie=1:length(eventmat_files)
    load(fullfile(workingdir,'eventmat',eventmat_files(ie).name));
	if isfield(event,'winpara')
		if ~isfield(event,'isgood')
			event.isgood = 1;
		end
		disp(eventmat_files(ie).name);
        disp('p for pick window, g for isgood, q for next event');
        disp('c for cheat sheet');
        event = exam_winpara_phase(event);
	end
    save(fullfile(workingdir,'eventmat',eventmat_files(ie).name),'event','-v7');
end
end % end of loop icomp
