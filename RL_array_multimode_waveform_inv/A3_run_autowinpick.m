% A3_RUN_AUTOWINPICK  Data preparation, step 3 (synthetic and real data):
% automatic surface-wave window of every event.
%
% Port of a3_a_run_autowinpick.m of ASWMS-Q (github.com/jbrussell/ASWMS-Q),
% run for each component of opt.aswms_components, with the ASWMS-Q picker
% auto_win_select.m (functions/aswms; narrow-band Gaussian envelopes at 20
% frequencies over parameters.periods, peaks tracked from the centre
% frequency, group delays fitted against distance, window from
% cycle_before (2) cycles before to cycle_after (5) cycles after the
% arrivals). The window of each event is
%       dist/v1 + t1 .. dist/v2 + t2  (s after event.otime),
% winpara = [v1 t1 v2 t2], saved in the event file (event.winpara) and in
% <opt.aswms_dir>winpara/<event id>_<component>.win; events without a
% window get a .bad file. Stations outside lalim x lolim or the main
% epicentral-distance range get stadata.isgood < 0. Figures (record
% section with the window, plot_win_select.m) with opt.aswms_plot = true,
% saved as winpara/<event id>_<component>.png.
% Input/output: <opt.aswms_dir>eventmat/<event id>_<component>.mat
% Next (optional, MATLAB, interactive): A3b_run_exam_winpara.m; then
% D2_step2_prep_mineos_idagrn.m.

parameters_RL;
isfigure = optget(opt, 'aswms_plot', false);
is_overwrite = 1;
for icomp = 1:numel(opt.aswms_components)
parameters = aswms_setup(opt, opt.aswms_components{icomp});
workingdir = parameters.workingdir;
eventmatpath = [workingdir,'eventmat/'];
outwinpath = [workingdir,'winpara/'];

if ~exist(outwinpath,'dir')
	mkdir(outwinpath)
end

% Setup Error Codes for Bad data
setup_ErrorCode

periods = parameters.periods;
comp = parameters.component;

matfiles = dir([eventmatpath,'/*_',comp,'.mat']);
for ie = 1:length(matfiles)

	clear event eventcs CS
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
	disp(event.id)

	if ~is_overwrite
		filename = [outwinpath,'/',event.id,'_',comp,'.bad'];
		if exist(filename,'file')
			disp(['Exist: ',filename,' Skip!']);
			continue;
		end
		filename = [outwinpath,'/',event.id,'_',comp,'.win'];
		if exist(filename,'file')
			disp(['Exist: ',filename,' Skip!']);
			continue;
		end
	end

	% set up some useful arrays
	stlas = [event.stadata(:).stla];
	stlos = [event.stadata(:).stlo];
	stnms = {event.stadata(:).stnm};
	dists = [event.stadata(:).dist];

	% automatically select the signal window by using ftan method
	disp('Start to picking the window');
	tic
		[winpara event] = auto_win_select(event);
	toc
	if isfigure && length(find([event.stadata.isgood]>0))>1
		plot_win_select(event,periods,winpara);
		print(gcf, '-dpng', '-r100', [outwinpath, event.id, '_', comp, '.png']);   % (saved)
	end
	if length(winpara) ~= 4
		filename = [outwinpath,'/',event.id,'_',comp,'.bad'];
		fp = fopen(filename,'w');
		fprintf(fp,'%f\n',0);
		fclose(fp);
		continue;
	end
	filename = [outwinpath,'/',event.id,'_',comp,'.win'];
	fp = fopen(filename,'w');
	fprintf(fp,'%s %f %f %f %f\n',event.id,winpara(1),winpara(2),winpara(3),winpara(4));
	fclose(fp);
	event.winpara = winpara;
	save([eventmatpath,matfiles(ie).name],'event','-v7');
	disp(winpara);

end % end of event
end % end of loop icomp
