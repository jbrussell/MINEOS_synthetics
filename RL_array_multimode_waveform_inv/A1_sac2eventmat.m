% A1_SAC2EVENTMAT  Data preparation, step 1 (synthetic and real data):
% event-based SAC files -> ASWMS-Q event files.
%
% Port of a1_a_sac2eventmat.m of ASWMS-Q (github.com/jbrussell/ASWMS-Q),
% same steps and criteria, run for each component of opt.aswms_components
% ({Love, Rayleigh}); parameters from opt through aswms_parameters.m
% (overrides in opt.aswms):
%   - events deeper than maxdepth (50 km) or smaller than minMw (5.5) are
%     skipped (from the first SAC file of the event)
%   - traces whose start differs from the event origin by > 1 s: error
%   - SNR (calc_SNR.m: band-passed with parameters.prefilter, energy in
%     the group-velocity window dist/(1.1*refv)..dist/(0.5*refv) over the
%     energy outside) <= snr_tol (3): station skipped
%   - resampled (anti-alias + spline) if resample_delta > DELTA
%   - events with mean epicentral distance outside min_dist_tol ..
%     max_dist_tol (20-115 deg) skipped
% Input : <opt.aswms_dir>sacdata/<event id>/*<component>.sac, event list
%         sacdata/evlist.txt (real data: your SAC files in this layout,
%         with event and station header fields; synthetic: A0)
% Output: <opt.aswms_dir>eventmat/<event id>_<component>.mat, structure
%         event (evla, evlo, evdp, otime, id, Mw, stadata(:): stla, stlo,
%         stel, dist, otime, delta, data, cmp, stnm, filename, snr)
% Next: A2_cleanup_events.m

parameters_RL;
isoverwrite = 1;
iscleanevs = 1;
isplotsnr = 0;
for icomp = 1:numel(opt.aswms_components)
comp = opt.aswms_components{icomp};
parameters = aswms_setup(opt, comp);
dbpath = parameters.dbpath;
eventfile = parameters.eventfile;
workingdir = parameters.workingdir;
outpath = [workingdir,'eventmat/'];
minMw = parameters.minMw;
maxdepth = parameters.maxdepth;
max_dist_tol = parameters.max_dist_tol;
min_dist_tol = parameters.min_dist_tol;
snr_tol = parameters.snr_tol;

if ~exist(outpath, 'dir')
	mkdir(outpath);
end
if iscleanevs                                % (this component's files only)
    delete([outpath, '*_', comp, '.mat']);
end

eventids = textread([dbpath,eventfile],'%s');
for ie = 1:length(eventids)
	matfilename = [outpath,char(eventids(ie)),'_',comp,'.mat'];
	if ~isoverwrite && exist(matfilename)
		disp(['Exist ',matfilename,', Skip!']);
		continue;
	end
	clear event
	datapath = [dbpath, char(eventids(ie)),'/'];
	disp(datapath);
	saclist = dir([datapath,'*',comp,'.sac']);
    if isempty(saclist)
        continue
    end
    iisac = 0;
    is_skip_mag_dep = 0; % Initialize to 0
	for isac = 1:length(saclist)
		% read sac file
		sacfilename = [datapath,saclist(isac).name];
		sac = readsac(sacfilename);
		days = (datenum(sac.NZYEAR-1,12,31)+sac.NZJDAY-datenum(1970,1,1));
		otime = days*24*3600 + sac.NZHOUR*3600+sac.NZMIN*60+sac.NZSEC+sac.NZMSEC/1000;
		% initial the event information by the first sac file
		if isac == 1
            if sac.EVDP/1000 > maxdepth || sac.MAG < minMw
                is_skip_mag_dep = 1;
                disp(['Magnitude or depth exceeded for ',matfilename,', Skip!']);
                break;
            end
			event.evla = sac.EVLA;
			event.evlo = sac.EVLO;
            event.evdp = sac.EVDP/1000;
			event.otime = otime;
			event.dbpath = datapath;
			event.id = char(eventids(ie));
			event.otimestr = datestr(datenum(sac.NZYEAR-1,12,31,sac.NZHOUR,sac.NZMIN,sac.NZSEC)+sac.NZJDAY);
            if isfield(sac,'MAG')
				event.Mw = sac.MAG;
			else
				event.Mw = 9999;
			end
        end

		% wbh check for timing errors
        if abs(event.otime - otime) > 1
            error(['Data start time more than 1 s different from event time for ',matfilename]);
        end

        % calculate snr
        snr = calc_SNR(sac,parameters,isplotsnr);
        if snr <= snr_tol || isnan(snr)
            continue
        end

        % resample the data if necessary
        data_delta = sac.DELTA;
        if parameters.resample_delta > data_delta
            stadata = sac.DATA1;
        	old_taxis = 0:data_delta:(length(stadata)-1)*data_delta;
        	resample_delta = parameters.resample_delta;
        	new_taxis = 0:resample_delta:(length(stadata)-1)*data_delta;
        	% apply anti-alias filter
        	fN = 1/2/data_delta;
        	w_c = 1./2/resample_delta/fN;
        	[b,a] = butter(2,w_c,'low');
        	stadata = filtfilt(b,a,stadata);
        	sac.DATA1 = interp1(old_taxis,stadata,new_taxis,'spline');
        	sac.DELTA = resample_delta;
        end

        iisac = iisac+1;
		% build up event.stadata structure
		event.stadata(iisac).stla = sac.STLA;
		event.stadata(iisac).stlo = sac.STLO;
		event.stadata(iisac).stel = sac.STEL;
		event.stadata(iisac).dist = vdist(sac.STLA,sac.STLO,sac.EVLA,sac.EVLO)/1e3;
		event.stadata(iisac).otime = otime+sac.B;
		event.stadata(iisac).delta = sac.DELTA;
		event.stadata(iisac).data = sac.DATA1;
		event.stadata(iisac).cmp = sac.KCMPNM;
		event.stadata(iisac).stnm = sac.KSTNM;
		event.stadata(iisac).filename = sacfilename;
        event.stadata(iisac).snr = snr;
    end
    try
		if mean([event.stadata(:).dist]) > max_dist_tol
		    disp(['Max distance exceeded for ',char(eventids(ie)),', Skip!']);
		    continue;
		elseif mean([event.stadata(:).dist]) < min_dist_tol
		    disp(['Min distance exceeded for ',char(eventids(ie)),', Skip!']);
		    continue;
		end
    end
    if ~is_skip_mag_dep
        matfilename = [outpath,char(eventids(ie)),'_',comp,'.mat'];
        save(matfilename,'event','-v7')   % (-v7: readable by MATLAB when written by Octave)
        disp(['Save to ',matfilename]);
    end
end % end of loop ie
end % end of loop icomp
