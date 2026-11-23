% A2_CLEANUP_EVENTS  Data preparation, step 2 (synthetic and real data):
% remove events that are close in time and may interfere with each other,
% and events with too few stations.
%
% Port of a2_cleanup_events.m of ASWMS-Q (github.com/jbrussell/ASWMS-Q)
% (written by Ge Jin), same criteria, run for each component of
% opt.aswms_components:
%   - an event is bad if it has no stations or fewer than min_sta_num (5)
%   - surface-wave time window of each event at the array: origin +
%     min(dist)/5 km/s .. origin + max(dist)/2 km/s (distances to the
%     corners and centre of lalim x lolim); two events whose windows
%     overlap are both bad
%   - the event files of bad events are deleted
% The event map of ASWMS-Q is drawn when the MATLAB Mapping Toolbox is
% available.
% Input/output: <opt.aswms_dir>eventmat/<event id>_<component>.mat
% Next: A3_run_autowinpick.m

parameters_RL;
for icomp = 1:numel(opt.aswms_components)
parameters = aswms_setup(opt, opt.aswms_components{icomp});
workingdir = parameters.workingdir;
eventmatpath = [workingdir,'eventmat/'];

% Setup Error Codes for Bad data
setup_ErrorCode

min_sta_num = parameters.min_sta_num;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize = gridsize*3;
xnode = [lalim(1),mean(lalim),lalim(2)];
ynode = [lolim(1),mean(lolim),lolim(2)];
[xi yi] = meshgrid(xnode,ynode);

clear evotimes evlas evlos evids isgood win_start win_end
matfiles = dir([eventmatpath,'/*_',parameters.component,'.mat']);
for ie = 1:length(matfiles)
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
	disp(event.id)
	evotimes(ie) = event.otime;
	evlas(ie) = event.evla;
	evlos(ie) = event.evlo;
	evids(ie) = {event.id};
	isgood(ie) = 1;
	dist = deg2km(distance(xi,yi,evlas(ie),evlos(ie)));
	win_start(ie) = evotimes(ie) + min(dist(:)/5);
	win_end(ie) = evotimes(ie) + max(dist(:)/2);
	if ~isfield(event,'stadata')
		isgood(ie)=0;
	elseif length(event.stadata)<min_sta_num
		isgood(ie)=0;
	end
end % end of event loop

for ie = 1:length(evlas)
	for je = 1:length(evlas)
		if ie == je
			continue;
		end
		if win_start(ie) > win_start(je) && win_start(ie) < win_end(je)
			isgood(ie) = 0;
			isgood(je) = 0;
		end
		if win_end(ie) > win_start(je) && win_end(ie) < win_end(je)
			isgood(ie) = 0;
			isgood(je) = 0;
		end
	end
end

disp('Bad events:')
badind = find(isgood == 0);
for ie = badind
	disp(evids(ie));
end

com = 'y';
if com == 'y'
for ie = badind
	delete([eventmatpath,char(evids(ie)),'_',parameters.component,'.mat'])   % (this component's file)
end
end

Ngood = sum(isgood == 1);
Nevs = numel(evids);
disp([parameters.component, ': ', num2str(Ngood),'/',num2str(Nevs),' events kept']);

%% Plot map (MATLAB Mapping Toolbox only)
if exist('axesm', 'file') && exist('shaperead', 'file')
figure(3); clf;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
box off;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
plotm(evlas(logical(isgood)),evlos(logical(isgood)),'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
isbad = ~logical(isgood);
if ~isempty(evlas(isbad))
    plotm(evlas(isbad),evlos(isbad),'o','color',[0.4 0 0],'MarkerSize',10,'linewidth',1);
end
plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
for ii = [30 60 90 120]
    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
    plotm(latc,longc,'-k','linewidth',2)
end
end
end % end of loop icomp
