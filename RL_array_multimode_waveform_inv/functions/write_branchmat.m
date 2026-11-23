function write_branchmat(P, opt)
% write_branchmat(P, opt)
%
% Write the single-branch reference synthetics of a prep structure
% (prep_mineos_idagrn.m: P.u_ref_T, P.u_ref_S, branches 0..P.nm-1 of the
% reference card P.refid) in the format of the ASWMS-Q event files
% (A1_sac2eventmat.m), one file per event and component:
%
%   <opt.aswms_dir>branchmat/<P.refid>/<event id>_<component>.mat
%
% with a structure "event" whose fields are those of the data event files,
%   id, evla, evlo, evdp (km), Mw, otime (epoch s of the CMT centroid time:
%   time zero of the synthetics), otimestr, stadata(j): stla, stlo, stel
%   (m), dist (km, vdist), otime (epoch s of the first sample), delta,
%   cmp, stnm, data
% except that stadata(j).data is [Nt x nm], column n+1 = branch n, plus
%   refid, branches (0..nm-1), Tdisp, phv_ref, grv_ref ([numel(Tdisp) x nm],
%   MINEOS phase and group velocities of the reference) and maxF.
% Components: opt.aswms_components ({Love, Rayleigh}).
wd = optget(opt, 'aswms_dir', ['ASWMS' filesep]);
if wd(end) ~= '/' && wd(end) ~= filesep, wd = [wd filesep]; end
comps = optget(opt, 'aswms_components', {'BHT', 'BHZ'});
outdir = [wd 'branchmat' filesep P.refid filesep];
if ~exist(outdir, 'dir'), mkdir(outdir); end
TYs = 'TS'; ep0 = datenum(1970, 1, 1);
for wv = 1:2
    TY = TYs(wv);
    U = P.(['u_ref_' TY]);
    for e = 1:numel(P.ev)
        ev = P.ev(e);
        event = struct();
        event.id = ev.id; event.evla = ev.lat; event.evlo = ev.lon; event.evdp = ev.dep;
        event.Mw = ev.Mw;
        event.otime = (ev.t_cent - ep0)*86400;
        event.otimestr = datestr(ev.t_cent);
        event.refid = P.refid; event.branches = 0:P.nm-1; event.maxF = P.maxF;
        event.Tdisp = P.Tdisp; event.phv_ref = P.(['phv_ref_' TY]); event.grv_ref = P.(['grv_ref_' TY]);
        for j = 1:numel(P.sta.name)
            event.stadata(j).stla = P.sta.lat(j);
            event.stadata(j).stlo = P.sta.lon(j);
            event.stadata(j).stel = 1000*P.sta.elev(j);
            event.stadata(j).dist = vdist(P.sta.lat(j), P.sta.lon(j), ev.lat, ev.lon)/1e3;
            event.stadata(j).otime = event.otime + P.t(1);
            event.stadata(j).delta = P.dt;
            event.stadata(j).cmp = comps{wv};
            event.stadata(j).stnm = P.sta.name{j};
            event.stadata(j).data = cell2mat(cellfun(@(x) x(:, j), U(e, :), 'UniformOutput', false));
        end
        save([outdir ev.id '_' comps{wv} '.mat'], 'event', '-v7');
    end
end
fprintf('reference branch synthetics (%s, branches 0-%d) written to %s\n', P.refid, P.nm - 1, outdir);
end
