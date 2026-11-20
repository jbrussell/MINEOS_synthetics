function S = idagrn_synth(mdir, cardid, TYPE, maxF, nbranch, events, sta, len_hr, dt, branches)
% S = idagrn_synth(mdir, cardid, TYPE, maxF, nbranch, events, sta, len_hr, dt, branches)
%
% Synthetic seismograms with idagrn6 from the MINEOS mode table of
% <mdir>/CARDS/<cardid>.card (made beforehand with mineos_modes.m with the
% same TYPE, maxF and nbranch).
%
% TYPE     - 'T' (toroidal; transverse component) or 'S' (spheroidal; vertical)
% events   - struct array: id, lat, lon, dep (km), fac (multiplication
%            factor), mt = [Mrr Mtt Mpp Mrt Mrp Mtp] (units of
%            1e20 dyne-cm / fac; see IDAGRN/EVT_FILES/README)
% sta      - struct: name (cell), lat, lon, elev (km)
% len_hr   - record length (hours); dt - sample interval (s)
% branches - [] for the full synthetics (all branches in the table,
%            idagrn6_sac), or a list of branches for single-branch
%            synthetics (idagrn6_mask)
%
% S.t : [Nt x 1] time after origin (s)
% S.u : {nev x nb} cell of [Nt x Nsta] displacement (m); nb = 1 for the full
%       synthetics, numel(branches) otherwise
%
% jbrussell 11/2026

here = pwd;
cleanupObj = onCleanup(@() restore(here));
cd(mdir);
setenv('MINEOS_CARDID', cardid);
setenv('MINEOS_TYPE', TYPE);
setenv('MINEOS_MAXF', sprintf('%.2f', maxF));
setenv('MINEOS_NMODES', sprintf('%d', nbranch));
setenv('GFORTRAN_STDIN_UNIT', '5');
setenv('GFORTRAN_STDOUT_UNIT', '6');
setenv('GFORTRAN_STDERR_UNIT', '0');
parameter_FRECHET;
setpath_idagrn;
CARDID = param.CARDID;
if TYPE == 'T', COMP = 'T'; else, COMP = 'Z'; end

% station file
stafile = fullfile(mdir, 'IDAGRN', 'STATION', 'RLinv_stations.stn');
fid = fopen(stafile, 'w');
fprintf(fid, '%d,\n', numel(sta.name));
for j = 1:numel(sta.name)
    fprintf(fid, '%-8s %10.4f %10.4f %8.4f\n', sta.name{j}, sta.lat(j), sta.lon(j), sta.elev(j));
end
fclose(fid);

% branch masks (plot_wk)
if isempty(branches)
    write_plotwk(TYPE, CARDID);
    sys(['cat run_plotwk.', lower(TYPE), ' | plot_wk > plot_wk.LOG'], 'plot_wk');
    nb = 1;
else
    for BR = branches(:)'
        write_plotwk_branch(TYPE, CARDID, BR);
        sys(['cat run_plotwk_branch.', lower(TYPE), ' | plot_wk > plot_wk.LOG'], 'plot_wk');
    end
    nb = numel(branches);
end

S.u = cell(numel(events), nb);
for e = 1:numel(events)
    ev = events(e);
    evtfile = fullfile(mdir, 'IDAGRN', 'EVT_FILES', ['RLinv_', ev.id]);
    fid = fopen(evtfile, 'w');
    fprintf(fid, '%s\n%.3f %.3f %.2f\n%g\n', ev.id, ev.lat, ev.lon, ev.dep, ev.fac);
    fprintf(fid, '%.5f\n', ev.mt);
    fclose(fid);
    for ib = 1:nb
        delete(['*.', COMP, '.sac']);
        if isempty(branches)
            write_idagrn(TYPE, CARDID, evtfile, stafile, len_hr, dt, COMP);
            copyfile(['run_idagrn.', lower(TYPE)], 'idagrn.in');   % idagrn6 reads idagrn.in
            sys('cat idagrn.in | idagrn6_sac > idagrn.LOG', 'idagrn6_sac');
        else
            write_idagrn_branch(TYPE, CARDID, branches(ib), evtfile, stafile, len_hr, dt, COMP);
            copyfile(['run_idagrn_branch.', lower(TYPE)], 'idagrn.in');
            sys('cat idagrn.in | idagrn6_mask > idagrn.LOG', 'idagrn6_mask');
        end
        U = [];
        for j = 1:numel(sta.name)
            X = rdsac([sta.name{j}, '.', COMP, '.sac']);
            if j == 1
                S.t = X.HEADER.B + (0:numel(X.d)-1)'*X.HEADER.DELTA;
                U = zeros(numel(X.d), numel(sta.name));
            end
            U(:, j) = X.d(:);
        end
        S.u{e, ib} = U;
    end
    delete(evtfile);
end
delete(['*.', COMP, '.sac'], stafile, 'idagrn.in', 'idagrn.LOG', 'plot_wk.LOG', ['run_*.', lower(TYPE)]);
end

function sys(com, name)
[status, log] = system(com);
if status ~= 0
    error('idagrn_synth: something is wrong at %s\n%s', name, log);
end
end

function restore(here)
cd(here);
for v = {'MINEOS_CARDID','MINEOS_TYPE','MINEOS_MAXF','MINEOS_NMODES'}
    setenv(v{1}, '');
end
setenv('GFORTRAN_STDIN_UNIT', '-1');
setenv('GFORTRAN_STDOUT_UNIT', '-1');
setenv('GFORTRAN_STDERR_UNIT', '-1');
end
