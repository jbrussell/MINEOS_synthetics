function out = mineos_modes(mdir, cardid, TYPE, nbranch, maxF, periods, kbranches)
% out = mineos_modes(mdir, cardid, TYPE, nbranch, maxF, periods, kbranches)
%
% Run MINEOS for the card <mdir>/CARDS/<cardid>.card (with its .qmod):
% mode table (a1_run_mineos_check.m), phase and group velocities of
% branches 0..max(kbranches) at the requested periods and, if kbranches is
% not empty, the phase-velocity kernels of those branches (plot_wk,
% frechet, frechet_cv and draw_frechet_gv, as in a2_mk_kernels_cv.m but for
% any branch).
%
% mdir      - path to the run_MINEOS directory
% cardid    - card name (without .card)
% TYPE      - 'S' (spheroidal, Rayleigh) or 'T' (toroidal, Love)
% nbranch   - number of mode branches in the table (N_modes; <0 = all)
% maxF      - maximum frequency of the table (mHz)
% periods   - periods (s) for the phase velocities and kernels
% kbranches - branches for which to compute kernels (e.g. 0:2), [] = none
%
% out.phv, out.grv : [Nper x nb] phase/group velocity (km/s) of branches
%                    0..nb-1 (nb = max(nbranch, max(kbranches)+1, 3)), elastic
%                    (column 6 and 7 of the .q file); NaN where the branch
%                    does not exist
% out.K{ib}        : kernels of branch kbranches(ib), fields vsv, vsh, vpv,
%                    vph, rho, eta [Nnode x Nper] and rad [Nnode x 1] (card
%                    order, center to surface). MINEOS units: dc (km/s) =
%                    sum K .* dr (m) .* dm (m/s, kg/m^3; eta dimensionless).
%                    Toroidal kernels cover the top Nnode card nodes.
% out.qfile, out.table : paths of the .q and .table files
%
% The settings are passed to parameter_FRECHET.m through environment
% variables (MINEOS_CARDID, MINEOS_TYPE, MINEOS_MAXF, MINEOS_NMODES,
% MINEOS_PERIODS); the MINEOS scripts run in <mdir>.
%
% jbrussell 11/2026

if nargin < 7, kbranches = []; end
periods = periods(:);
here = pwd;
cleanupObj = onCleanup(@() restore(here));
cd(mdir);
setenv('MINEOS_CARDID', cardid);
setenv('MINEOS_TYPE', TYPE);
setenv('MINEOS_MAXF', sprintf('%.2f', maxF));
setenv('MINEOS_NMODES', sprintf('%d', nbranch));
setenv('MINEOS_PERIODS', sprintf('%.6f ', periods));

% mode table
run_script('a1_run_mineos_check');

parameter_FRECHET;
if TYPE == 'S', TYPEID = param.STYPEID; else, TYPEID = param.TTYPEID; end
out.qfile = [CARDTABLE, param.CARDID, '.', TYPEID, '.q'];
out.table = [CARDTABLE, param.CARDID, '.', TYPEID, '.table'];
nb = max([3, max(kbranches)+1, nbranch]);           % all branches in the table
out.phv = nan(numel(periods), nb); out.grv = out.phv;
for n = 0:nb-1
    try
        [p, g] = readMINEOS_qfile_per(out.qfile, periods, n);
        out.phv(:, n+1) = p(:); out.grv(:, n+1) = g(:);
    catch
    end
end
out.periods = periods;
out.K = {};
if isempty(kbranches), return; end

% kernels
setenv('GFORTRAN_STDIN_UNIT', '5');
setenv('GFORTRAN_STDOUT_UNIT', '6');
setenv('GFORTRAN_STDERR_UNIT', '0');
CARDID = param.CARDID;
write_plotwk(TYPE, CARDID);
sys(['cat run_plotwk.', lower(TYPE), ' | plot_wk > plot_wk.LOG'], 'plot_wk');
[~, eig_fils] = system(['ls ', CARDTABLE, CARDID, '.', TYPEID, '_1.eig_fix 2>&1 | cat']);
if ~isempty(strfind(eig_fils, 'No such file'))
    write_frechet(TYPE, CARDID, 0, []);
else
    write_frech_chk(0);
end
sys(['cat run_frechet.', lower(TYPE), ' | frechet > frechet.LOG'], 'frechet');
fp = [param.frechet, CARDID, '/'];
for ib = 1:numel(kbranches)
    BR = kbranches(ib);
    delete([fp, CARDID, '.', TYPEID, '.fcv.', num2str(BR)]);
    write_frechcv(TYPE, CARDID, BR);
    sys(['cat run_frechcv.', lower(TYPE), ' | frechet_cv > frechet_cv.LOG'], 'frechet_cv');
    FR = frechcv_asc(TYPE, CARDID, BR);
    K.rad = FR(1).rad(:);
    for f = {'vsv','vsh','vpv','vph','rho','eta'}
        K.(f{1}) = zeros(numel(K.rad), numel(FR));
        if isfield(FR, f{1})
            for ip = 1:numel(FR), K.(f{1})(:, ip) = FR(ip).(f{1})(:); end
        end
    end
    if K.rad(1) > K.rad(end)              % card order: center to surface
        K.rad = flipud(K.rad);
        for f = {'vsv','vsh','vpv','vph','rho','eta'}, K.(f{1}) = flipud(K.(f{1})); end
    end
    out.K{ib} = K;
    delete([fp, CARDID, '.', TYPEID, '.', num2str(BR), '.*.cv']);
end
delete([fp, '*.fcv.*'], [fp, '*.frech']);
delete(['run_plotwk.', lower(TYPE)], ['run_frechet.', lower(TYPE)], ['run_frechcv.', lower(TYPE)], ...
       ['run_frechcv_asc.', lower(TYPE)], '*.LOG');
end

function sys(com, name)
[status, log] = system(com);
if status ~= 0
    error('mineos_modes: something is wrong at %s\n%s', name, log);
end
end

function run_script(name)
% run a script in this function's own workspace (its "clear" only clears this)
eval(name);
end

function restore(here)
cd(here);
for v = {'MINEOS_CARDID','MINEOS_TYPE','MINEOS_MAXF','MINEOS_NMODES','MINEOS_PERIODS'}
    setenv(v{1}, '');
end
setenv('GFORTRAN_STDIN_UNIT', '-1');
setenv('GFORTRAN_STDOUT_UNIT', '-1');
setenv('GFORTRAN_STDERR_UNIT', '-1');
end
