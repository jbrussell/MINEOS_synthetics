function [phv, K] = run_mineos_phv_kernels(cardid, wavetype, periods, nmode, do_kernels, maxF)
% [phv, K] = run_mineos_phv_kernels(cardid, wavetype, periods, nmode, do_kernels, maxF)
%
% Run MINEOS (a1_run_mineos_check.m, and a2_mk_kernels_cv.m if kernels are
% requested) for the card ./CARDS/<cardid>.card (with ./CARDS/<cardid>.qmod)
% and return the phase velocities and phase-velocity kernels of mode
% branch nmode at the requested periods.
%
% The card, wave type and periods are passed to parameter_FRECHET.m via
% environment variables (MINEOS_CARDID, MINEOS_TYPE, MINEOS_PERIODS,
% MINEOS_MAXF), so parameter_FRECHET.m does not need to be edited. The a1
% and a2 scripts are run inside a helper function so that their "clear"
% commands do not clear the caller's workspace, and the caller's figures
% are hidden from the "close all" in a2 (a2 still draws its kernel figures
% 62 and 99). Must be called from the run_MINEOS directory.
%
% INPUTS
% cardid     - card name (without .card)
% wavetype   - 'R' (Rayleigh, spheroidal) or 'L' (Love, toroidal)
% periods    - periods (s) [N x 1]
% nmode      - mode branch (0 = fundamental). Kernels: only nmode = 0
%              (the branch is set in a2_mk_kernels_cv.m)
% do_kernels - 1 = also compute kernels (OPTIONAL, default 1)
% maxF       - maximum frequency (mHz) of the mode calculation (OPTIONAL,
%              default: 10% above the highest frequency of interest)
%
% OUTPUTS
% phv - phase velocity (km/s) [N x 1], interpolated from the mode table
%       (elastic phase velocity, column 6 of the .q file)
% K   - kernels at the card nodes, ordered as in read_model_card (center
%       to surface), [Nnode x N] each:
%       K.vsv, K.vsh, K.vpv, K.vph, K.rho, K.eta : dc/dm(r) (MINEOS units,
%       see run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m for the conversion
%       to layer kernels), K.rad: radius (m) [Nnode x 1], K.per: period of
%       the mode used for each kernel (closest mode to each period)
%       Love waves: K.vpv, K.vph, K.eta are zero.
%
% jbrussell 11/2026

if nargin < 5 || isempty(do_kernels)
    do_kernels = 1;
end
if nargin < 6 || isempty(maxF)
    maxF = ceil(1.1*1000/min(periods)) + 0.05;
end
periods = periods(:);
switch upper(wavetype)
    case {'R','S'}
        TYPE = 'S';
    case {'L','T'}
        TYPE = 'T';
    otherwise
        error('wavetype must be R or L');
end
if do_kernels && nmode ~= 0
    error('run_mineos_phv_kernels: kernels are only set up for the fundamental mode (branch in a2_mk_kernels_cv.m)');
end

% Pass the settings to parameter_FRECHET.m
setenv('MINEOS_CARDID', cardid);
setenv('MINEOS_TYPE', TYPE);
setenv('MINEOS_PERIODS', sprintf('%.6f ', periods));
setenv('MINEOS_MAXF', sprintf('%.2f', maxF));
cleanupObj = onCleanup(@() reset_env());

% Dispersion
run_script('a1_run_mineos_check');
parameter_FRECHET;
if TYPE == 'S'
    TYPEID = param.STYPEID;
else
    TYPEID = param.TTYPEID;
end
qfile = [CARDTABLE, param.CARDID, '.', TYPEID, '.q'];
[phv, ~, ~] = readMINEOS_qfile_per(qfile, periods, nmode);
phv = phv(:);

% Kernels
K = [];
if do_kernels
    % a2 starts with "close all": hide the caller's figures from it
    figs = findall(0, 'type', 'figure');
    vis = cell(size(figs));
    for ii = 1:length(figs)
        vis{ii} = get(figs(ii), 'HandleVisibility');
        set(figs(ii), 'HandleVisibility', 'off');
    end
    try
        out = run_script('a2_mk_kernels_cv');
    catch err
        restore_figs(figs, vis);
        rethrow(err);
    end
    restore_figs(figs, vis);
    if TYPE == 'S'
        FRECH = out.FRECH_S;
    else
        FRECH = out.FRECH_T;
    end
    nper = length(FRECH);
    K.rad = FRECH(1).rad(:);
    K.per = [FRECH(:).per];
    flds = {'vsv','vsh','vpv','vph','rho','eta'};
    for k = 1:length(flds)
        K.(flds{k}) = zeros(length(K.rad), nper);
        for ip = 1:nper
            if isfield(FRECH, flds{k})
                K.(flds{k})(:,ip) = FRECH(ip).(flds{k})(:);
            end
        end
    end
    % MINEOS writes the kernels from the surface down; flip to the card order
    if K.rad(1) > K.rad(end)
        K.rad = flipud(K.rad);
        for k = 1:length(flds)
            K.(flds{k}) = flipud(K.(flds{k}));
        end
    end
    % Remove the ascii kernel files
    delete([param.frechet, param.CARDID, '/', param.CARDID, '.', TYPEID, '.*.cv']);
end

end

function out = run_script(name)
% Run a script in this function's own workspace, so that its "clear" only
% clears this workspace, and return the kernel structures it made (if any)
eval(name);
out = struct();
if exist('FRECH_S','var'); out.FRECH_S = FRECH_S; end
if exist('FRECH_T','var'); out.FRECH_T = FRECH_T; end
end

function restore_figs(figs, vis)
for ii = 1:length(figs)
    if ishandle(figs(ii))
        set(figs(ii), 'HandleVisibility', vis{ii});
    end
end
end

function reset_env()
setenv('MINEOS_CARDID', '');
setenv('MINEOS_TYPE', '');
setenv('MINEOS_PERIODS', '');
setenv('MINEOS_MAXF', '');
end
