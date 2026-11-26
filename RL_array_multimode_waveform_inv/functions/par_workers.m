function M = par_workers(mode, nwork, nruns)
% M = par_workers(mode, nwork, nruns)
%
% Number of workers for the parfor loops over chains/ensembles
% (D3_step2_waveform_inv.m: parfor (q = 1:nruns, M)), and the MATLAB
% parallel pool (Parallel Computing Toolbox) for them:
%   mode = [] (default, opt.parallel): M = Inf, i.e. MATLAB's own rule - an
%          open pool is used, or one is started if the preference "create a
%          parallel pool automatically" is on, otherwise serial
%   mode = true  : a pool of min(nwork, nruns) workers, at most the
%          cluster's NumWorkers (pool_size.m; the rest of the runs queue), is
%          started if none is open (an open pool is reused);
%          M = min(nwork, pool size)
%   mode = false : M = 0 (serial in the client)
% nwork : maximum number of workers (opt.par_workers; [] = nruns)
% In Octave (parfor runs serially) and without the toolbox, M = 0.
if isempty(nwork), nwork = nruns; end
M = 0;
if exist('OCTAVE_VERSION', 'builtin') || nruns < 1, return; end
if isempty(mode), M = Inf; return; end
if ~mode, return; end
if exist('parpool', 'file') == 2 && license('test', 'Distrib_Computing_Toolbox')
    p = gcp('nocreate');
    if isempty(p), p = parpool(pool_size(min(nwork, nruns))); end
    M = min([nwork, nruns, p.NumWorkers]);
    fprintf('parallel: %d chains/ensembles on %d workers%s\n', nruns, M, ...
            repmat(' (the rest queue)', 1, M < nruns));
else
    warning('par_workers: Parallel Computing Toolbox not available; running serially');
end
end
