function R = s1b_run_chains(nch, seeds, nworkers, CLd, CRd, sig, cfg, mc, Tdisp, fb, nm, zplot, nit, burn, thin, w)
% R = s1b_run_chains(nch, seeds, nworkers, CLd, CRd, sig, cfg, mc, Tdisp, fb, nm, zplot, nit, burn, thin, w)
%
% The nch independent rjMcMC chains of the Bayesian dispersion step 1
% (run_disp_inv.m, chain k with random seed seeds(k)), in parallel over a
% MATLAB parallel pool (Parallel Computing Toolbox) when nworkers > 0,
% otherwise one after the other. Used by D1_step1_Bayesian_disp.m
% (opt.s1b_parallel; D1_step1_Bayesian_disp_parallel.m).
%
% nworkers : 0 = serial; > 0 = at most that many workers (a pool is
%            started with min(nworkers, nch) workers if none is open; an
%            open pool is used as it is). Without the toolbox, or in
%            Octave (parfor runs serially there), the chains run serially.
% The chains do not share state (each seeds its own generator, rng(seed),
% and the forward model is the linearized MINEOS dispersion in cfg.lin),
% so the results are the same as in serial runs with the same seeds.
% R : {nch x 1} outputs of run_disp_inv.
M = 0;                                               % parfor (k = .., 0): serial
if nworkers > 0 && ~exist('OCTAVE_VERSION', 'builtin')
    if exist('parpool', 'file') == 2 && license('test', 'Distrib_Computing_Toolbox')
        p = gcp('nocreate');
        if isempty(p), p = parpool(min(nworkers, nch)); end
        M = min(nworkers, p.NumWorkers);
    else
        warning('s1b_run_chains: Parallel Computing Toolbox not available; running the %d chains serially', nch);
    end
end
R = cell(nch, 1);
parfor (k = 1:nch, M)
    R{k} = run_disp_inv(CLd, CRd, sig, cfg, mc, Tdisp, fb, nm, zplot, nit, burn, thin, seeds(k), w);
end
end
