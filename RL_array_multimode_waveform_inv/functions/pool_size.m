function n = pool_size(n)
% n = pool_size(n)
% Size for parpool: the requested n, but at most the NumWorkers of the
% default cluster profile (for the 'local'/'Processes' cluster, by default
% the number of cores), so that parpool does not fail with "You requested
% a minimum of N workers, but the cluster ... allow[s] a maximum of M".
% Runs beyond the pool size queue in parfor. To allow more workers, raise
% NumWorkers in the cluster profile (Parallel > Create and Manage Clusters).
try
    c = parcluster;
    if isfinite(c.NumWorkers) && c.NumWorkers > 0 && n > c.NumWorkers
        fprintf('parallel: %d workers requested, the cluster allows %d: using %d\n', n, c.NumWorkers, c.NumWorkers);
        n = c.NumWorkers;
    end
catch
end
end
