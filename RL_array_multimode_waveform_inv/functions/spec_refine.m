function linc = spec_refine(linc, rec, Cf, sig, rho, cfg, mc, nm, ngn)
% linc = spec_refine(linc, rec, Cf, sig, rho, cfg, mc, nm, ngn)
% Refine the linearization points of the nuisance parameters (tau_n(w),
% dth_n(w), B_n(w)) of every record for a FIXED structure Cf = {Love,
% Rayleigh} band dispersion by damped Gauss-Newton steps (nuis_relin).
% linc : {nrec} structs tau0 (nm x Kt), dth0 (nm x Ka), b0 (nm x K); pass
%        [] to start from tau = dth = 0, B = 1.
% ngn  : number of Gauss-Newton steps (default 30).
%
% Continuation: when starting from [] with a frequency-dependent time shift
% (Kt > 1), the first ngn steps keep the higher tau_n(w) terms at zero
% (constant tau_n), and ngn more steps then free them. Freeing all terms
% at once lets the overlapping branches trade delays and Gauss-Newton ends
% in poorer local minima (with PA5 data and the step-1 reference: median
% chi-square of the Love fundamental-mode windows 2.0 vs 1.13 with the
% continuation, at the true structure).
if nargin < 9 || isempty(ngn), ngn = 30; end
nrec = numel(rec);
Kt = size(cfg.Psit, 2);
stages = {mc};
if isempty(linc)
    linc = nuis_zero(nrec, nm, cfg);
    if Kt > 1
        mc1 = mc; mc1.taup = [mc.taup(1), 1e-3*ones(1, Kt - 1)];   % higher terms held at 0
        stages = {mc1, mc};
    end
end
for k = 1:numel(stages)
    for r = 1:nrec
        for itr = 1:ngn
            linc{r} = nuis_relin(rec(r), Cf{rec(r).wave}, linc{r}, cfg, sig, rho, stages{k}, nm, true);
        end
    end
end
end
