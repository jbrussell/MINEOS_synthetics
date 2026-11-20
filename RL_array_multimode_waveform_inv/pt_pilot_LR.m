function lin = pt_pilot_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, seed, nit)
% PT_PILOT_LR  Starting linearization point for run_pt_LR (version 4): a
% short ordinary chain (run_chain_LR, with the spectral nuisance model)
% gives tau, the constant arrival-angle term dth0 and the higher angle
% terms dthh for every record, plus the noise sigma and r; spec_refine
% then polishes tau0, dth(w) and the spectral amplitudes by Gauss-Newton.
if nargin < 9, nit = 1500; end
mcp = mc; mcp.nit = nit; mcp.burn = round(2*nit/3); mcp.thin = 5;
[svp, stp] = run_chain_LR(rec, cfg, mcp, Tdisp, fb, nm, zplot, seed, 0);
lin.tau0 = mean(svp.tau, 3);
D = mean(svp.dthc, 4);                       % nrec x nm x Ka
lin.dth0 = D(:, :, 1); lin.dthh = D(:, :, 2:end);
lin.b0 = ones(numel(rec), nm);
lin.pilot_CL = stp.CL; lin.pilot_CR = stp.CR;
lin.sig = stp.sig; lin.rho = stp.rho;
end
