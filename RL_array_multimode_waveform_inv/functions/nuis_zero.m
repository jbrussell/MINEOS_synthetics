function linc = nuis_zero(nrec, nm, cfg)
% linc = nuis_zero(nrec, nm, cfg)
% Nuisance linearization points tau = 0, dth = 0, B = 1 for nrec records.
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2);
l = struct('tau0', zeros(nm, Kt), 'dth0', zeros(nm, Ka), 'b0', ones(nm, K));
linc = repmat({l}, nrec, 1);
end
