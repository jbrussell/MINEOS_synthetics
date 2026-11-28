function L = ridge_wave_L(rec, Cf, lin0, cfg, sig, rho, mc, nm)
L = 0;
for r = 1:numel(rec)
    [~, Lr] = nuis_relin(rec(r), Cf{rec(r).wave}, lin0{r}, cfg, sig, rho, mc, nm, false);
    L = L + Lr;
end
end
