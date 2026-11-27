% LOBE_SCAN  logL of the line data vs a uniform relative shift of one branch's phase velocity
% (all periods), the other branches at the true dispersion. Two setups:
%  oracle : setup_line1.mat  (reference = truth, B/tau/dth fixed at the truth)
%  2b     : setup_line2b.mat (step-1 reference, tau fixed at the truth, B free, dth ~0)
% sigma = sigma_true and sigma*3 (the posterior noise level), rho = empirical.
addpath('functions');
eps_ = (-3:0.025:3)'/100; ne = numel(eps_);
setups = {'setup_line1.mat', 'oracle'; 'setup_line2b.mat', '2b'};
R = struct();
for is = 1:2
    S = load(setups{is, 1});
    rec = S.rec; cfg = S.cfg; mc = S.mc; nm = S.nm; Tdisp = S.Tdisp; fb = S.fb;
    C0 = {S.CtrueL, S.CtrueR};
    wave = [rec.wave];
    for sf = [1 3]
        sig = S.sig_true*sf; rho = S.rho_emp;
        L = NaN(ne, nm, 2);
        for wv = 1:2
            r = find(wave == wv);
            for n = 1:nm
                for ie = 1:ne
                    C = C0{wv}; C(:, n) = C(:, n)*(1 + eps_(ie));
                    Cf = disp_to_band(C, Tdisp, fb);
                    Lr = 0;
                    for q = r
                        [~, Lq] = nuis_relin(rec(q), Cf, S.lin0{q}, cfg, sig, rho, mc, nm, false);
                        Lr = Lr + Lq;
                    end
                    L(ie, n, wv) = Lr;
                end
                fprintf('%s sigma x%d %s branch %d done\n', setups{is, 2}, sf, 'LR'(wv), n-1);
            end
        end
        R.(sprintf('%s_s%d', strrep(setups{is,2}, '2b', 'b2'), sf)) = L;
    end
end
R.eps = eps_;
save('lobe_scan.mat', '-struct', 'R');
