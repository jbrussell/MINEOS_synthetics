function U = lin_update(burnfiles, updf, refid, prepfile, cfg_cur, cfg, rec, opt, mc, zplot, Tdisp, fb, nm, lin0, Ltrue)
% U = lin_update(burnfiles, updf, refid, prepfile, cfg_cur, cfg, rec, opt, mc, zplot, Tdisp, fb, nm, lin0, Ltrue)
%
% One linearization-only reference update (opt.update_ref = 'lin', D3),
% shared by all chains/ensembles: when the burn-in (stage) files of all of
% them exist, ONE process (lock directory <updf>.lock) runs MINEOS (no
% idagrn6: opt.lin_only) for the card of the reference model with Vsv and
% Vsh between zmin and zmax replaced by the mean of their saved late
% samples (T = 1), and saves U in updf; the other processes wait for updf.
% The waveform reference (rec.Uref, rec.cref) is not changed, so every
% chain keeps its own nuisance linearization points.
%
% burnfiles : {nchains} burn-in (stage) files (save_burn.m: sv1, info1, Sf)
% cfg_cur   : configuration with the linearization used in this stage (the
%             change of the states' dispersion and log-likelihood from it
%             to the new one is reported, not tested)
% cfg       : configuration of the setup (its reference card is updated)
% U         : vsv, vsh, eta, nsamples, refid, prepfile, mode = 'lin',
%             dL, dCrms, dCmax (per chain), reverted = false, lin0, Ltrue
lockd = [updf '.lock'];
while exist(updf, 'file') ~= 2
    allb = all(cellfun(@(f) exist(f, 'file') == 2, burnfiles));
    if allb && ~(exist(lockd, 'dir') == 7) && system(sprintf('mkdir "%s" 2>/dev/null', lockd)) == 0
        V = 0; Vh = 0; Ve = 0; ns = 0;
        for i = 1:numel(burnfiles)
            B = load(burnfiles{i});
            V = V + sum(B.sv1.vsv, 2); Vh = Vh + sum(B.sv1.vsh, 2); Ve = Ve + sum(B.sv1.eta, 2);
            ns = ns + size(B.sv1.vsv, 2);
        end
        U = struct();
        U.vsv = V/ns; U.vsh = Vh/ns; U.eta = []; U.nsamples = ns;
        if isfield(mc, 'eta_free') && mc.eta_free, U.eta = Ve/ns; end
        U.refid = refid; U.prepfile = prepfile; U.mode = 'lin';
        fprintf('linearization update: mean of %d late samples -> %s (MINEOS)\n', ns, refid);
        o = opt; o.lin_only = true;
        Pu = update_reference(U.vsv, U.vsh, zplot, cfg, o, U.refid, U.prepfile, U.eta);
        cfg_u = apply_linearization(cfg, Pu);
        dL = zeros(numel(burnfiles), 1); dCmax = dL; dCrms = dL;
        for i = 1:numel(burnfiles)
            B = load(burnfiles{i}); s1 = B.Sf{1};
            [sk, okk] = lr_common('disp', struct(), s1.m, cfg_u, Tdisp, fb, nm, 'both');
            [s0, ok0] = lr_common('disp', struct(), s1.m, cfg_cur, Tdisp, fb, nm, 'both');
            Ln = -Inf;
            if okk && ok0
                Ln = 0;
                for r = 1:numel(rec)
                    [~, Lr] = nuis_relin(rec(r), sk.Cf{rec(r).wave}, s1.lin{r}, cfg_u, s1.sig, s1.rho, mc, nm, false);
                    Ln = Ln + Lr;
                end
                dC = 100*[(sk.CL - s0.CL)./s0.CL, (sk.CR - s0.CR)./s0.CR];
                dC = dC(isfinite(dC)); dCmax(i) = max(abs(dC)); dCrms(i) = sqrt(mean(dC.^2));
            end
            dL(i) = Ln - B.info1.Ltrace(end, 1);
        end
        U.dL = dL; U.dCmax = dCmax; U.dCrms = dCrms; U.reverted = false;
        U.lin0 = lin0; U.Ltrue = Ltrue;
        fprintf('linearization update (%s): dispersion change of the states (%%) rms %s, max %s; logL change %s\n', ...
            refid, mat2str(round(100*dCrms')/100), mat2str(round(100*dCmax')/100), mat2str(round(dL')));
        tmpf = tmp_name('refupd'); save(tmpf, 'U'); movefile(tmpf, updf);
        rmdir(lockd);
    else
        fprintf('waiting for all chains/ensembles and the linearization update (%s)\n', updf);
        pause(60);
    end
end
U = load(updf); U = U.U;
end
