function [sv, st, info] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, seed, lin, pt, ie)
% RUN_PT_LR  Parallel-tempering rjMcMC for the joint Love/Rayleigh
% inversion (version 4: independent B-spline profiles for Vsv and xi,
% kx <= kv, no discontinuities, bottom tied to the reference).
% Per-record amplitudes, time shifts and arrival angles are all
% marginalized analytically, by linearizing tau and dtheta about the pilot
% estimates in LIN (event_G_lin; or event_G_spec, with smooth spectral
% amplitudes B_n(w) and frequency-dependent arrival angles dth_n(w), when
% LIN comes from spec_refine), so only structure and the 8 noise
% parameters (sigma and AR(1) correlation r per window type) are tempered.
% Only the T = 1 chain is saved.

rng(seed);
nrec = numel(rec); wave = [rec.wave]'; K = numel(pt.T);
mue = cell(nrec, 1); linc = cell(nrec, 1);
spec = isfield(lin, 'K');      % event_G_spec: B_n(w) and frequency-dependent dth_n(w)
if spec
    [cfg.Phi, cfg.Psi] = spec_basis(cfg.w, lin.K, lin.Ka);
    Gfun = @event_G_spec;
    for r = 1:nrec
        linc{r} = struct('tau0', lin.tau0(r,:), 'dth0', reshape(lin.dth0(r,:,:), nm, lin.Ka), ...
                         'b0', reshape(lin.b0(r,:,:), nm, lin.K));
        [mue{r}, sb2, itau, ith] = spec_layout(linc{r}, nm, lin.K, lin.Ka, mc.sb, mc.tau_prior, lin.dth_prior);
    end
else
    Gfun = @event_G_lin;
    sb2 = [mc.sb^2*ones(2*nm, 1); mc.tau_prior^2*ones(nm, 1); mc.dth_prior^2*ones(nm, 1)];
    itau = 2*nm + (1:nm); ith = (3*nm + (1:nm))';
    for r = 1:nrec
        mue{r} = [repmat([1; 0], nm, 1); -lin.tau0(r,:)'; -lin.dth0(r,:)'];
        linc{r} = struct('tau0', lin.tau0(r,:), 'dth0', lin.dth0(r,:), 'b0', lin.b0(r,:));
    end
end
lk = @(r, Cf, sig, rho) like_rec(Gfun, rec(r), Cf, sig, rho, linc{r}, mue{r}, sb2, cfg, mc.sb);

% ---------- independent prior draws for every temperature ----------
S = cell(K, 1);
for k = 1:K
    while true
        m = lr_common('draw_prior', cfg, mc);
        [s, ok] = lr_common('disp', struct(), m, cfg, Tdisp, fb, nm, 'both');
        if ok, break; end
    end
    s.sig = zeros(1, 4);
    for wv = 1:2, s.sig(2*wv-1:2*wv) = std(rec(find(wave == wv, 1)).d); end
    s.rho = 0.5*ones(1, 4);
    s.G = cell(nrec, 1); s.Le = zeros(nrec, 1);
    for r = 1:nrec, [s.G{r}, s.Le(r)] = lk(r, s.Cf{wave(r)}, s.sig, s.rho); end
    S{k} = s;
end
st0 = rmfield(S{1}, {'M', 'G', 'Le'}); st0.prof = lr_common('profiles', S{1}, zplot, cfg);
st0.tau = lin.tau0; st0.dth = lin.dth0(:,:,1);
if spec, st0.dthc = lin.dth0; end

nsave = floor((pt.nit - pt.burn)/pt.thin);
sv.st0 = st0;
sv.vsv = zeros(numel(zplot), nsave); sv.vsh = sv.vsv; sv.xi = sv.vsv;
sv.CL = zeros(numel(Tdisp), nm, nsave); sv.CR = sv.CL;
sv.kv = zeros(nsave, 1); sv.kx = sv.kv; sv.sig = zeros(nsave, 4); sv.rho = sv.sig;
sv.tau = zeros(nrec, nm, nsave); sv.dth = sv.tau;
Ka = size(ith, 2);                      % sv.dthc: all arrival-angle coefficients
sv.dthc = zeros(nrec, nm, Ka, nsave);
sv.L = zeros(pt.nit, 1);
info.T = pt.T; info.swap = zeros(K-1, 2); info.Ltrace = zeros(pt.nit, K);
info.acc = zeros(2, 4, 2);
isv = 0; tic;
fprintf('PT ensemble %d: %d temperatures, T = %s\n', ie, K, mat2str(round(pt.T*100)/100));

for it = 1:pt.nit
    for k = 1:K
        s = S{k}; Tk = pt.T(k);
        [mn, okp, set, mv] = lr_common('propose', s.m, cfg, mc, Tk);
        is = 1 + strcmp(set, 'xi');
        if okp
            if k == 1, info.acc(is, mv, 1) = info.acc(is, mv, 1) + 1; end
            if is == 2, which = 'love'; rr = find(wave == 1)'; else, which = 'both'; rr = 1:nrec; end
            [Sn, ok] = lr_common('disp', s, mn, cfg, Tdisp, fb, nm, which);
            if ok
                Gn = s.G; Ln = s.Le;
                for r = rr, [Gn{r}, Ln(r)] = lk(r, Sn.Cf{wave(r)}, s.sig, s.rho); end
                if log(rand) < (sum(Ln) - sum(s.Le)) / Tk
                    Sn.G = Gn; Sn.Le = Ln; s = Sn;
                    if k == 1, info.acc(is, mv, 2) = info.acc(is, mv, 2) + 1; end
                end
            end
        end
        % noise move: level sigma or AR(1) correlation r (cached operators)
        sg = s.sig; rh = s.rho; iw = randi(4); okn = true;
        if rand < 0.5, sg(iw) = sg(iw)*exp(mc.step_lsig*randn);
        else, rh(iw) = rh(iw) + mc.step_rho*randn; okn = rh(iw) >= 0 && rh(iw) < mc.rhomax; end
        if okn
            Ln = s.Le;
            for r = find(wave == ceil(iw/2))'
                Ln(r) = event_loglike_ar(rec(r).d, s.G{r}, rec(r).wt, rec(r).first, sg, rh, mue{r}, sb2);
            end
            if log(rand) < (sum(Ln) - sum(s.Le)) / Tk, s.sig = sg; s.rho = rh; s.Le = Ln; end
        end
        S{k} = s; info.Ltrace(it, k) = sum(s.Le);
    end
    % swaps between neighbouring temperatures (alternating parity)
    for p = (1 + mod(it, 2)):2:(K-1)
        info.swap(p, 1) = info.swap(p, 1) + 1;
        if log(rand) < (sum(S{p+1}.Le) - sum(S{p}.Le)) * (1/pt.T(p) - 1/pt.T(p+1))
            tmp = S{p}; S{p} = S{p+1}; S{p+1} = tmp; info.swap(p, 2) = info.swap(p, 2) + 1;
        end
    end
    s = S{1}; sv.L(it) = sum(s.Le);
    if it > pt.burn && mod(it - pt.burn, pt.thin) == 0
        isv = isv + 1;
        pr = lr_common('profiles', s, zplot, cfg);
        sv.vsv(:, isv) = pr(:,1); sv.vsh(:, isv) = pr(:,2); sv.xi(:, isv) = pr(:,3);
        sv.CL(:,:,isv) = s.CL; sv.CR(:,:,isv) = s.CR;
        sv.kv(isv) = numel(s.m.V.zc); sv.kx(isv) = numel(s.m.X.zc); sv.sig(isv,:) = s.sig; sv.rho(isv,:) = s.rho;
        for r = 1:nrec
            [~, xh] = event_loglike_ar(rec(r).d, s.G{r}, rec(r).wt, rec(r).first, s.sig, s.rho, mue{r}, sb2);
            sv.tau(r,:,isv) = lin.tau0(r,:) + xh(itau)';
            dc = reshape(linc{r}.dth0, nm, Ka) + xh(ith);         % nm x Ka
            sv.dthc(r,:,:,isv) = reshape(dc, 1, nm, Ka);
            sv.dth(r,:,isv) = dc(:, 1)';                          % band-mean angle
        end
    end
    if mod(it, 100) == 0
        fprintf('PT %d  it %5d  logL(T=1) %.1f  kv=%d kx=%d  swap%% %s  %.0f s\n', ie, it, ...
            sum(S{1}.Le), numel(S{1}.m.V.zc), numel(S{1}.m.X.zc), ...
            mat2str(round(100*info.swap(:,2)./max(info.swap(:,1),1))'), toc);
    end
end
info.Lmean = mean(info.Ltrace(pt.burn+1:end, :), 1)';
sv.acc = info.acc;
st = rmfield(S{1}, {'M', 'G', 'Le'});
st.tau = mean(sv.tau(:,:,max(1,isv-20):isv), 3); st.dth = mean(sv.dth(:,:,max(1,isv-20):isv), 3);
st.dthc = mean(sv.dthc(:,:,:,max(1,isv-20):isv), 4);
fprintf('PT %d done: swap acceptance %s %%\n', ie, ...
    mat2str(round(100*info.swap(:,2)./max(info.swap(:,1),1))'));
end

function [G, L] = like_rec(Gfun, rc, Cf, sig, rho, lin, mu, sb2, cfg, sb)
G = Gfun(rc, Cf, lin, cfg, sig(2*rc.wave-1:2*rc.wave), sb);
L = event_loglike_ar(rc.d, G, rc.wt, rc.first, sig, rho, mu, sb2);
end
