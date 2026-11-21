function [sv, st] = run_chain_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, seed, ic, lin0)
% RUN_CHAIN_LR  One rjMcMC chain for the joint Love/Rayleigh inversion
% with independent transdimensional B-spline profiles for Vsv and xi
% (kx <= kv, no discontinuities, bottom tied to the reference), with the
% time shifts and arrival angles SAMPLED with the exact (non-linearized)
% waveforms (opt.nuisance = 'sample'; the default, 'marginal', uses
% run_pt_LR.m with one temperature instead).
%
% Per-record, per-mode nuisance model (event_G_spec, basis-only form):
%   B_n(w)   smooth complex excitation correction on K = mc.spec_K hat
%            functions in log w -- marginalized analytically;
%   tau_n(w) time shift on Kt = mc.spec_Kt Legendre terms in normalized
%            log w -- sampled, prior std mc.taup (s);
%   dth_n(w) arrival-angle anomaly on Ka = mc.spec_Ka Legendre terms --
%            sampled, prior std mc.dthp (deg).
% Noise is AR(1)-correlated with sigma and r per window type, sampled
% hierarchically. The starting model is a random draw from the prior; tau
% and dth start from lin0 (nrec x nm x Kt and nrec x nm x Ka) if given, or 0.

rng(seed);
nrec = numel(rec); wave = [rec.wave]';
K = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2); nb = 2*nm*K;
mu_b = repmat([1; 0], nm*K, 1); sb2 = mc.sb^2*ones(nb, 1);
dthp = mc.dthp(1:Ka); taup = mc.taup(1:Kt);
lp_nuis = @(tau, dc) -0.5*sum(sum(tau.^2 ./ taup.^2)) - 0.5*sum(sum(dc.^2 ./ dthp.^2));
Gfun = @(r, Cf, tau, dc) event_G_spec(rec(r), Cf, struct('tau0', tau, 'dth0', dc, ...
           'b0', ones(nm, K)), cfg, [], [], true);
llk = @(r, G, sig, rho) event_loglike_ar(rec(r).d, G, rec(r).wt, rec(r).first, sig, rho, mu_b, sb2);

% ---------- random start from the prior ----------
while true
    m = lr_common('draw_prior', cfg, mc);
    [st, ok] = lr_common('disp', struct(), m, cfg, Tdisp, fb, nm, 'both');
    if ok, break; end
end
st.tauc = zeros(nrec, nm, Kt); st.dthc = zeros(nrec, nm, Ka);
if nargin >= 10 && ~isempty(lin0)
    for r = 1:nrec
        st.tauc(r,:,:) = reshape(lin0{r}.tau0, 1, nm, Kt); st.dthc(r,:,:) = reshape(lin0{r}.dth0, 1, nm, Ka);
    end
end
tauo = @(r, s_) reshape(s_.tauc(r,:,:), nm, Kt);
st.sig = zeros(1, 4);
for wv = 1:2, st.sig(2*wv-1:2*wv) = std(rec(find(wave == wv, 1)).d); end
st.rho = 0.5*ones(1, 4);
G = cell(nrec, 1); Le = zeros(nrec, 1);
for r = 1:nrec
    G{r} = Gfun(r, st.Cf{wave(r)}, tauo(r, st), reshape(st.dthc(r,:,:), nm, Ka));
    Le(r) = llk(r, G{r}, st.sig, st.rho);
end
st0 = rmfield(st, 'M'); st0.prof = lr_common('profiles', st, zplot, cfg);
st0.dth = st.dthc(:,:,1); st0.tau = st.tauc(:,:,1);

nsave = floor((mc.nit - mc.burn)/mc.thin);
sv.st0 = st0;
sv.vsv = zeros(numel(zplot), nsave); sv.vsh = sv.vsv; sv.xi = sv.vsv;
sv.CL = zeros(numel(Tdisp), nm, nsave); sv.CR = sv.CL;
sv.kv = zeros(nsave, 1); sv.kx = sv.kv; sv.sig = zeros(nsave, 4); sv.rho = sv.sig;
sv.tau = zeros(nrec, nm, nsave); sv.dth = sv.tau;
sv.tauc = zeros(nrec, nm, Kt, nsave);   % all time-shift coefficients
sv.dthc = zeros(nrec, nm, Ka, nsave);   % all arrival-angle coefficients
sv.L = zeros(mc.nit, 1);
acc = zeros(2, 4, 2);                   % [vsv; xi] x [b d m v] x [tried acc]
accn = [0 0]; isv = 0; tic;
fprintf('chain %d: start with %d interior Vsv and %d xi control points\n', ic, numel(m.V.zc), numel(m.X.zc));

for it = 1:mc.nit
    % ---------- structural move ----------
    [mn, okp, set, mv] = lr_common('propose', st.m, cfg, mc, 1);
    is = 1 + strcmp(set, 'xi');
    if okp
        acc(is, mv, 1) = acc(is, mv, 1) + 1;
        % xi changes Vsh, which both Love and Rayleigh waves see: recompute both
        which = 'both'; rr = 1:nrec;
        [Sn, ok] = lr_common('disp', st, mn, cfg, Tdisp, fb, nm, which);
        if ok
            Gn = G; Ln = Le;
            for r = rr
                Gn{r} = Gfun(r, Sn.Cf{wave(r)}, tauo(r, st), reshape(st.dthc(r,:,:), nm, Ka));
                Ln(r) = llk(r, Gn{r}, st.sig, st.rho);
            end
            if log(rand) < sum(Ln) - sum(Le)
                Sn.tauc = st.tauc; Sn.dthc = st.dthc; Sn.sig = st.sig; Sn.rho = st.rho;
                st = Sn; G = Gn; Le = Ln; acc(is, mv, 2) = acc(is, mv, 2) + 1;
            end
        end
    end
    % ---------- per-record nuisance sweep: one tau_n(w) or dth_n(w) coefficient ----------
    for r = 1:nrec
        tau = tauo(r, st); dc = reshape(st.dthc(r,:,:), nm, Ka); n = randi(nm);
        if rand < 0.5, a = randi(Kt); tau(n, a) = tau(n, a) + mc.step_tau*randn;
        else, a = randi(Ka); dc(n, a) = dc(n, a) + mc.step_dth*randn; end
        Gt = Gfun(r, st.Cf{wave(r)}, tau, dc);
        Lt = llk(r, Gt, st.sig, st.rho);
        accn(1) = accn(1) + 1;
        if log(rand) < Lt + lp_nuis(tau, dc) - Le(r) - lp_nuis(tauo(r, st), reshape(st.dthc(r,:,:), nm, Ka))
            st.tauc(r,:,:) = reshape(tau, 1, nm, Kt); st.dthc(r,:,:) = reshape(dc, 1, nm, Ka); G{r} = Gt; Le(r) = Lt;
            accn(2) = accn(2) + 1;
        end
    end
    % ---------- hierarchical noise: level sigma or correlation r ----------
    sg = st.sig; rh = st.rho; iw = randi(4); okn = true;
    if rand < 0.5, sg(iw) = sg(iw)*exp(mc.step_lsig*randn);
    else, rh(iw) = rh(iw) + mc.step_rho*randn; okn = rh(iw) >= 0 && rh(iw) < mc.rhomax; end
    if okn
        rr = find(wave == ceil(iw/2))'; Ln = Le;
        for r = rr, Ln(r) = llk(r, G{r}, sg, rh); end
        if log(rand) < sum(Ln) - sum(Le), st.sig = sg; st.rho = rh; Le = Ln; end
    end
    % ---------- save ----------
    sv.L(it) = sum(Le);
    if it > mc.burn && mod(it - mc.burn, mc.thin) == 0
        isv = isv + 1;
        p = lr_common('profiles', st, zplot, cfg);
        sv.vsv(:, isv) = p(:,1); sv.vsh(:, isv) = p(:,2); sv.xi(:, isv) = p(:,3);
        sv.CL(:,:,isv) = st.CL; sv.CR(:,:,isv) = st.CR;
        sv.kv(isv) = numel(st.m.V.zc); sv.kx(isv) = numel(st.m.X.zc);
        sv.sig(isv,:) = st.sig; sv.rho(isv,:) = st.rho; sv.tau(:,:,isv) = st.tauc(:,:,1);
        sv.tauc(:,:,:,isv) = st.tauc;
        sv.dthc(:,:,:,isv) = st.dthc; sv.dth(:,:,isv) = st.dthc(:,:,1);
    end
    if mod(it, 500) == 0
        fprintf('chain %d  it %5d  logL %.1f  kv=%d kx=%d  r=%s  %.0f s\n', ic, it, sum(Le), ...
            numel(st.m.V.zc), numel(st.m.X.zc), mat2str(round(st.rho*100)/100), toc);
    end
end
sv.acc = acc; sv.accn = accn;
st = rmfield(st, 'M'); st.dth = st.dthc(:,:,1); st.tau = st.tauc(:,:,1);
fprintf('chain %d done: acceptance (b/d/m/v) Vsv %s %%, xi %s %%, nuisance %.0f%%\n', ic, ...
    mat2str(round(100*squeeze(acc(1,:,2)./max(acc(1,:,1),1)))), ...
    mat2str(round(100*squeeze(acc(2,:,2)./max(acc(2,:,1),1)))), 100*accn(2)/accn(1));
end
