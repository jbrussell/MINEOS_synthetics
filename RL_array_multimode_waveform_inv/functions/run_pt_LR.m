function [sv, st, info, S] = run_pt_LR(rec, cfg, mc, Tdisp, fb, nm, zplot, seed, lin0, pt, ie, S0)
% RUN_PT_LR  Parallel-tempering rjMcMC for the joint Love/Rayleigh
% inversion (independent B-spline profiles for Vsv and xi, kx <= kv, no
% discontinuities, bottom tied to the reference). With a single
% temperature (pt.T = 1) it is an ordinary rjMcMC chain.
%
% Nuisance parameters (Laplace marginalization): for every record and
% branch, the smooth spectral amplitudes B_n(w) (exactly linear), the
% frequency-dependent time shift tau_n(w) and the frequency-dependent
% arrival angle dth_n(w) are integrated out analytically, tau and dth by
% linearizing the waveforms about a linearization point (nuis_relin.m).
% Every state of every temperature carries its own linearization points,
% which are moved to the posterior mean of tau, dth and B for the state's
% current structure every mc.relin iterations (default 1: before every
% proposal) by one Gauss-Newton step, and the state's likelihood is then
% re-evaluated there. The posterior mean used for the update is cached
% from the likelihood evaluation that produced the state's current
% log-likelihood (no extra solve). Only Vsv, xi and the 8 noise parameters (sigma and
% AR(1) correlation r per window type) are sampled and tempered.
%
% lin0 : {nrec} initial linearization points (e.g. fit_nuisance_LR at the
%        reference dispersion)
% pt   : T (temperatures, T(1) = 1), nit, burn, thin; optional jump_last:
%        cycle-skip checks (cycle_jump.m, every mc.jump_every iterations,
%        all temperatures) up to this iteration (default pt.burn)
% pt.ckpt, pt.ckpt_every : OPTIONAL checkpoint file and interval
%        (iterations, default 100): see the checkpoint section below
% S0   : OPTIONAL {K} initial states (fields m, sig, rho and, optionally,
%        lin: {nrec} linearization points), e.g. the final states S of an
%        earlier call (continuation after a reference update); default:
%        independent prior draws with linearization points lin0
% Only the T = 1 chain is saved (sv); S returns the final states of all
% temperatures (without the stored operators).

rng(seed);
nrec = numel(rec); wave = [rec.wave]'; K = numel(pt.T);
Kb = size(cfg.Phi, 2); Ka = size(cfg.Psi, 2); Kt = size(cfg.Psit, 2);
relin = max(1, round(mc.relin));
jev = 0; if isfield(mc, 'jump_every') && ~isempty(mc.jump_every), jev = round(mc.jump_every); end
jlast = pt.burn; if isfield(pt, 'jump_last'), jlast = pt.jump_last; end
info.jumps = zeros(0, 3);                         % [iteration, temperature index, accepted changes]
lay = @(lr) spec_layout(lr, nm, Kb, Ka, mc.sb, mc.taup, mc.dthp);

% ---------- independent prior draws for every temperature ----------
if nargin < 12, S0 = {}; end
S = cell(K, 1);
for k = 1:K
    if ~isempty(S0)
        [s, ok] = lr_common('disp', struct(), S0{k}.m, cfg, Tdisp, fb, nm, 'both');
        if ~ok, error('run_pt_LR: initial state %d has no valid dispersion', k); end
        s.sig = S0{k}.sig; s.rho = S0{k}.rho;
        if isfield(S0{k}, 'lin'), s.lin = S0{k}.lin; else, s.lin = lin0; end
    else
        while true
            m = lr_common('draw_prior', cfg, mc);
            [s, ok] = lr_common('disp', struct(), m, cfg, Tdisp, fb, nm, 'both');
            if ok, break; end
        end
        s.sig = zeros(1, 4);
        for wv = 1:2, r1_ = find(wave == wv, 1);
            if isempty(r1_), s.sig(2*wv-1:2*wv) = 1; continue; end          % (wave type not used, opt.waves)
            nl_ = rec_nlev(rec(r1_)); if isempty(nl_), nl_ = 1; end
            s.sig(2*wv-1:2*wv) = std(rec(r1_).d ./ nl_); end
        s.rho = 0.5*ones(1, 4);
        s.lin = lin0;
    end
    s.G = cell(nrec, 1); s.Le = zeros(nrec, 1); s.xh = cell(nrec, 1);
    for r = 1:nrec
        [~, s.Le(r), s.G{r}, s.xh{r}] = nuis_relin(rec(r), s.Cf{wave(r)}, s.lin{r}, cfg, s.sig, s.rho, mc, nm);
    end
    S{k} = s;
end
st0 = rmfield(S{1}, {'M', 'G', 'Le', 'lin', 'xh'}); st0.prof = lr_common('profiles', S{1}, zplot, cfg);

nsave = floor((pt.nit - pt.burn)/pt.thin);
sv.st0 = st0;
sv.vsv = zeros(numel(zplot), nsave); sv.vsh = sv.vsv; sv.xi = sv.vsv; sv.eta = sv.vsv;
sv.CL = zeros(numel(Tdisp), nm, nsave); sv.CR = sv.CL;
sv.kv = zeros(nsave, 1); sv.kx = sv.kv; sv.sig = zeros(nsave, 4); sv.rho = sv.sig; sv.xs = NaN(nsave, 1); sv.ss = NaN(nsave, 2);
sv.tau = zeros(nrec, nm, nsave); sv.dth = sv.tau;          % band-mean (constant) terms
sv.tauc = zeros(nrec, nm, Kt, nsave);                      % all time-shift coefficients
sv.dthc = zeros(nrec, nm, Ka, nsave);                      % all arrival-angle coefficients
sv.L = zeros(pt.nit, 1);
info.T = pt.T; info.swap = zeros(max(K-1, 1), 2); info.Ltrace = zeros(pt.nit, K);
info.acc = zeros(3, 4, 2);                % [vsv; xi; eta] x [b d m v] x [tried acc]
isv = 0; tic;
fprintf('PT %d: %d temperature(s), T = %s, relinearization every %d iteration(s)\n', ...
    ie, K, mat2str(round(pt.T*100)/100), relin);

% ---------- checkpoints (pt.ckpt, every pt.ckpt_every iterations) ----------
% the states of all temperatures (without the stored operators, recomputed
% here), the samples so far, info, the iteration and the random-number
% state; a run started with an existing checkpoint resumes from it
ckf = ''; ckev = 100; it0 = 1;
if isfield(pt, 'ckpt') && ~isempty(pt.ckpt), ckf = pt.ckpt; end
if isfield(pt, 'ckpt_every') && ~isempty(pt.ckpt_every), ckev = pt.ckpt_every; end
if ~isempty(ckf) && exist(ckf, 'file') == 2
    C = load(ckf);
    S = C.S; sv = C.sv; info = C.info; isv = C.isv; it0 = C.it + 1;
    for k = 1:K
        s = S{k}; s.G = cell(nrec, 1); s.Le = zeros(nrec, 1); s.xh = cell(nrec, 1);
        for r = 1:nrec
            [~, s.Le(r), s.G{r}, s.xh{r}] = nuis_relin(rec(r), s.Cf{wave(r)}, s.lin{r}, cfg, s.sig, s.rho, mc, nm);
        end
        S{k} = s;
    end
    rng(C.rngst);
    fprintf('PT %d: resumed from %s at iteration %d of %d\n', ie, ckf, C.it, pt.nit);
    clear C
end
for it = it0:pt.nit
    for k = 1:K
        s = S{k}; Tk = pt.T(k);
        % ---------- relinearize tau, dth and B about the current structure ----------
        % s.xh{r} is the posterior mean of [beta; dtau; dth] cached by the
        % likelihood evaluation that produced s.Le(r) (same G, sig, rho, lin)
        if mod(it - 1, relin) == 0
            for r = 1:nrec
                [~, ~, itau, ith] = lay(s.lin{r});
                xh = s.xh{r}; lr = s.lin{r};
                lr.b0   = reshape(xh(1:2:2*nm*Kb) + 1i*xh(2:2:2*nm*Kb), Kb, nm).';
                lr.tau0 = lr.tau0 + max(min(reshape(xh(itau), size(itau)), mc.max_dtau), -mc.max_dtau);
                lr.dth0 = lr.dth0 + max(min(reshape(xh(ith), size(ith)), mc.max_dth), -mc.max_dth);
                s.lin{r} = lr;
                [~, s.Le(r), s.G{r}, s.xh{r}] = nuis_relin(rec(r), s.Cf{wave(r)}, lr, cfg, s.sig, s.rho, mc, nm);
            end
        end
        % ---------- cycle-skip check of the time shifts (burn-in) ----------
        if jev > 0 && it <= jlast && mod(it, jev) == 0
            nj = 0;
            for r = 1:nrec
                [s.lin{r}, s.Le(r), s.G{r}, s.xh{r}, a] = cycle_jump(rec(r), s.Cf{wave(r)}, s.lin{r}, ...
                    s.Le(r), s.G{r}, s.xh{r}, cfg, s.sig, s.rho, mc, nm);
                nj = nj + a;
            end
            info.jumps(end+1, :) = [it, k, nj];
            if k == 1 && nj > 0
                fprintf('PT %d  it %5d  cycle-skip check: %d time-shift change(s) kept (T = 1)\n', ie, it, nj);
            end
        end
        % ---------- structural move ----------
        [mn, okp, set, mv] = lr_common('propose', s.m, cfg, mc, Tk);
        is = find(strcmp(set, {'vsv', 'xi', 'eta'}));
        if okp
            if k == 1, info.acc(is, mv, 1) = info.acc(is, mv, 1) + 1; end
            % xi changes Vsh, which both Love and Rayleigh waves see: recompute both
            which = 'both'; rr = 1:nrec;
            [Sn, ok] = lr_common('disp', s, mn, cfg, Tdisp, fb, nm, which);
            if ok
                Gn = s.G; Ln = s.Le; Xn = s.xh;
                for r = rr
                    [~, Ln(r), Gn{r}, Xn{r}] = nuis_relin(rec(r), Sn.Cf{wave(r)}, s.lin{r}, cfg, s.sig, s.rho, mc, nm);
                end
                dlp = lr_common('logprior', mn, cfg, mc) - lr_common('logprior', s.m, cfg, mc);
                if log(rand) < (sum(Ln) - sum(s.Le)) / Tk + dlp      % smoothness prior untempered
                    Sn.G = Gn; Sn.Le = Ln; Sn.xh = Xn; Sn.lin = s.lin; s = Sn;
                    if k == 1, info.acc(is, mv, 2) = info.acc(is, mv, 2) + 1; end
                end
            end
        end
        % ---------- sampled prior strengths (opt.xi_ref_hier, opt.smooth_hier; prior only, untempered) ----------
        if ((isfield(mc, 'xi_ref_hier') && mc.xi_ref_hier) || (isfield(mc, 'smooth_hier') && mc.smooth_hier)), s.m = lr_common('hyper', s.m, cfg, mc); end
        % ---------- noise move: level sigma or AR(1) correlation r ----------
        sg = s.sig; rh = s.rho; iw = randi(4); okn = true;
        if rand < 0.5 || (isfield(mc, 'fixed_corr') && mc.fixed_corr), sg(iw) = sg(iw)*exp(mc.step_lsig*randn);   % (fixed correlation: sigma only)
        else, rh(iw) = rh(iw) + mc.step_rho*randn; okn = rh(iw) >= 0 && rh(iw) < mc.rhomax; end
        if isfield(mc, 'single_window') && mc.single_window   % one window per wave type: tied slots
            pr_ = 2*ceil(iw/2) - [1 0]; sg(pr_) = sg(iw); rh(pr_) = rh(iw);
        end
        if okn
            Ln = s.Le; Xn = s.xh;
            for r = find(wave == ceil(iw/2))'
                [mu, sb2] = lay(s.lin{r});
                [Ln(r), Xn{r}] = event_loglike_ar(rec(r).d, s.G{r}, rec(r).wt, rec(r).first, sg, rh, mu, sb2, rec_nlev(rec(r)), rec_W(rec(r)), rec_wld(rec(r)));
            end
            if log(rand) < (sum(Ln) - sum(s.Le)) / Tk, s.sig = sg; s.rho = rh; s.Le = Ln; s.xh = Xn; end
        end
        S{k} = s; info.Ltrace(it, k) = sum(s.Le);
    end
    % ---------- swaps between neighbouring temperatures (alternating parity) ----------
    for p = (1 + mod(it, 2)):2:(K-1)
        info.swap(p, 1) = info.swap(p, 1) + 1;
        if log(rand) < (sum(S{p+1}.Le) - sum(S{p}.Le)) * (1/pt.T(p) - 1/pt.T(p+1))
            tmp = S{p}; S{p} = S{p+1}; S{p+1} = tmp; info.swap(p, 2) = info.swap(p, 2) + 1;
        end
    end
    % ---------- save the T = 1 state ----------
    s = S{1}; sv.L(it) = sum(s.Le);
    if it > pt.burn && mod(it - pt.burn, pt.thin) == 0
        isv = isv + 1;
        pr = lr_common('profiles', s, zplot, cfg);
        sv.vsv(:, isv) = pr(:,1); sv.vsh(:, isv) = pr(:,2); sv.xi(:, isv) = pr(:,3); sv.eta(:, isv) = pr(:,4);
        sv.CL(:,:,isv) = s.CL; sv.CR(:,:,isv) = s.CR;
        sv.kv(isv) = numel(s.m.V.zc); sv.kx(isv) = numel(s.m.X.zc); sv.sig(isv,:) = s.sig; sv.rho(isv,:) = s.rho;
        if isfield(s.m, 'xs'), sv.xs(isv) = s.m.xs; end
        if isfield(s.m, 'ss'), sv.ss(isv, :) = s.m.ss; end
        for r = 1:nrec
            [~, ~, itau, ith] = lay(s.lin{r});
            xh = s.xh{r};                                                   % cached posterior mean
            tc = s.lin{r}.tau0 + reshape(xh(itau), size(itau)); dc = s.lin{r}.dth0 + reshape(xh(ith), size(ith));   % nm x Kt, nm x Ka
            sv.tauc(r,:,:,isv) = reshape(tc, 1, nm, Kt); sv.tau(r,:,isv) = tc(:, 1)';
            sv.dthc(r,:,:,isv) = reshape(dc, 1, nm, Ka); sv.dth(r,:,isv) = dc(:, 1)';
        end
    end
    if mod(it, 100) == 0
        fprintf('PT %d  it %5d  logL(T=1) %.1f  kv=%d kx=%d  r=%s  swap%% %s  %.0f s\n', ie, it, ...
            sum(S{1}.Le), numel(S{1}.m.V.zc), numel(S{1}.m.X.zc), mat2str(round(S{1}.rho*100)/100), ...
            mat2str(round(100*info.swap(:,2)./max(info.swap(:,1),1))'), toc);
    end
    if ~isempty(ckf) && mod(it, ckev) == 0 && it < pt.nit
        ckpt_save(ckf, cellfun(@(x) rmfield(x, {'G', 'xh'}), S, 'uniformoutput', false), sv, info, isv, it, rng);
    end
end
info.Lmean = mean(info.Ltrace(min(pt.burn+1, pt.nit):end, :), 1)';
sv.acc = info.acc;
st = rmfield(S{1}, {'M', 'G', 'Le', 'xh'});
i0 = max(1, isv-20);
st.tauc = mean(sv.tauc(:,:,:,i0:isv), 4); st.dthc = mean(sv.dthc(:,:,:,i0:isv), 4);
st.tau = st.tauc(:,:,1); st.dth = st.dthc(:,:,1);
fprintf('PT %d done: acceptance (b/d/m/v) Vsv %s %%, xi %s %%, swap acceptance %s %%\n', ie, ...
    mat2str(round(100*squeeze(info.acc(1,:,2)./max(info.acc(1,:,1),1)))), ...
    mat2str(round(100*squeeze(info.acc(2,:,2)./max(info.acc(2,:,1),1)))), ...
    mat2str(round(100*info.swap(:,2)./max(info.swap(:,1),1))'));
if any(info.acc(3,:,1))
    fprintf('PT %d: eta acceptance (b/d/m/v) %s %%\n', ie, mat2str(round(100*squeeze(info.acc(3,:,2)./max(info.acc(3,:,1),1)))));
end
S = cellfun(@(x) rmfield(x, {'G', 'xh', 'M'}), S, 'uniformoutput', false);
end
