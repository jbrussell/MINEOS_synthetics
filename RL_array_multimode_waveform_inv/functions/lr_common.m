function varargout = lr_common(action, varargin)
% LR_COMMON  Shared pieces of the joint Love/Rayleigh samplers (smooth
% B-spline profiles, no discontinuities, bottom tied to the reference;
% MINEOS dispersion).
%
%   m = lr_common('draw_prior', cfg, mc)            random model (start)
%   [m, ok, set, mv] = lr_common('propose', m, cfg, mc, T)
%   [S, ok] = lr_common('disp', S, m, cfg, Tdisp, fb, nm, which)
%   v = lr_common('profiles', S, zplot, cfg)        [vsv vsh xi eta] on zplot
%   lp = lr_common('logprior', m, cfg, mc)          smoothness, Vsv and xi log priors
%   [m, lp] = lr_common('hyper', m, cfg, mc)        update the sampled xi prior strength (xi_ref_hier)
%
% Each of m.V (Vsv), m.X (xi) and, with mc.eta_free, m.E (eta) is one
% clamped cubic B-spline from the Moho (cfg.zmin) to cfg.zmax:
%   zc, vc : interior control points, U(zmin, zmax) x U(vmin, vmax)
%   ev     : [top; bottom] end values; top ~ U(vmin, vmax), bottom FIXED
%            to the reference model (cfg.vbot, cfg.xbot)
%            (cfg.vbot, cfg.xbot, cfg.ebot)
% Prior: uniform on the allowed numbers of interior points with the HARD
% CONSTRAINTS kx <= kv and ke <= kv (xi and eta never have more control
% points than Vsv). Without mc.eta_free, eta is the reference model's.
% With cfg.xi_fixed = true (Vsv-only inversions, step 1) xi is the
% reference model's everywhere (set mc.p_xi = 0, mc.kxmax = 0).
% Moves (xi with prob p_xi, eta with p_eta if free, Vsv otherwise):
% 1 birth, 2 death, 3 move of an interior point, 4 value change of an
% interior point or the top end.
% Births draw from the prior and deaths pick uniformly, so without the
% smoothness prior every allowed move is accepted on the likelihood ratio
% alone.
%
% Optional smoothness prior (mc.smooth_sig = [sig_1 sig_2]; 0 = off): a
% Gaussian penalty on the curvature of two profiles between the Moho and
% zmax, chosen by mc.smooth_param:
%   'vsv_vsh' (default): X1 = Vsv, X2 = Vsh = Vsv*sqrt(xi) (sig in km/s/km^2);
%                        xi is smoothed only through Vsh
%   'vsv_xi'           : X1 = Vsv (km/s/km^2), X2 = xi (1/km^2); Vsv and xi
%                        smoothed independently
%   log p_s(m) = -1/2 sum_{c = 1,2} int (Xc''(z)/sig_c)^2 dz / mc.smooth_len,
% i.e. independent N(0, sig_c^2) on the curvature every mc.smooth_len km
% (default 10 km); X'' by finite differences of the B-spline profile on a
% 1 km grid. It is a function of the profile, not of the control points,
% so it applies the same way for any number of control points; it
% multiplies the uniform priors (the joint prior is no longer uniform in
% k: it favours fewer, smoother control points). The acceptance ratio of
% the structural moves includes its change (untempered, so it cancels in
% the parallel-tempering swaps).
%
% Optional prior toward the reference Vsv (mc.vsv_ref_sig, km/s; 0 = off):
% a Gaussian penalty on the departure of Vsv from the reference model of
% the prep file (cfg.refp: the step-1 model, or the true model in an
% oracle test; fixed through the linearization updates), between the Moho
% and zmax,
%   log p_r(m) = -1/2 int ((Vsv(z) - Vsv_ref(z))/vsv_ref_sig)^2 dz / mc.vsv_ref_len,
% i.e. independent N(Vsv_ref, vsv_ref_sig^2) every mc.vsv_ref_len km
% (default 10 km), on the same 1 km grid. LARGER vsv_ref_sig = WEAKER
% prior. xi is not affected. It multiplies the uniform prior on Vsv
% (mc.vsmin..vsmax still bounds it) and enters the acceptance ratios like
% the smoothness prior.
%
% Optional prior on xi toward a fixed value (mc.xi_ref_sig; 0 = off, the
% default): log p_x(m) = -1/2 int ((xi(z) - xi_ref_mean)/xi_ref_sig)^2 dz /
% mc.xi_ref_len, with xi_ref_mean default 1 (isotropic; a fixed value, not
% the reference model's xi, so an oracle test does not pull toward the
% truth). It multiplies the uniform prior mc.ximin..ximax.
% Hierarchical version (mc.xi_ref_hier = true): the strength sigma_xi is
% not set but sampled with the model (m.xs), with a scale-invariant
% (log-uniform) hyperprior on mc.xi_ref_sig_range, so the data decide how
% strongly xi is pulled toward xi_ref_mean; the prior then includes its
% normalization, -N log(sigma_xi), N = (zmax - zmin)/xi_ref_len (the number
% of independent Gaussian terms), and mc.xi_ref_sig is not used.
%   [m, lp] = lr_common('hyper', m, cfg, mc)   Metropolis updates of m.xs
%                                              (no forward calculation)
% Hierarchical smoothing (mc.smooth_hier = true): the same for the two
% curvature strengths (m.ss = [sig_1 sig_2], log-uniform hyperprior on
% mc.smooth_sig_range, normalization -N_c log(sig_c) with N_c the number of
% curvature terms per smooth_len km); mc.smooth_sig is then not used. With
% cfg.xi_fixed only the Vsv curvature is used.
% Depth-dependent Vsv prior: mc.vsv_ref_sigz (km/s on the 1 km grid
% zmin:zmax, e.g. from the step-1 posterior std, opt.vsv_ref_sig = 'step1';
% vsv_ref_sigz.m) replaces the scalar mc.vsv_ref_sig.
switch action
    case 'draw_prior'
        [cfg, mc] = varargin{:};
        kv = randi([mc.kvmin mc.kvmax]);
        kx = randi([mc.kxmin min(mc.kxmax, kv)]);
        m.V = draw_one(cfg, kv, mc.vsmin, mc.vsmax, cfg.vbot);
        m.X = draw_one(cfg, kx, mc.ximin, mc.ximax, cfg.xbot);
        if xi_hier(mc), m.xs = draw_xs(mc); end
        if sm_hier(mc), m.ss = draw_ss(mc); end
        if eta_on(mc)
            ke = randi([mc.kemin min(mc.kemax, kv)]);
            m.E = draw_one(cfg, ke, mc.etamin, mc.etamax, cfg.ebot);
        end
        varargout = {m};
    case 'propose'
        [m, cfg, mc, T] = varargin{:};
        sc = T^0.25; ok = true;
        pe = 0; if eta_on(mc), pe = mc.p_eta; end
        u = rand;
        if u < mc.p_xi
            set = 'xi'; P = m.X; vmin = mc.ximin; vmax = mc.ximax; sv = mc.step_xi;
            kmin = mc.kxmin; kmax = min(mc.kxmax, numel(m.V.zc));
        elseif u < mc.p_xi + pe
            set = 'eta'; P = m.E; vmin = mc.etamin; vmax = mc.etamax; sv = mc.step_eta;
            kmin = mc.kemin; kmax = min(mc.kemax, numel(m.V.zc));
        else
            set = 'vsv'; P = m.V; vmin = mc.vsmin; vmax = mc.vsmax; sv = mc.step_vs;
            kmin = max(mc.kvmin, numel(m.X.zc));
            if eta_on(mc), kmin = max(kmin, numel(m.E.zc)); end
            kmax = mc.kvmax;
        end
        k = numel(P.zc); mv = randi(4);
        switch mv
            case 1  % birth (rejected if it would give kx > kv)
                if k >= kmax, ok = false; else
                    P.zc = [P.zc; cfg.zmin + (cfg.zmax - cfg.zmin)*rand];
                    P.vc = [P.vc; vmin + (vmax - vmin)*rand]; end
            case 2  % death (Vsv death rejected if it would give kv < kx)
                if k <= kmin, ok = false; else
                    i = randi(k); P.zc(i) = []; P.vc(i) = []; end
            case 3  % move
                if k == 0, ok = false; else
                    i = randi(k); P.zc(i) = P.zc(i) + mc.step_z*sc*randn;
                    ok = P.zc(i) > cfg.zmin && P.zc(i) < cfg.zmax; end
            case 4  % value of an interior point or of the (free) top end
                i = randi(k + 1);
                if i <= k, P.vc(i) = P.vc(i) + sv*sc*randn; ok = P.vc(i) > vmin && P.vc(i) < vmax;
                else, P.ev(1) = P.ev(1) + sv*sc*randn; ok = P.ev(1) > vmin && P.ev(1) < vmax; end
        end
        switch set, case 'xi', m.X = P; case 'eta', m.E = P; otherwise, m.V = P; end
        varargout = {m, ok, set, mv};
    case 'disp'
        % MINEOS dispersion (linearized about the reference card; see
        % mineos_dispersion.m). Tdisp, fb, nm are kept for compatibility.
        [S, m, cfg, Tdisp, fb, nm, which] = varargin{:};
        if strcmp(which, 'love')
            [CL, ~, vsv, vsh] = mineos_dispersion(m, cfg, 'love'); CR = S.CR;
        else
            [CL, CR, vsv, vsh, eta] = mineos_dispersion(m, cfg, 'both');
        end
        ok = all(isfinite(CL(cfg.need))) && all(isfinite(CR(cfg.need)));
        S.M = struct('vsv', vsv, 'vsh', vsh);
        if exist('eta', 'var'), S.M.eta = eta; end
        S.CL = CL; S.CR = CR; S.m = m;
        if ok, S.Cf = {disp_to_band(CL, Tdisp, fb), disp_to_band(CR, Tdisp, fb)}; end
        varargout = {S, ok};
    case 'logprior'
        [m, cfg, mc] = varargin{:};
        lp = 0; sg = [0 0]; sr = 0; sx = 0;
        if isfield(mc, 'smooth_sig') && ~isempty(mc.smooth_sig), sg = mc.smooth_sig(:)' .* [1 1]; end
        hs = sm_hier(mc);
        if hs                                       % hierarchical smoothing: sampled strengths
            if isfield(m, 'ss') && numel(m.ss) == 2, sg = m.ss(:)'; else, sg = sqrt(prod(ss_range(mc)))*[1 1]; end
        end
        if isfield(mc, 'vsv_ref_sigz') && ~isempty(mc.vsv_ref_sigz)      % depth-dependent (e.g. step-1 posterior std)
            sr = mc.vsv_ref_sigz(:);
        elseif isfield(mc, 'vsv_ref_sig') && isnumeric(mc.vsv_ref_sig) && ~isempty(mc.vsv_ref_sig)
            sr = mc.vsv_ref_sig;
        end
        if isfield(mc, 'xi_ref_sig') && ~isempty(mc.xi_ref_sig), sx = mc.xi_ref_sig; end
        hx = xi_hier(mc);
        if hx                                       % hierarchical: sampled strength
            if isfield(m, 'xs') && ~isempty(m.xs), sx = m.xs; else, sx = sqrt(prod(xs_range(mc))); end   % (state without m.xs yet)
        end
        if ~any(sg > 0) && ~any(sr > 0) && ~(sx > 0), varargout = {0}; return; end
        dz = 1; zg = (cfg.zmin:dz:cfg.zmax)';
        vs = bspline_profile([cfg.zmin; m.V.zc(:); cfg.zmax], [m.V.ev(1); m.V.vc(:); m.V.ev(2)], zg);
        xi = xi_of(m, cfg, zg);
        if sx > 0                                   % prior on xi toward a fixed value
            xm = 1; if isfield(mc, 'xi_ref_mean') && ~isempty(mc.xi_ref_mean), xm = mc.xi_ref_mean; end
            lx = 10; if isfield(mc, 'xi_ref_len') && ~isempty(mc.xi_ref_len), lx = mc.xi_ref_len; end
            lp = lp - 0.5*sum(((xi - xm)/sx).^2)*dz/lx;
            if hx, lp = lp - numel(zg)*dz/lx * log(sx); end   % normalization (sigma_xi sampled)
        end
        if any(sr > 0)                              % prior toward the reference Vsv
            vr = interp1(cfg.refp.z, cfg.refp.vsv, zg, 'linear', 'extrap');
            lr = 10; if isfield(mc, 'vsv_ref_len') && ~isempty(mc.vsv_ref_len), lr = mc.vsv_ref_len; end
            if numel(sr) > 1 && numel(sr) ~= numel(zg), error('lr_common: vsv_ref_sigz has %d values, the 1 km grid %d', numel(sr), numel(zg)); end
            lp = lp - 0.5*sum(((vs - vr)./sr).^2)*dz/lr;
        end
        if any(sg > 0)
            sp = 'vsv_vsh'; if isfield(mc, 'smooth_param') && ~isempty(mc.smooth_param), sp = mc.smooth_param; end
            switch sp
                case 'vsv_vsh', P = [vs, vs .* sqrt(xi)];
                case 'vsv_xi',  P = [vs, xi];
                otherwise, error('lr_common: unknown smooth_param ''%s'' (vsv_vsh or vsv_xi)', sp);
            end
            ell = 10; if isfield(mc, 'smooth_len') && ~isempty(mc.smooth_len), ell = mc.smooth_len; end
            for c = find(sg > 0)
                if c == 2 && isfield(cfg, 'xi_fixed') && cfg.xi_fixed && hs, continue; end   % (Vsv only: one profile)
                d2 = diff(P(:, c), 2) / dz^2;
                lp = lp - 0.5*sum((d2/sg(c)).^2)*dz/ell;
                if hs, lp = lp - numel(d2)*dz/ell * log(sg(c)); end   % normalization (sigma sampled)
            end
        end
        varargout = {lp};
    case 'hyper'
        % Metropolis updates of the sampled prior strengths (log random
        % walks; log-uniform hyperpriors, so the acceptance ratio is the ratio
        % of the priors incl. their normalization; no forward calculation):
        % m.xs (xi toward xi_ref_mean, mc.xi_ref_hier) and m.ss (curvature of
        % the two smoothed profiles, mc.smooth_hier)
        [m, cfg, mc] = varargin{:};
        hx = xi_hier(mc); hs = sm_hier(mc);
        if hx && (~isfield(m, 'xs') || isempty(m.xs)), m.xs = draw_xs(mc); end
        if hs && (~isfield(m, 'ss') || numel(m.ss) ~= 2), m.ss = draw_ss(mc); end
        lp = lr_common('logprior', m, cfg, mc);
        if ~hx && ~hs, varargout = {m, lp}; return; end
        st = 0.3; if isfield(mc, 'xi_ref_step') && ~isempty(mc.xi_ref_step), st = mc.xi_ref_step; end
        nss = 2; if isfield(cfg, 'xi_fixed') && cfg.xi_fixed, nss = 1; end
        for j = 1:5
            if hx
                mn = m; mn.xs = m.xs * exp(st*randn); rg = xs_range(mc);
                if mn.xs >= rg(1) && mn.xs <= rg(2)
                    lpn = lr_common('logprior', mn, cfg, mc);
                    if log(rand) < lpn - lp, m = mn; lp = lpn; end
                end
            end
            if hs
                rg = ss_range(mc);
                for c = 1:nss
                    mn = m; mn.ss(c) = m.ss(c) * exp(st*randn);
                    if mn.ss(c) < rg(1) || mn.ss(c) > rg(2), continue; end
                    lpn = lr_common('logprior', mn, cfg, mc);
                    if log(rand) < lpn - lp, m = mn; lp = lpn; end
                end
            end
        end
        varargout = {m, lp};
    case 'profiles'
        % [vsv vsh xi eta] on zplot: the B-spline profiles between zmin and
        % zmax, the reference model elsewhere (and eta at the reference
        % unless it is a model parameter)
        [S, zplot, cfg] = varargin{:};
        vsv = interp1(cfg.refp.z, cfg.refp.vsv, zplot);
        vsh = interp1(cfg.refp.z, cfg.refp.vsh, zplot);
        in = zplot >= cfg.zmin & zplot < cfg.zmax; m = S.m;
        vs = bspline_profile([cfg.zmin; m.V.zc(:); cfg.zmax], [m.V.ev(1); m.V.vc(:); m.V.ev(2)], zplot(in));
        xi = xi_of(m, cfg, zplot(in));
        vsv(in) = vs; vsh(in) = vs .* sqrt(xi);
        if isfield(cfg.refp, 'eta'), eta = interp1(cfg.refp.z, cfg.refp.eta, zplot); else, eta = ones(size(zplot)); end
        if isfield(m, 'E') && ~isempty(m.E)
            eta(in) = bspline_profile([cfg.zmin; m.E.zc(:); cfg.zmax], [m.E.ev(1); m.E.vc(:); m.E.ev(2)], zplot(in));
        end
        varargout = {[vsv, vsh, (vsh./vsv).^2, eta]};
end
end

function on = eta_on(mc)
on = isfield(mc, 'eta_free') && mc.eta_free;
end

function P = draw_one(cfg, k, vmin, vmax, vbot)
P.zc = cfg.zmin + (cfg.zmax - cfg.zmin)*rand(k, 1);
P.vc = vmin + (vmax - vmin)*rand(k, 1);
P.ev = [vmin + (vmax - vmin)*rand; vbot];          % bottom tied to reference
end

function h = xi_hier(mc)
% hierarchical xi prior (strength sampled) on?
h = isfield(mc, 'xi_ref_hier') && ~isempty(mc.xi_ref_hier) && mc.xi_ref_hier;
end

function rg = xs_range(mc)
% range of the log-uniform hyperprior on the xi prior strength
rg = [0.005 0.2];
if isfield(mc, 'xi_ref_sig_range') && numel(mc.xi_ref_sig_range) == 2, rg = mc.xi_ref_sig_range; end
end

function xs = draw_xs(mc)
% draw from the log-uniform hyperprior
rg = xs_range(mc); xs = exp(log(rg(1)) + (log(rg(2)) - log(rg(1)))*rand);
end

function xi = xi_of(m, cfg, z)
% xi profile at depths z: the model's B-spline, or the reference model's
% xi when it is held fixed (cfg.xi_fixed, Vsv-only inversions)
if isfield(cfg, 'xi_fixed') && cfg.xi_fixed
    xi = interp1(cfg.refp.z, (cfg.refp.vsh ./ cfg.refp.vsv).^2, z, 'linear', 'extrap');
else
    xi = bspline_profile([cfg.zmin; m.X.zc(:); cfg.zmax], [m.X.ev(1); m.X.vc(:); m.X.ev(2)], z);
end
end

function h = sm_hier(mc)
% hierarchical smoothing (curvature prior strengths sampled) on?
h = isfield(mc, 'smooth_hier') && ~isempty(mc.smooth_hier) && mc.smooth_hier;
end

function rg = ss_range(mc)
% range of the log-uniform hyperprior on the curvature prior strengths
rg = [5e-5 5e-3];
if isfield(mc, 'smooth_sig_range') && numel(mc.smooth_sig_range) == 2, rg = mc.smooth_sig_range; end
end

function ss = draw_ss(mc)
% draw both strengths from the log-uniform hyperprior
rg = ss_range(mc); ss = exp(log(rg(1)) + (log(rg(2)) - log(rg(1)))*rand(1, 2));
end
