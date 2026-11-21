function varargout = lr_common(action, varargin)
% LR_COMMON  Shared pieces of the joint Love/Rayleigh samplers (smooth
% B-spline profiles, no discontinuities, bottom tied to the reference;
% MINEOS dispersion).
%
%   m = lr_common('draw_prior', cfg, mc)            random model (start)
%   [m, ok, set, mv] = lr_common('propose', m, cfg, mc, T)
%   [S, ok] = lr_common('disp', S, m, cfg, Tdisp, fb, nm, which)
%   v = lr_common('profiles', S, zplot, cfg)        [vsv vsh xi eta] on zplot
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
% Moves (xi with prob p_xi, eta with p_eta if free, Vsv otherwise):
% 1 birth, 2 death, 3 move of an interior point, 4 value change of an
% interior point or the top end.
% Births draw from the prior and deaths pick uniformly, so every allowed
% move is accepted on the likelihood ratio alone.
switch action
    case 'draw_prior'
        [cfg, mc] = varargin{:};
        kv = randi([mc.kvmin mc.kvmax]);
        kx = randi([mc.kxmin min(mc.kxmax, kv)]);
        m.V = draw_one(cfg, kv, mc.vsmin, mc.vsmax, cfg.vbot);
        m.X = draw_one(cfg, kx, mc.ximin, mc.ximax, cfg.xbot);
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
        if exist('eta', 'var'), S.M.eta = eta; end S.CL = CL; S.CR = CR; S.m = m;
        if ok, S.Cf = {disp_to_band(CL, Tdisp, fb), disp_to_band(CR, Tdisp, fb)}; end
        varargout = {S, ok};
    case 'profiles'
        % [vsv vsh xi eta] on zplot: the B-spline profiles between zmin and
        % zmax, the reference model elsewhere (and eta at the reference
        % unless it is a model parameter)
        [S, zplot, cfg] = varargin{:};
        vsv = interp1(cfg.refp.z, cfg.refp.vsv, zplot);
        vsh = interp1(cfg.refp.z, cfg.refp.vsh, zplot);
        in = zplot >= cfg.zmin & zplot < cfg.zmax; m = S.m;
        vs = bspline_profile([cfg.zmin; m.V.zc(:); cfg.zmax], [m.V.ev(1); m.V.vc(:); m.V.ev(2)], zplot(in));
        xi = bspline_profile([cfg.zmin; m.X.zc(:); cfg.zmax], [m.X.ev(1); m.X.vc(:); m.X.ev(2)], zplot(in));
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
