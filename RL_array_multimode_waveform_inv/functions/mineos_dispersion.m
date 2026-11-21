function [CL, CR, vsv, vsh, eta] = mineos_dispersion(m, cfg, which)
% [CL, CR, vsv, vsh, eta] = mineos_dispersion(m, cfg, which)
%
% Love and Rayleigh phase velocities [Nper x nm] of a model m (B-spline
% profiles of Vsv, xi and, if m.E exists, eta; see lr_common.m) from the
% linearized MINEOS dispersion cfg.lin (mineos_lin_setup.m,
% mineos_lin_disp.m): every parameter's kernel times its perturbation from
% the reference card at the nodes between cfg.zmin and cfg.zmax below the
% seafloor, with Vpv, Vph and density scaled to Vsv and Vsh by the
% reference ratios and eta at the reference unless it is a model
% parameter. which = 'both' (default), 'love' or 'rayleigh' ([] for the
% other). vsv, vsh (km/s), eta: the model at the nodes. m = [] gives the
% reference.
%
% jbrussell 11/2026
if nargin < 3, which = 'both'; end
L = cfg.lin;
if isempty(m)
    vsv = L.vsv0; vsh = L.vsh0; eta = L.eta0;
else
    vsv = bspline_profile([cfg.zmin; m.V.zc(:); cfg.zmax], [m.V.ev(1); m.V.vc(:); m.V.ev(2)], L.z);
    xi  = bspline_profile([cfg.zmin; m.X.zc(:); cfg.zmax], [m.X.ev(1); m.X.vc(:); m.X.ev(2)], L.z);
    vsh = vsv .* sqrt(xi);
    if isfield(m, 'E') && ~isempty(m.E)
        eta = bspline_profile([cfg.zmin; m.E.zc(:); cfg.zmax], [m.E.ev(1); m.E.vc(:); m.E.ev(2)], L.z);
    else
        eta = L.eta0;
    end
end
[CL, CR] = mineos_lin_disp(L, vsv, vsh, eta, which);
end
