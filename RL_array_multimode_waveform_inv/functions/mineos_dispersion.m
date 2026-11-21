function [CL, CR, vsv, vsh] = mineos_dispersion(m, cfg, which)
% [CL, CR, vsv, vsh] = mineos_dispersion(m, cfg, which)
%
% Love and Rayleigh phase velocities [Nper x nm] of a model m (two
% B-spline profiles, see lr_common.m) from the linearized MINEOS
% dispersion cfg.lin (mineos_lin_setup.m): c = c_ref + G * [dVsv; dVsh] at
% the card nodes between cfg.zmin and cfg.zmax below the seafloor.
% Rayleigh waves see both Vsv and Vsh (radially anisotropic MINEOS
% kernels). which = 'both' (default), 'love' or 'rayleigh' ([] for the other).
% vsv, vsh: the model at the nodes (km/s). m = [] gives the reference.
%
% jbrussell 11/2026
if nargin < 3, which = 'both'; end
L = cfg.lin;
if isempty(m)
    vsv = L.vsv0; vsh = L.vsh0;
else
    vsv = bspline_profile([cfg.zmin; m.V.zc(:); cfg.zmax], [m.V.ev(1); m.V.vc(:); m.V.ev(2)], L.z);
    xi  = bspline_profile([cfg.zmin; m.X.zc(:); cfg.zmax], [m.X.ev(1); m.X.vc(:); m.X.ev(2)], L.z);
    vsh = vsv .* sqrt(xi);
end
dm = [vsv - L.vsv0; vsh - L.vsh0];
nm = numel(L.GL);
CL = []; CR = [];
if ~strcmp(which, 'rayleigh')
    CL = L.cL;
    for n = 1:nm, CL(:, n) = CL(:, n) + L.GL{n} * dm; end
end
if ~strcmp(which, 'love')
    CR = L.cR;
    for n = 1:nm, CR(:, n) = CR(:, n) + L.GR{n} * dm; end
end
end
