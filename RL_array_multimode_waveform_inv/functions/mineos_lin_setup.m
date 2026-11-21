function lin = mineos_lin_setup(card, KT, KS, cT, cS, zw, zmin, zmax)
% lin = mineos_lin_setup(card, KT, KS, cT, cS, zw, zmin, zmax)
%
% Linearized MINEOS dispersion about a reference card: for the card nodes
% between zmin and zmax km below the seafloor (seafloor at zw km below
% sea level), the Love (toroidal) and Rayleigh (spheroidal) phase
% velocities of branches n = 0..nm-1 at the dispersion nodes are
%
%   c_n(T) = c_ref,n(T) + sum_nodes sum_X  G^X_n(T, node) dX(node),
%   X = Vsv, Vsh, Vpv, Vph, rho, eta,
%
% with each parameter's own MINEOS kernel and its own perturbation
% dX = X - X_ref (mineos_lin_disp.m). The kernels are
%   G^X = K_X .* dr * 1000   for the velocities (km/s per km/s) and density
%                            (km/s per g/cm^3),
%   G^eta = K_eta .* dr      (km/s per unit eta),
% (K: MINEOS kernels from mineos_modes.m, dr: node weights in m; MINEOS
% units dc (km/s) = sum K .* dr (m) .* dX (m/s, kg/m^3; eta dimensionless)).
% Love kernels of Vpv, Vph and eta are zero.
%
% The perturbations of Vpv, Vph and density follow those of Vsv and Vsh
% through the reference ratios (Vpv/Vsv, Vph/Vsh and rho/Vsv held fixed;
% lin.rpv, lin.rph, lin.rrho), and eta is the reference value unless it is
% a model parameter (mc.eta_free); see mineos_lin_disp.m.
%
% card   - reference card (read_model_card.m)
% KT, KS - {nm} kernels of branches 0..nm-1 (mineos_modes.m out.K)
% cT, cS - [Nper x nm] reference phase velocities (km/s)
% zw     - water depth (km); zmin, zmax - depths below the seafloor (km)
%
% lin.z     : depths below the seafloor of the model nodes [Nn x 1]
% lin.vsv0, vsh0, vpv0, vph0 (km/s), rho0 (g/cm^3), eta0 : reference values
% lin.KL{n}, lin.KR{n} : structs of [Nper x Nn] kernels G^X (fields vsv,
%                        vsh, vpv, vph, rho, eta), Love and Rayleigh
% lin.cL, lin.cR : reference phase velocities [Nper x nm]
% lin.rpv, lin.rph, lin.rrho : the ratios used for the Vp and rho scaling
%
% Nodes: the card nodes with zmin <= z - zw <= zmax, excluding the
% crust-side node at the Moho (the model starts at the mantle side).
%
% jbrussell 11/2026

nnode = numel(card.rad);
zsf = card.z - zw;                                   % depth below the seafloor
dr = gradient(card.rad);
in = find(zsf >= zmin - 1e-6 & zsf <= zmax + 1e-6);
% at the Moho keep only the mantle-side node (the deeper of the pair)
atmoho = in(abs(zsf(in) - zmin) < 1e-6);
if numel(atmoho) > 1
    in = setdiff(in, atmoho(2:end));                 % card order: deeper node first
end
lin.idx = in(:);
lin.z = zsf(in);
lin.vsv0 = card.vsv(in)/1000;
lin.vsh0 = card.vsh(in)/1000;
lin.vpv0 = card.vpv(in)/1000;
lin.vph0 = card.vph(in)/1000;
lin.rho0 = card.rho(in)/1000;                        % g/cm^3
lin.eta0 = card.eta(in);
lin.rpv = card.vpv(in) ./ card.vsv(in);
lin.rph = card.vph(in) ./ card.vsh(in);
lin.rrho = card.rho(in) ./ card.vsv(in);              % (g/cm^3) per (km/s)
nm = size(cS, 2);
for n = 1:nm
    lin.KL{n} = node_K(KT{n}, in, nnode, dr);
    lin.KR{n} = node_K(KS{n}, in, nnode, dr);
end
lin.cL = cT; lin.cR = cS;
end

function G = node_K(K, in, nnode, dr)
% per-parameter kernels [Nper x Nn] at the nodes in (NaN -> 0)
iK = in - (nnode - numel(K.rad));                    % toroidal kernels: top nodes only
f = {'vsv', 'vsh', 'vpv', 'vph', 'rho', 'eta'};
sc = [1000 1000 1000 1000 1000 1];                   % card units -> km/s, g/cm^3; eta dimensionless
for i = 1:numel(f)
    if isfield(K, f{i})
        g = (K.(f{i})(iK, :) .* dr(in) * sc(i))';
    else
        g = zeros(size(K.vsv, 2), numel(in));
    end
    g(isnan(g)) = 0;
    G.(f{i}) = g;
end
end
