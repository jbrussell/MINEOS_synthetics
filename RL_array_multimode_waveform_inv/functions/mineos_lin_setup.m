function lin = mineos_lin_setup(card, KT, KS, cT, cS, zw, zmin, zmax)
% lin = mineos_lin_setup(card, KT, KS, cT, cS, zw, zmin, zmax)
%
% Linearized MINEOS dispersion about a reference card: for the card nodes
% between zmin and zmax km below the seafloor (seafloor at zw km below
% sea level), the Love (toroidal) and Rayleigh (spheroidal) phase
% velocities of branches n = 0..nm-1 at the dispersion nodes are
%
%   c_n(T) = c_ref,n(T) + sum_nodes [ Gsv_n(T, node) dVsv + Gsh_n(T, node) dVsh ]
%
% Vpv, Vph and density follow Vsv and Vsh with the ratios of the reference
% card (Vpv/Vsv, Vph/Vsh and rho/Vsv held fixed), so the effective kernels
% are
%   Gsv = Gvsv + (Vpv/Vsv)_0 Gvpv + (rho/Vsv)_0 Grho
%   Gsh = Gvsh + (Vph/Vsh)_0 Gvph
% with G = K .* dr * 1000 (K: MINEOS kernels from mineos_modes.m, dr: node
% weights in m; velocities in km/s, density in g/cm^3). eta is held at the
% reference. Used by mineos_dispersion.m.
%
% card   - reference card (read_model_card.m)
% KT, KS - {nm} kernels of branches 0..nm-1 (mineos_modes.m out.K)
% cT, cS - [Nper x nm] reference phase velocities (km/s)
% zw     - water depth (km); zmin, zmax - depths below the seafloor (km)
%
% lin.z     : depths below the seafloor of the model nodes [Nn x 1]
% lin.vsv0, lin.vsh0 : reference values at the nodes (km/s)
% lin.GL{n}, lin.GR{n} : [Nper x 2Nn] = [Gsv Gsh] (Love, Rayleigh)
% lin.cL, lin.cR : reference phase velocities [Nper x nm]
% lin.rpv, lin.rph, lin.rrho : the ratios used for the scaling
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
lin.rpv = card.vpv(in) ./ card.vsv(in);
lin.rph = card.vph(in) ./ card.vsh(in);
lin.rrho = card.rho(in) ./ card.vsv(in);              % (g/cm^3) per (km/s)
nm = size(cS, 2);
for n = 1:nm
    lin.GL{n} = node_G(KT{n}, in, nnode, dr, lin);
    lin.GR{n} = node_G(KS{n}, in, nnode, dr, lin);
end
lin.cL = cT; lin.cR = cS;
end

function G = node_G(K, in, nnode, dr, lin)
% [Nper x 2Nn] = [dc/dVsv dc/dVsh] (km/s per km/s) at the nodes in, with
% Vp and density scaled to Vs
iK = in - (nnode - numel(K.rad));                    % toroidal kernels: top nodes only
g = @(f) (K.(f)(iK, :) .* dr(in) * 1000)';
Gsv = g('vsv') + g('vpv') .* lin.rpv' + g('rho') .* lin.rrho';
Gsh = g('vsh') + g('vph') .* lin.rph';
Gsv(isnan(Gsv)) = 0; Gsh(isnan(Gsh)) = 0;
G = [Gsv, Gsh];
end
