function c = ref_card_update(c, zw, zmin, zmax, zp, vsv, vsh, eta)
% c = ref_card_update(c, zw, zmin, zmax, zp, vsv, vsh, eta)
%
% New MINEOS card (read_model_card.m structure, m/s and kg/m^3) from the
% card c with Vsv and Vsh replaced by the profiles vsv(zp), vsh(zp) (km/s;
% zp in km below the seafloor) at the nodes between zmin and zmax below the
% seafloor (seafloor at zw km below sea level; at zmin only the mantle-side
% node). Vpv, Vph and density follow with the Vpv/Vsv, Vph/Vsh and rho/Vsv
% ratios of c, exactly as in the linearized forward model
% (mineos_lin_setup.m). eta is replaced by eta(zp) if given (and not
% empty), otherwise kept; everything outside [zmin, zmax] is kept.
zsf = c.z - zw;
in = find(zsf >= zmin - 1e-6 & zsf <= zmax + 1e-6 & c.vsv > 0);
atm = in(abs(zsf(in) - zmin) < 1e-6);            % node pair at the Moho:
if numel(atm) > 1, in = setdiff(in, atm(2:end)); end   % keep the mantle side (card order: deeper first)
ok = isfinite(vsv(:)) & isfinite(vsh(:)) & zp(:) >= zmin - 1e-6;
rpv = c.vpv(in) ./ c.vsv(in); rph = c.vph(in) ./ c.vsh(in); rrho = c.rho(in) ./ c.vsv(in);
c.vsv(in) = 1000*interp1(zp(ok), vsv(ok), zsf(in), 'linear', 'extrap');
c.vsh(in) = 1000*interp1(zp(ok), vsh(ok), zsf(in), 'linear', 'extrap');
c.vpv(in) = rpv .* c.vsv(in);
c.vph(in) = rph .* c.vsh(in);
c.rho(in) = rrho .* c.vsv(in);
if nargin >= 8 && ~isempty(eta)
    c.eta(in) = interp1(zp(ok), eta(ok), zsf(in), 'linear', 'extrap');
end
end
