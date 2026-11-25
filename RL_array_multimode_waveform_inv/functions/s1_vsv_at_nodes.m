function v = s1_vsv_at_nodes(zplot, vsv, zn, zmin, zmax)
% v = s1_vsv_at_nodes(zplot, vsv, zn, zmin, zmax)
%
% A profile (Vsv of a step-1 inversion, or any Vsv, Vsh, eta profile; on
% the plot grid zplot, km below the seafloor) at the inversion nodes zn of
% the card (lin.z), for s1_card_from_vsv.m and the linearized dispersion
% (D5_step2_check_linearization.m). Only the plot-grid points inside the model domain
% [zmin, zmax] are used (linear extrapolation to nodes just inside it): the
% plot grid also covers the crust, and interpolating across the Moho would
% give the uppermost-mantle node (at zmin, the mantle side of the
% discontinuity) a value halfway between the crust and the mantle.
in = zplot >= zmin - 1e-6 & zplot <= zmax + 1e-6;
v = interp1(zplot(in), vsv(in), zn, 'linear', 'extrap');
end
