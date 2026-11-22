function c = pa5_iso_ref_card(ct, zmoho, zlin)
% c = pa5_iso_ref_card(ct, zmoho, zlin)
%
% Isotropic version of the true model ct (read_model_card.m; PA5) without
% its lid/LVZ structure, as the starting model of
% SURF96_Rayl_1Dinv/test_kernel_inversion_vsvh_vpvh_rho_eta_PA5_twostage.m:
%   - between the Moho (zmoho, km below sea level, mantle side) and zlin km
%     below sea level, Vsv = Vsh = a straight line in depth between the
%     isotropic (Voigt) Vs of ct at the Moho and at zlin, i.e. no lid, no
%     low-velocity zone and no discontinuities in Vs (both nodes of the
%     68.5- and 166-km node pairs get the same value);
%   - there, Vpv = (Vpv/Vsv)_ct Vs, Vph = (Vph/Vsh)_ct Vs, rho = (rho/Vsv)_ct Vs
%     at every node (the ratios of ct, including their jumps at its
%     discontinuities), so Vp and rho jump where the ratios of ct jump;
%   - eta, Q and the node depths are those of ct everywhere, and ct is used
%     unchanged above the Moho and below zlin.
% With the true Vsv and Vsh in place of the line, the model is ct.
%
% Used by D1_step1_rayleigh_isotropic_inv.m (opt.s1_start_model = 'pa5_iso') as the step-1 starting model.
z = ct.z(:);
vs_iso = sqrt((2*ct.vsv.^2 + ct.vsh.^2)/3);
% mantle side of the Moho: the deeper node of the pair at zmoho
im = find(abs(z - zmoho) < 1e-6);
[~, k] = max(vs_iso(im)); i1 = im(k);
v1 = vs_iso(i1);
v2 = interp1(z(z > zmoho + 1e-6), vs_iso(z > zmoho + 1e-6), zlin);   % (no discontinuity at zlin)
in = find(z >= zmoho - 1e-6 & z <= zlin & vs_iso > 0);
in = in(~(abs(z(in) - zmoho) < 1e-6 & in ~= i1));                     % mantle side only
vs = v1 + (v2 - v1)*(z(in) - zmoho)/(zlin - zmoho);
c = ct;
c.vpv(in) = vs .* ct.vpv(in)./ct.vsv(in);
c.vph(in) = vs .* ct.vph(in)./ct.vsh(in);
c.rho(in) = vs .* ct.rho(in)./ct.vsv(in);
c.vsv(in) = vs; c.vsh(in) = vs;
end
