function c = smooth_lid_lvz_card(ct, zmoho, sig, zblend)
% c = smooth_lid_lvz_card(ct, zmoho, sig, zblend)
%
% Starting model for step 1 built from the true model ct (read_model_card.m):
% an ISOTROPIC Vs profile (Voigt average of ct) whose upper-mantle structure
% (lid, low-velocity zone and the discontinuities between them) is smoothed
% away, with the Vp/Vs, rho/Vs and eta of the true model.
%   - mantle below the Moho (zmoho, km below sea level): the isotropic Vs of
%     ct is averaged with a Gaussian of standard deviation sig (km),
%     normalized over the mantle nodes above the 410 (so the Moho and the
%     410 do not leak into it); the result is continuous across the lid/LVZ
%     discontinuities of ct (one value at both nodes of a node pair);
%   - from zblend(1) to zblend(2) km it tapers (cosine) into the isotropic
%     Vs of ct, which is used unchanged below zblend(2) and above the Moho;
%   - Vpv = (Vpv/Vsv)_true Vs, Vph = (Vph/Vsh)_true Vs, rho = (rho/Vsv)_true Vs,
%     eta = eta_true at every solid node (the true ratios, including their
%     jumps at the true discontinuities), so these match the true model when
%     Vs does; Q and the node depths are those of ct.
% The water layer and the crust are those of ct.
%
% Used by D1_step1_rayleigh_isotropic_inv.m with opt.s1_start_model = 'true_smooth'.
z = ct.z(:);                                   % depth below sea level (km)
vs_iso = sqrt((2*ct.vsv.^2 + ct.vsh.^2)/3);
solid = ct.vsv > 0;
rpv = ct.vpv./max(ct.vsv, eps); rph = ct.vph./max(ct.vsh, eps); rrho = ct.rho./max(ct.vsv, eps);
% mantle nodes between the Moho (mantle side) and the 410 (upper side)
z410 = min(z(z > 300 & [diff(z) == 0; false] | z > 300 & [false; diff(z) == 0]));
if isempty(z410), z410 = 2*zblend(2); end
im = find(z > zmoho - 1e-6 & z < z410 + 1e-6 & solid);
im = im(~(abs(z(im) - zmoho) < 1e-6 & vs_iso(im) < 4000));   % mantle side of the Moho pair
im = im(~(abs(z(im) - z410) < 1e-6 & vs_iso(im) > min(vs_iso(abs(z - z410) < 1e-6))));
% one value per depth (average of a node pair), Gaussian average in depth
[zu, ~, ju] = unique(z(im));
vu = accumarray(ju, vs_iso(im), [], @mean);
du = diff([zu(1); (zu(1:end-1) + zu(2:end))/2; zu(end)]);   % depth weight of each node
vsm = zeros(size(zu));
for k = 1:numel(zu)
    w = exp(-0.5*((zu - zu(k))/sig).^2) .* du;
    vsm(k) = sum(w.*vu)/sum(w);
end
w = ones(size(zu));
tb = zu > zblend(1) & zu < zblend(2);
w(tb) = 0.5*(1 + cos(pi*(zu(tb) - zblend(1))/(zblend(2) - zblend(1))));
w(zu >= zblend(2)) = 0;
vsn = w.*vsm + (1 - w).*vu;
vs = vs_iso;
vs(im) = vsn(ju);
c = ct;
c.vsv = vs; c.vsh = vs;
c.vpv(solid) = rpv(solid).*vs(solid);
c.vph(solid) = rph(solid).*vs(solid);
c.rho(solid) = rrho(solid).*vs(solid);
c.eta = ct.eta;
end
