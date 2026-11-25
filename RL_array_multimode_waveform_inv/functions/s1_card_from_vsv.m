function card = s1_card_from_vsv(card_ref, idx, vsv)
% card = s1_card_from_vsv(card_ref, idx, vsv)
%
% Card of a Vsv-only (step-1) model: Vsv (km/s) at the card nodes idx (the
% inversion nodes, lin.idx of mineos_lin_setup.m) replaced by vsv, and
% Vsh, Vpv, Vph and density scaled with it so that the reference card's
% xi = (Vsh/Vsv)^2, Vpv/Vsv, Vph/Vsh and rho/Vsv are kept exactly; eta, Q
% and every other node unchanged.
card = card_ref; i = idx(:); v = 1000*vsv(:);         % card units (m/s)
card.vsv(i) = v;
card.vsh(i) = v .* card_ref.vsh(i) ./ card_ref.vsv(i);
card.vpv(i) = v .* card_ref.vpv(i) ./ card_ref.vsv(i);
card.vph(i) = card.vsh(i) .* card_ref.vph(i) ./ card_ref.vsh(i);
card.rho(i) = v .* card_ref.rho(i) ./ card_ref.vsv(i);
end
