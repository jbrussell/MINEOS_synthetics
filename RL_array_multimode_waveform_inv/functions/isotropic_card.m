function c = isotropic_card(c)
% c = isotropic_card(c)
% Isotropic version of a MINEOS card (read_model_card.m): at every node
% Vsv = Vsh = Voigt-average Vs = sqrt((2 Vsv^2 + Vsh^2)/3),
% Vpv = Vph = Voigt-average Vp = sqrt((Vpv^2 + 4 Vph^2)/5), eta = 1.
vs = sqrt((2*c.vsv.^2 + c.vsh.^2)/3);
vp = sqrt((c.vpv.^2 + 4*c.vph.^2)/5);
c.vsv = vs; c.vsh = vs; c.vpv = vp; c.vph = vp; c.eta = ones(size(c.eta));
end
