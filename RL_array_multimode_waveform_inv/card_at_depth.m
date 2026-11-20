function v = card_at_depth(card, f, zq, side)
% v = card_at_depth(card, f, zq, side)
% Value of field f of a MINEOS card (read_model_card.m) at depths zq (km
% below sea level), linear between nodes. At a discontinuity (two nodes at
% the same depth) side = 'below' (default) takes the deeper (lower) side
% and 'above' the shallower side.
if nargin < 4, side = 'below'; end
z = flipud(card.z(:)); y = flipud(card.(f)(:));     % surface first, depth increasing
v = zeros(size(zq));
for i = 1:numel(zq)
    k = find(z <= zq(i));                            % nodes at or above zq
    if isempty(k), v(i) = y(1); continue; end
    if strcmp(side, 'below'), k = k(end); else, k = find(z >= zq(i), 1, 'first'); k = k(1); end
    if z(k) == zq(i) || k == numel(z)
        v(i) = y(k);
    else
        if strcmp(side, 'above'), k = k - 1; end
        v(i) = y(k) + (y(k+1) - y(k)) * (zq(i) - z(k)) / (z(k+1) - z(k));
    end
end
end
