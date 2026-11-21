function Cf = disp_to_band(C, Tdisp, fb)
% DISP_TO_BAND  Interpolate dispersion nodes (periods Tdisp) to band
% frequencies fb, mode by mode, using only nodes where the mode exists.
% Frequencies outside a mode's existing nodes are NaN (masked in event_G).
x = log(Tdisp(:)); xq = log(1./fb(:));
Cf = NaN(numel(fb), size(C, 2));
for n = 1:size(C, 2)
    ok = ~isnan(C(:, n));
    if nnz(ok) >= 2
        Cf(:, n) = interp1(x(ok), C(ok, n), xq, 'pchip', NaN);
    end
end
end
