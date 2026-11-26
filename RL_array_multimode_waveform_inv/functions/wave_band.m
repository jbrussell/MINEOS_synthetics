function b = wave_band(fc, wv)
% b = wave_band(fc, wv)
%
% Band-pass corners [f1 f2 f3 f4] (Hz) of wave type wv from opt.fc:
% opt.fc is one band [f1 f2 f3 f4] for both wave types, or two rows,
% [Love; Rayleigh]. wv = 1 (Love, transverse), 2 (Rayleigh, vertical), or
% 0: the band that covers both (lowest f1, f2; highest f3, f4), used for
% the common frequency grid and the dispersion nodes.
if isempty(fc), fc = [1/112 1/100 1/20 1/17]; end
if size(fc, 2) ~= 4 || size(fc, 1) > 2, error('wave_band: opt.fc must be [f1 f2 f3 f4] or [Love; Rayleigh] (2 x 4), Hz'); end
if wv == 0
    b = [min(fc(:, 1)) min(fc(:, 2)) max(fc(:, 3)) max(fc(:, 4))];
else
    b = fc(min(wv, size(fc, 1)), :);
end
end
