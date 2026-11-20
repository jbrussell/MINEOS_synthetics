function tp = bp_taper(f, fc)
% BP_TAPER  Cosine band-pass taper with corners fc = [f1 f2 f3 f4] (Hz).
tp = double(f > fc(1) & f < fc(4));
lo = f > fc(1) & f < fc(2); tp(lo) = 0.5 - 0.5*cos(pi*(f(lo)-fc(1))/(fc(2)-fc(1)));
hi = f > fc(3) & f < fc(4); tp(hi) = 0.5 + 0.5*cos(pi*(f(hi)-fc(3))/(fc(4)-fc(3)));
end
