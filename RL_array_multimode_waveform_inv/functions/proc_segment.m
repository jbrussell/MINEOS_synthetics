function u = proc_segment(U, segtap, taper, cfg)
% u = proc_segment(U, segtap, taper, cfg)
% Processing of a record segment ([cfg.N x Nst] traces) as for the data and
% the reference spectra: Tukey taper segtap, one-sided band-pass taper on
% the band cfg.ib, back to real traces.
X = fft(U .* segtap);
Y = zeros(cfg.N, size(U, 2)); Y(cfg.ib, :) = X(cfg.ib, :) .* taper;
u = 2*real(ifft(Y));
end
