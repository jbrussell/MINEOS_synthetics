function [rec, cfg, Cref_f, CrefL, CrefR] = apply_reference(rec, cfg, P, Tdisp, fb)
% [rec, cfg, Cref_f, CrefL, CrefR] = apply_reference(rec, cfg, P, Tdisp, fb)
%
% Install the reference model of the prep structure P (prep_mineos_idagrn.m)
% in the records and the configuration of the step-2 waveform inversion:
%   rec(r).Uref   : reference spectra (idagrn6_mask branches 0..nm-1 of
%                   P.card_ref), processed like the data (time segment of
%                   the record, Tukey window, band-pass taper) and divided by
%                   the record's data normalization rec(r).scale
%   rec(r).cref   : reference phase velocities in the band
%   rec(r).refsyn : reference prediction (sum of the branches) in the windows
%   cfg.lin       : MINEOS phase velocities and kernels of P.card_ref
%                   (mineos_lin_setup.m), cfg.need, cfg.vbot, cfg.xbot, cfg.ebot,
%                   cfg.card_ref (the reference card)
% The data, the time segments (rec.tred), the windows (rec.idx, rec.wt)
% and the normalizations (rec.scale) are NOT changed, so the likelihoods of
% different references are evaluated on the same data.
%
% Used by sub_setup_records.m (initial reference) and for reference
% updates (end of burn-in in D3, D6).
%
% Cref_f       : {Love, Rayleigh} reference phase velocities at the band frequencies
% CrefL, CrefR : reference phase velocities at the dispersion nodes Tdisp
nm = P.nm;
cfg.lin = mineos_lin_setup(P.card_ref, P.K_ref_T, P.K_ref_S, P.phv_ref_T, P.phv_ref_S, ...
                           cfg.zw, cfg.zmin, cfg.zmax);
cfg.need = isfinite(P.phv_ref_T) & isfinite(P.phv_ref_S);
[~, ib_] = max(cfg.lin.z);
cfg.vbot = cfg.lin.vsv0(ib_); cfg.xbot = (cfg.lin.vsh0(ib_)/cfg.lin.vsv0(ib_))^2;
cfg.ebot = cfg.lin.eta0(ib_);
cfg.card_ref = P.card_ref;
CrefL = P.phv_ref_T; CrefR = P.phv_ref_S;
Cref_f = {disp_to_band(CrefL, Tdisp, fb), disp_to_band(CrefR, Tdisp, fb)};

uname = {'T', 'S'};
segtap = tukeywin_(cfg.N, 0.1);
tapers = {bp_taper(fb, wave_band(cfg.fc, 1)), bp_taper(fb, wave_band(cfg.fc, 2))};   % band of each wave type
Nst = size(cfg.xy, 1);
for r = 1:numel(rec)
    wv = rec(r).wave; e = rec(r).event;
    t  = (0:cfg.N-1)'*cfg.dt + rec(r).tred;
    it = round((t - P.t(1))/cfg.dt) + 1;
    Uref = cell(1, nm);
    Y = zeros(cfg.N, Nst);
    for n = 1:nm
        X = fft(P.(['u_ref_' uname{wv}]){e, n}(it, :) .* segtap);
        Uref{n} = X(cfg.ib, :) .* tapers{wv} / rec(r).scale;
        Y(cfg.ib, :) = Y(cfg.ib, :) + Uref{n};
    end
    rec(r).Uref = Uref;
    rec(r).cref = Cref_f{wv};
    u = 2*real(ifft(Y)); rec(r).refsyn = u(rec(r).idx);
end
end
