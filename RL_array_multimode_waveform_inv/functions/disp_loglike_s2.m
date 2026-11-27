function L = disp_loglike_s2(CR, cfg, lam)
% L = disp_loglike_s2(CR, cfg, lam)
%
% Log-likelihood of the fundamental-mode Rayleigh phase velocities measured
% for step 1 (opt.s2_disp: the step-1 data of opt.s1_file inside the step-2
% dispersion periods) for a step-2 state, added to the waveform
% log-likelihood (joint dispersion + waveform inversion):
%
%   L = w * sum_i [ -log(lam*sd_i) - 1/2 ((c_i(m) - cobs_i)/(lam*sd_i))^2 ]
%
% c_i(m): the state's branch-0 Rayleigh phase velocity CR(:, 1) (at
% cfg.Tdisp, linearized MINEOS) interpolated (pchip in log period) to the
% data periods cfg.s2d.T; sd_i the data errors (cfg.s2d.sd); lam the sampled
% error scale (log-uniform hyperprior on cfg.s2d.lrange); w = cfg.s2d.w
% (opt.s2_disp_weight, default 1; < 1 down-weights the dispersion, e.g. when
% it is measured from the same seismograms as the waveforms).
% Without cfg.s2d (opt.s2_disp false): L = 0.
if ~isfield(cfg, 's2d') || isempty(cfg.s2d), L = 0; return; end
D = cfg.s2d;
c = interp1(log(cfg.Tdisp(:)), CR(:, 1), log(D.T(:)), 'pchip');
r = (c - D.c(:)) ./ (lam*D.sd(:));
L = D.w * (-sum(log(lam*D.sd(:))) - 0.5*sum(r.^2));
end
