function sigz = vsv_ref_sigz(opt, cfg)
% sigz = vsv_ref_sigz(opt, cfg)
%
% Depth-dependent width (km/s) of the prior toward the reference Vsv on the
% 1 km grid cfg.zmin:cfg.zmax (lr_common 'logprior', mc.vsv_ref_sigz), or []
% when opt.vsv_ref_sig is a number (scalar width, or 0 = off).
%
% opt.vsv_ref_sig = 'step1': the posterior std of the step-1 Vsv
% (s1_unc.vsv_std in opt.s1_file, from D1_step1_Bayesian_disp.m or
% D1_step1_Bayesian_waveform.m; the linearized D1 also writes one, but it
% is not a posterior and is typically too narrow), times
% opt.vsv_ref_scale (default 1), with a floor opt.vsv_ref_sig_min
% (default 0.01 km/s: the step-1 std goes to 0 where its model is tied to
% the starting model, at the bottom). The reference model of step 2 should
% then be the step-1 model (opt.refid = opt.s1_final), so the prior is the
% step-1 posterior carried forward (Gaussian approximation, independent
% every opt.vsv_ref_len km).
sigz = [];
v = optget(opt, 'vsv_ref_sig', 0);
if isnumeric(v), return; end
if ~strcmp(v, 'step1'), error('vsv_ref_sigz: opt.vsv_ref_sig must be a number or ''step1'''); end
S = load(opt.s1_file);
if ~isfield(S, 's1_unc'), error('vsv_ref_sigz: %s has no s1_unc (rerun step 1)', opt.s1_file); end
if isfield(S, 'method') && strcmp(S.method, 'linearized')
    warning('vsv_ref_sigz: the step-1 uncertainty of %s is linearized (not a posterior; usually too narrow)', opt.s1_file);
end
U = S.s1_unc; ok = isfinite(U.vsv_std);
zg = (cfg.zmin:1:cfg.zmax)';
sigz = interp1(U.z(ok), U.vsv_std(ok), zg, 'linear', 'extrap') * optget(opt, 'vsv_ref_scale', 1);
sigz = max(sigz, optget(opt, 'vsv_ref_sig_min', 0.01));
end
