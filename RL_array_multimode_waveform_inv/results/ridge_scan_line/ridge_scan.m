% RIDGE_SCAN  Ridge or separate modes? log-likelihood of the line-test data
% along straight lines between the posterior-mean models of the ensembles of
% tests 2b (random starts) and 2f (refT1_TNjitter + joint dispersion term)
% and the true model, in the 2b setup (step-1 isotropic waveform reference,
% tau fixed at the truth, B free and marginalized (nuis_relin), dth ~0).
%
% Point t on the line from A (t = 0) to B (t = 1): dispersion
% C(t) = (1-t) C_A + t C_B (Love and Rayleigh, all branches), where C_A is
% the posterior mean of the ensemble's T = 1 samples (equivalent to a line
% between the mean profiles, as the dispersion is linearized about the
% reference card: the profiles at t are (1-t) m_A + t m_B).
% Terms at each point:
%   L   waveform log-likelihood (Love + Rayleigh), noise sigma/rho = medians of
%       the posterior samples of all 8 ensembles
%   Ld  joint dispersion term of 2e/2f (error scale 1)
%   Lp  Vsv prior toward the step-1 model (vsv_ref_sig 'step1'), from the
%       mean Vsv profiles interpolated linearly
% A barrier (L dropping below both ends inside 0 < t < 1) means separate
% modes; a monotonic or concave curve means one connected ridge/basin.
addpath('functions');
S = load('setup_line2b.mat'); F = load('setup_line2f.mat');
rec = S.rec; cfg = S.cfg; mc = S.mc; nm = S.nm; Tdisp = S.Tdisp; fb = S.fb; zplot = S.zplot;
cfgd = F.cfg;                                     % has the joint dispersion data (cfg.s2d)
tests = {'line2b', 'line2f'}; E = struct('name', {}, 'CL', {}, 'CR', {}, 'vsv', {}, 'xi', {}, 'Lchain', {});
sig = []; rho = [];
for it = 1:2
    for i = 1:4
        C = load(sprintf('chain_%s_%02d.mat', tests{it}, i)); sv = C.sv;
        e = struct('name', sprintf('%s_E%d', strrep(tests{it}, 'line', ''), i), ...
                   'CL', mean(sv.CL, 3), 'CR', mean(sv.CR, 3), 'vsv', mean(sv.vsv, 2), 'xi', mean(sv.xi, 2), ...
                   'Lchain', mean(sv.L(end-numel(sv.kv)+1:end)));
        E(end+1) = e; %#ok<SAGROW>
        sig = [sig; sv.sig]; rho = [rho; sv.rho]; %#ok<AGROW>
    end
end
sig = median(sig, 1); rho = median(rho, 1);
tru = struct('name', 'truth', 'CL', S.CtrueL, 'CR', S.CtrueR, 'vsv', S.tru.vsv(:), 'xi', S.tru.xi(:), 'Lchain', NaN);
E2f = E(5:8);                                     % 2f: one converged mode, its mean
m2f = struct('name', '2f_mean', 'CL', mean(cat(3, E2f.CL), 3), 'CR', mean(cat(3, E2f.CR), 3), ...
             'vsv', mean([E2f.vsv], 2), 'xi', mean([E2f.xi], 2), 'Lchain', mean([E2f.Lchain]));
P = [{tru}, num2cell(E), {m2f}];                  % 1 truth, 2-5 2b, 6-9 2f, 10 2f mean
pairs = [ones(8, 1), (2:9)'];                     % truth -> each ensemble
for i = 2:4, for j = i+1:5, pairs(end+1, :) = [i j]; end, end   %#ok<SAGROW> 2b - 2b
for i = 2:5, pairs(end+1, :) = [10 i]; end        %#ok<SAGROW> 2f mean - 2b
t = (-0.25:0.025:1.25)'; nt = numel(t);
% Vsv prior toward step 1 on the 1 km grid zmin:zmax
zg = (cfg.zmin:cfg.zmax)'; vref = interp1(cfg.refp.z, cfg.refp.vsv, zg); sz = mc.vsv_ref_sigz(:);
lpv = @(v) -0.5*sum(((interp1(zplot, v, zg) - vref)./sz).^2)/mc.vsv_ref_len;
evalL = @(CL, CR) ridge_wave_L(rec, {disp_to_band(CL, Tdisp, fb), disp_to_band(CR, Tdisp, fb)}, S.lin0, cfg, sig, rho, mc, nm);
R = struct('t', t, 'names', {cellfun(@(p) p.name, P, 'UniformOutput', false)}, 'pairs', pairs, 'sig', sig, 'rho', rho);
R.L = NaN(nt, size(pairs, 1)); R.Ld = R.L; R.Lp = R.L; R.Ltrue = S.Ltrue;
for k = 1:size(pairs, 1)
    A = P{pairs(k, 1)}; B = P{pairs(k, 2)};
    for q = 1:nt
        CL = (1 - t(q))*A.CL + t(q)*B.CL; CR = (1 - t(q))*A.CR + t(q)*B.CR;
        R.L(q, k) = evalL(CL, CR);
        R.Ld(q, k) = disp_loglike_s2(CR, cfgd, 1);
        R.Lp(q, k) = lpv((1 - t(q))*A.vsv + t(q)*B.vsv);
    end
    i0 = find(abs(t) < 1e-9); i1 = find(abs(t - 1) < 1e-9); in = t > 0 & t < 1;
    bar = min(R.L(in, k)) - min(R.L([i0 i1], k));
    [lm, im] = max(R.L(:, k));
    fprintf('%-8s -> %-8s  L(0) %9.1f  L(1) %9.1f  max %9.1f at t=%5.2f  barrier %7.1f   Ld %6.1f..%6.1f  Lp %7.1f..%7.1f\n', ...
        A.name, B.name, R.L(i0, k), R.L(i1, k), lm, t(im), bar, R.Ld(i0, k), R.Ld(i1, k), R.Lp(i0, k), R.Lp(i1, k));
end
% xi bands of the end points (for the figure)
zb = [7 60; 60 160; 160 250; 250 350];
R.xib = NaN(numel(P), 4);
for i = 1:numel(P)
    for b = 1:4, R.xib(i, b) = mean(P{i}.xi(zplot >= zb(b, 1) & zplot < zb(b, 2))); end
    fprintf('%-8s xi bands %s   mean-model L %.1f (chain mean L %.1f)\n', P{i}.name, mat2str(round(R.xib(i, :)*1000)/1000), ...
        evalL(P{i}.CL, P{i}.CR), P{i}.Lchain);
end
R.prof = cellfun(@(p) [p.vsv, p.xi], P, 'UniformOutput', false); R.zplot = zplot;
save('ridge_scan.mat', '-struct', 'R');

