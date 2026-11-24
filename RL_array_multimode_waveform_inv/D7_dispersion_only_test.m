% D7_DISPERSION_ONLY_TEST  Diagnostic (synthetic test): invert the TRUE
% phase velocities directly, without waveforms or nuisance parameters, with
% the model space, prior, proposals and MINEOS linearization of a step-2
% run (setup<opt.disp_tag>.mat, made by D3; functions/run_disp_inv.m).
%
% Two data sets, Love and Rayleigh branches 0..nm-1 at the dispersion
% periods, relative standard deviation opt.disp_sig (default 0.25%):
%   'lin'    : the linearized dispersion of the true model (true Vsv, Vsh
%              between zmin and zmax, the rest at the reference) - exactly
%              representable by the linearization: tests the resolution of
%              the dispersion and the parameterization only
%   'mineos' : the exact MINEOS dispersion of the true model - adds the
%              linearization error
% If xi is recovered here but not by the waveform inversion, the problem is
% in the waveform/nuisance part (time shifts, branch interference); if not,
% in the resolution, the parameterization or the linearization.
%
% Options: disp_tag ('_plain8lin'), disp_sig (0.0025), disp_nit (20000),
% disp_burn (10000), disp_chains (4), disp_w ([1 1]: weights of Love and
% Rayleigh), disp_cases ({'lin', 'mineos'}), disp_branches ([] = all; e.g.
% 0 or [0 1]: only these branches are data), disp_out (suffix of the
% output folder); summary.txt is appended to. Output: figs_disp<disp_tag>/ (disp_<case>.png, summary.txt).

parameters_RL;
tg_ = optget(opt, 'disp_tag', '_plain8lin');
S0 = load(['setup' tg_ '.mat'], 'cfg', 'mc', 'Tdisp', 'fb', 'nm', 'zplot', 'tru', 'refp', ...
           'CtrueL', 'CtrueR', 'ClinL', 'ClinR');
sigd = optget(opt, 'disp_sig', 0.0025);
nit = optget(opt, 'disp_nit', 20000); burn = optget(opt, 'disp_burn', 10000);
nch = optget(opt, 'disp_chains', 4); wts = optget(opt, 'disp_w', [1 1]); thin = 20;
br_ = optget(opt, 'disp_branches', []);       % branches fitted ([] = all 0..nm-1)
od = ['figs_disp' tg_ optget(opt, 'disp_out', '')]; if ~exist(od, 'dir'), mkdir(od); end
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
cases = optget(opt, 'disp_cases', {'lin', 'mineos'});
z = S0.zplot; tru = S0.tru; refp = S0.refp;
bands = [7 60; 60 160; 160 250; 250 350];
fid = fopen(fullfile(od, 'summary.txt'), 'a');
for ic = 1:numel(cases)
    if strcmp(cases{ic}, 'lin'), CLd = S0.ClinL; CRd = S0.ClinR; else, CLd = S0.CtrueL; CRd = S0.CtrueR; end
    if ~isempty(br_)                               % only these branches are data
        off = setdiff(0:S0.nm-1, br_) + 1; CLd(:, off) = NaN; CRd(:, off) = NaN;
    end
    R = cell(nch, 1); tic;
    for k = 1:nch
        R{k} = run_disp_inv(CLd, CRd, sigd, S0.cfg, S0.mc, S0.Tdisp, S0.fb, S0.nm, z, nit, burn, thin, 100*find(strcmp(cases{ic}, {'lin', 'mineos'})) + k, wts);
    end
    V = struct('vsv', [], 'vsh', [], 'xi', []);
    for f_ = {'vsv', 'vsh', 'xi'}, V.(f_{1}) = cell2mat(cellfun(@(x) x.(f_{1}), R', 'UniformOutput', false)); end
    % summary
    for o = [1 fid]
        fprintf(o, 'Case %s (data: %s dispersion of the true model, branches %s), sig %.2f%%, %d chains x %d it (%.0f s); acceptance %s\n', ...
            cases{ic}, cases{ic}, mat2str(ifelse_br(br_, S0.nm)), 100*sigd, nch, nit, toc, mat2str(round(100*cellfun(@(x) x.acc, R'))/100));
        fprintf(o, '   final log L per chain %s (0 = exact fit)\n', mat2str(round(cellfun(@(x) x.L(end), R'))));
        fprintf(o, '   xi band averages: posterior mean [5-95%%] | true | chain means\n');
        for b = 1:size(bands, 1)
            ib = z >= bands(b, 1) & z < bands(b, 2);
            xs = mean(V.xi(ib, :), 1); q = sort(xs); n = numel(q);
            cm = cellfun(@(x) mean(mean(x.xi(ib, :))), R');
            fprintf(o, '   %3d-%3d km: %.3f [%.3f-%.3f] | %.3f | %s\n', bands(b, :), mean(xs), ...
                q(max(1, round(0.05*n))), q(round(0.95*n)), mean(tru.xi(ib)), mat2str(round(cm*1000)/1000));
        end
        for f_ = {'vsv', 'vsh', 'xi'}
            X = V.(f_{1}); in = z >= S0.cfg.zmin & z <= S0.cfg.zmax; q = sort(X(in, :), 2); n = size(q, 2);
            lo = q(:, max(1, round(0.05*n))); hi = q(:, round(0.95*n)); t_ = tru.(f_{1})(in);
            fprintf(o, '   %-4s truth in 5-95%%: %3.0f%% | mean|err| %.3f\n', f_{1}, 100*mean(t_ >= lo & t_ <= hi), mean(abs(mean(X(in, :), 2) - t_)));
        end
    end
    % figure: posteriors (density, 5-95%, mean, chain means, truth, reference)
    fig = figure('position', [0 0 1500 650]);
    nm_ = {'vsv', 'vsh', 'xi'}; lab = {'V_{SV} (km/s)', 'V_{SH} (km/s)', '\xi = (V_{SH}/V_{SV})^2'};
    xl = {[4 5], [4 5.1], [0.88 1.22]}; col = lines(nch);
    for q_ = 1:3
        subplot(1, 3, q_); hold on; box on;
        X = V.(nm_{q_}); s = sort(X, 2); n = size(s, 2);
        plot(s(:, max(1, round(0.05*n))), z, 'b--', s(:, round(0.95*n)), z, 'b--');
        for k = 1:nch, plot(mean(R{k}.(nm_{q_}), 2), z, '-', 'color', col(k, :), 'linewidth', 0.6); end
        h1 = plot(mean(X, 2), z, 'b-', 'linewidth', 2);
        h2 = plot(tru.(nm_{q_}), z, 'r-', 'linewidth', 1.5);
        h3 = plot(refp.(nm_{q_}), z, 'g-', 'linewidth', 1);
        set(gca, 'ydir', 'reverse'); ylim([0 400]); xlim(xl{q_}); xlabel(lab{q_}); ylabel('Depth below seafloor (km)');
        plot(xl{q_}, S0.cfg.zmax*[1 1], 'k:');
        if q_ == 1, legend([h2 h3 h1], {'true (PA5)', 'reference (step 1)', 'posterior mean'}, 'location', 'southwest'); end
    end
    axes('position', [0 0.95 1 0.04], 'visible', 'off');
    text(0.5, 0.5, sprintf(['Dispersion-only inversion of the %s dispersion of the true model (Love + Rayleigh branches %s, sig %.2f%%); ' ...
        'thin colour lines: chain means, dashed: 5-95%%'], cases{ic}, mat2str(ifelse_br(br_, S0.nm)), 100*sigd), 'horizontalalignment', 'center');
    print(fig, '-dpng', '-r110', fullfile(od, ['disp_' cases{ic} '.png'])); close(fig);
    save(fullfile(od, ['disp_' cases{ic} '.mat']), 'R', '-v7');
end
fclose(fid);
fprintf('figures and summary in %s/\n', od);
