function compare_runs(tags, names, outfile)
% compare_runs(tags, names, outfile)
% Posterior mean and 5-95% range of Vsv, Vsh and xi for several runs
% (setup<tag>.mat, chain<tag>_XX.mat), with the true model (PA5) and the
% reference model of each run, in one figure.
%   compare_runs({'_plain','_pt'}, {'plain rjMcMC','parallel tempering'}, 'figs_compare.png')
if exist('OCTAVE_VERSION', 'builtin'), set(0, 'defaultfigurevisible', 'off'); end
nr = numel(tags);
flds = {'vsv', 'vsh', 'xi'};
labs = {'V_{SV} (km/s)', 'V_{SH} (km/s)', '\xi = (V_{SH}/V_{SV})^2'};
xr = {[3.9 5.1], [3.9 5.3], [0.88 1.22]};
fig = figure('position', [50 50 300*nr 900]);
for ir = 1:nr
    S = load(['setup' tags{ir} '.mat'], 'tru', 'refp', 'mc', 'zplot', 'cfg');
    X = struct('vsv', [], 'vsh', [], 'xi', []);
    for ic = 1:S.mc.nchains
        f = sprintf('chain%s_%02d.mat', tags{ir}, ic);
        if exist(f, 'file') ~= 2, continue; end
        D = load(f);
        for q = 1:3, X.(flds{q}) = [X.(flds{q}), D.sv.(flds{q})]; end
    end
    z = S.zplot;
    for q = 1:3
        subplot(3, nr, (q-1)*nr + ir); hold on; box on;
        lo = prctile_(X.(flds{q}), 5); hi = prctile_(X.(flds{q}), 95);
        fill([lo; flipud(hi)], [z; flipud(z)], [0.75 0.82 1], 'edgecolor', 'none');
        h2 = plot(S.refp.(flds{q}), z, '-', 'color', [0 0.6 0], 'linewidth', 2.5);
        h1 = plot(S.tru.(flds{q}), z, 'r-', 'linewidth', 1.2);
        h3 = plot(mean(X.(flds{q}), 2), z, 'b-', 'linewidth', 1.5);
        plot(xr{q}, S.cfg.zmax*[1 1], 'k:');
        set(gca, 'ydir', 'reverse'); ylim([0 380]); xlim(xr{q}); xlabel(labs{q});
        if ir == 1, ylabel('Depth below seafloor (km)'); end
        if q == 1, title(names{ir}); end
        if q == 1 && ir == 1
            legend([h1 h2 h3], {'true (PA5)', 'reference', 'posterior mean (5-95%)'}, ...
                'location', 'southwest', 'fontsize', 7);
        end
    end
end
print(fig, '-dpng', '-r110', outfile);
end
