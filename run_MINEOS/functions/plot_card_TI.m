function plot_card_TI(cards, styles, names, zmax, plotQ)
% plot_card_TI(cards, styles, names, zmax, plotQ)
% Plot Vsv, Vsh, Vpv, Vph, rho and eta of MINEOS card structures
% (read_model_card.m) above zmax (km) in a row of subplots of the current
% figure. With plotQ = 1 (OPTIONAL, default 0), Qmu and Qkappa are plotted
% in a 7th panel.
% cards: cell array of cards; styles: cell array of line styles;
% names: cell array of legend entries
%
% jbrussell 11/2026
if nargin < 5 || isempty(plotQ)
    plotQ = 0;
end
npan = 6 + (plotQ ~= 0);
flds = {'vsv','vsh','vpv','vph','rho','eta'};
labels = {'V_{SV} (km/s)','V_{SH} (km/s)','V_{PV} (km/s)','V_{PH} (km/s)','\rho (g/cm^3)','\eta'};
for k = 1:6
    subplot(1,npan,k); box on; hold on;
    for ic = 1:length(cards)
        c = cards{ic}; Iz = c.z <= zmax;
        v = c.(flds{k})(Iz);
        if k <= 5; v = v/1000; end
        plot(v, c.z(Iz), styles{ic}, 'linewidth', 2);
    end
    xlabel(labels{k});
    if k==1; ylabel('Depth (km)'); title('Starting and true models'); legend(names,'location','southwest'); end
    ylim([0 zmax]);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
end
if plotQ
    subplot(1,npan,npan); box on; hold on;
    c = cards{1}; Iz = c.z <= zmax;
    h1 = plot(c.qmu(Iz), c.z(Iz), '-r', 'linewidth', 2);
    h2 = plot(c.qkap(Iz), c.z(Iz), '-b', 'linewidth', 2);
    xlabel('Q');
    legend([h1 h2], {'Q_\mu','Q_\kappa'}, 'location', 'southeast');
    ylim([0 zmax]);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','xscale','log');
end
end
