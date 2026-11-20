function plot_card_TI(cards, styles, names, zmax)
% plot_card_TI(cards, styles, names, zmax)
% Plot Vsv, Vsh, Vpv, Vph, rho and eta of MINEOS card structures
% (read_model_card.m) above zmax (km) in a 1 x 6 row of subplots of the
% current figure.
% cards: cell array of cards; styles: cell array of line styles;
% names: cell array of legend entries
%
% jbrussell 11/2026
flds = {'vsv','vsh','vpv','vph','rho','eta'};
labels = {'V_{SV} (km/s)','V_{SH} (km/s)','V_{PV} (km/s)','V_{PH} (km/s)','\rho (g/cm^3)','\eta'};
for k = 1:6
    subplot(1,6,k); box on; hold on;
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
end
