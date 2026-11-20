% Synthetic test: inversion of Rayleigh and/or Love wave phase velocities
% for a radially anisotropic (transversely isotropic, TI) model using
% MINEOS kernels and phase velocities.
%
% This is the MINEOS version of test_kernel_inversion_vsvh_vpvh_rho_eta.m
% in SURF96_Rayl_1Dinv (https://github.com/jbrussell/SURF96_Rayl_1Dinv),
% which uses the CPS TI programs. The model, synthetic data, inversion
% parameters and figures are the same, so the two can be compared directly.
%
% The inversion (run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m) solves for
% Vsv, Vsh, Vpv, Vph, Rho, and Eta at each card node shallower than zmax,
% enforcing Vpv/Vsv, Vph/Vsh, and Rho/Vsv of the starting model and damping
% xi = (Vsh/Vsv)^2 and eta toward the starting model. The kernels include
% all 6 parameters. MINEOS is rerun to recalculate the phase velocities and
% kernels after each iteration (nit_recalc_c, nit_recalc_G).
%
% Rayleigh + Love, Rayleigh-only, or Love-only phase velocities can be
% inverted (data_types below). Rayleigh waves are mostly sensitive to Vsv
% and Love waves to Vsh, so xi is only constrained by the data when both
% are used; with a single wave type xi is held close to the starting model
% (eps_xi large).
%
% Synthetic test: the true model is the starting model with Vsv, Vsh, Vpv,
% Vph, and Rho increased by 5% above 80 km (below the water) and radial
% anisotropy added between 20 and 200 km by raising Vsh and Vph by
% sqrt(1.15) (xi = (Vsh/Vsv)^2 = 1.15, phi = (Vpv/Vph)^2 = 1/1.15).
% The perturbations are applied to the card nodes in these depth ranges
% (the SURF96 version applies them to its layers, which are the averages
% of adjacent card nodes, so the edges of the perturbations differ by up
% to one node spacing).
%
% Differences from the SURF96 version:
% - The model is the full spherical MINEOS card; the SURF96 version uses
%   an earth-flattened layered model that ends in a half-space at zmax.
%   Only the nodes shallower than zmax are inverted for here.
% - Attenuation (use_Q = 1) is applied the same way as in the SURF96
%   version: the causal (physical dispersion) correction is applied to the
%   elastic MINEOS phase velocities and kernels to first order
%   (mineos_apply_Q.m), with the card velocities defined at fref. Q is
%   held fixed. However, Qmu and Qkappa are taken from the .qmod file (the
%   Q model MINEOS uses), whereas the SURF96 version uses the Q columns of
%   the card, which differ for this card (e.g. Qmu = 110 at 70-200 km in
%   the .qmod, 300 in the card).
%
% No need to comment out "clear" in a1_run_mineos_check.m or
% a2_mk_kernels_cv.m: run_mineos_phv_kernels.m runs them in their own
% workspace and passes the card, wave type and periods to
% parameter_FRECHET.m through environment variables.
%
% jbrussell 11/2026
%
clear; close all;

parameter_FRECHET; % paths to the functions and FORTRAN binaries

%% Data to invert
data_types = 'RL'; % 'RL' = Rayleigh + Love, 'R' = Rayleigh only, 'L' = Love only

%% Attenuation
% use_Q = 1: include Qp and Qs (from Qmu and Qkappa of the .qmod file) in
% the true and starting models, so that the synthetic data, the forward
% calculations and the kernels include the causal (physical dispersion)
% correction. Q is held fixed (not inverted). use_Q = 0: elastic.
use_Q = 0;
fref = 1; % [Hz] reference frequency of the Q model
if use_Q
    fref_Q = fref;
else
    fref_Q = []; % elastic (see run_mineos_phv_kernels.m)
end

%% Inversion regularization parameters
eps_data = 1; % data fit
eps_H = 0.02; % norm damping
eps_J = 0.001; % first derivative smoothing
eps_F = 2; % second derivative smoothing
eps_vpvs = 10; % enforce Vpv/Vsv and Vph/Vsh ratios
eps_rhovs = 10; % enforce rho/Vsv ratio
eps_eta = 10; % damp eta to starting model
if length(data_types) == 2
    eps_xi = 0.01; % damp xi=(Vsh/Vsv)^2 to starting model (weakly: R + L constrain xi)
else
    eps_xi = 100; % with a single wave type, xi is not constrained by the data
end

vsv_thresh_xi = 4; % [km/s] hold xi at the starting model where Vsv < this ("crust"); 0 = off

z_dampbot = 250; % [km] damp to starting model below this depth
zmax = 400; % [km] maximum depth of model space

nmode = 0; % mode branch number (0=fund; the kernels are only set up for the fundamental mode)

% Other inversion parameters
nit = 2; % total number of iterations
nit_recalc_c = 1; % number of iterations after which to recalculate the phase velocities
nit_recalc_G = nit_recalc_c; % number of iterations after which to recalculate the kernels

%% Generate the synthetic dataset for this test
% MINEOS card (in ./CARDS/, with its .qmod file)
cardname = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_Line1_dist100.00km';
cardid_true = 'c2_true'; % card written for the true model
cardid_work = 'c2_inv'; % card written for the recalculations during the inversion

card_start = read_model_card(['./CARDS/',cardname,'.card']);
% discontinuities (pairs of nodes at the same depth) in the inverted region
zc = card_start.z;
discs = unique(zc(find(diff(zc)==0 & zc(1:end-1)<=zmax)));
% Break the smoothing at the discontinuities of the card (water bottom,
% sediments, Moho). Smoothing across the Moho prevents the inversion from
% fitting the short-period (crust-sensitive) data, and the misfit then
% leaks into the mantle anisotropy. Set discs = [] to smooth across them.

% perturb true model
card_true = card_start;
I_pert = find(card_true.z <= 80 & card_true.vsv > 0);
for fld = {'vpv','vsv','rho','vph','vsh'}
    card_true.(fld{1})(I_pert) = card_true.(fld{1})(I_pert)*1.05;
end
% add radial anisotropy between 20 and 200 km: Vsh and Vph are both raised
% by the factor sqrt(xi_true), so that xi = (Vsh/Vsv)^2 = xi_true and
% phi = (Vpv/Vph)^2 = 1/xi_true. (The inversion ties Vph to Vsh, so the P
% anisotropy of the true model is consistent with that parameterization.)
xi_true = 1.15;
I_ani = find(card_true.z >= 20 & card_true.z <= 200 & card_true.vsv > 0);
card_true.vsh(I_ani) = card_true.vsh(I_ani)*sqrt(xi_true); % Vsh
card_true.vph(I_ani) = card_true.vph(I_ani)*sqrt(xi_true); % Vph
write_MINEOS_mod(card_true, ['./CARDS/',cardid_true,'.card']);
copyfile(['./CARDS/',cardname,'.qmod'], ['./CARDS/',cardid_true,'.qmod']);

% GENERATE SYNTHETIC DATASET
% Calculate dispersion for the true model (which is like the data)
vec_T = logspace(log10(10),log10(150),15);
periods = []; wavetype = '';
for wt = data_types
    periods = [periods; vec_T(:)];
    wavetype = [wavetype; repmat(wt, length(vec_T), 1)];
end
cobs = zeros(size(periods));
I_R = wavetype=='R'; I_L = wavetype=='L';
if any(I_R); cobs(I_R) = run_mineos_phv_kernels(cardid_true,'R',periods(I_R),nmode,0,[],fref_Q); end % "observations"
if any(I_L); cobs(I_L) = run_mineos_phv_kernels(cardid_true,'L',periods(I_L),nmode,0,[],fref_Q); end
cstd = cobs * 0.01; % observation uncertainties

%% Starting model
figure(1); clf; set(gcf,'position',[50 50 1200 600]);
plot_card_TI({card_start, card_true}, {'-b','-k'}, {'start','true'}, zmax, use_Q, ['./CARDS/',cardname,'.qmod']);

%% Do linearized inversion using MINEOS
[card_final, cpre, unc, Kfinal] = run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,wavetype,card_start,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,zmax,nit,nit_recalc_c,nmode,cardid_work,1,vsv_thresh_xi,nit_recalc_G,use_Q,fref);
cstart = unc.cstart; % phase velocities of the starting model

%% Plot final kernels
lgd = {};
Npers = length(vec_T);
clr = jet(Npers);
flds = {'dcdvsv','dcdvsh','dcdvpv','dcdvph','dcdrho','dcdeta'};
xl = {'dc/dVsv','dc/dVsh','dc/dVpv','dc/dVph','dc/d\rho','dc/d\eta'};
I = unc.I;
dz_node = abs(gradient(card_final.rad))/1000; % km
dz_node = dz_node(I);
for iw = 1:length(data_types)
    wt = data_types(iw);
    figure(100+iw); clf;
    set(gcf,'position',[100 100 1500 450],'color','w')
    for k = 1:6
        subplot(1,6,k); box on; hold on;
        for ip = 1:Npers
            plot(Kfinal.(wt).(flds{k})(:,ip)./dz_node,Kfinal.z,'-','color',clr(ip,:),'linewidth',2); % per km
            lgd{ip} = [num2str(vec_T(ip),'%.1f'),' s'];
        end
        xlabel(xl{k});
        if k==1; ylabel('Depth (km)'); title([wt,' kernels, final model']); end
        ylim([0 zmax]);
        set(gca,'fontsize',13,'linewidth',1.5,'ydir','reverse');
    end
    pos = get(gca,'Position');
    legend(lgd,'location','eastoutside');
    set(gca,'position',pos);
end

%% Plot data fit
figure(103); clf; set(gcf,'position',[300 300 1000 420]);
for iw = 1:length(data_types)
    wt = data_types(iw);
    I_w = wavetype == wt;
    subplot(1,length(data_types),iw); box on; hold on;
    plot(periods(I_w),cstart(I_w),'-ob','linewidth',2);
    plot(periods(I_w),cpre(I_w),'-or','linewidth',2);
    he = errorbar(periods(I_w),cobs(I_w),2*cstd(I_w),'sk'); set(he,'markersize',8,'markerfacecolor','k','linewidth',2);
    legend({'c start','c final','c obs'},'Location','southeast')
    xlabel('Period (s)');
    ylabel([wt,' phase velocity (km/s)']);
    set(gca,'FontSize',15,'linewidth',1.5);
end
fprintf('Reduced chi^2: start = %.2f, final = %.2f\n', ...
    mean(((cobs-cstart)./cstd).^2), mean(((cobs-cpre)./cstd).^2));

%% Plot true and recovered models: Vsv, Vsh, Vpv, Vph, rho, eta, xi
nlayer = length(I);
xi = @(c) (c.vsh(I)./c.vsv(I)).^2;
% xi uncertainty from the model covariance: xi = Vsh^2/Vsv^2,
%   dxi/dVsv = -2 xi/Vsv, dxi/dVsh = 2 xi/Vsh
xi_final = xi(card_final);
xi_std = nan(nlayer,1);
vsv_f = card_final.vsv(I)/1000; vsh_f = card_final.vsh(I)/1000;
for k = 1:nlayer
    if vsv_f(k) > 0
        Jx = [-2*xi_final(k)/vsv_f(k), 2*xi_final(k)/vsh_f(k)];
        ik = [k, nlayer+k];
        xi_std(k) = sqrt(Jx*unc.Cm(ik,ik)*Jx');
    end
end

vals = {@(c) c.vsv(I)/1000, @(c) c.vsh(I)/1000, @(c) c.vpv(I)/1000, @(c) c.vph(I)/1000, @(c) c.rho(I)/1000, @(c) c.eta(I), xi};
stds = {unc.vsv_std, unc.vsh_std, unc.vpv_std, unc.vph_std, unc.rho_std, unc.eta_std, xi_std};
labels = {'V_{SV} (km/s)','V_{SH} (km/s)','V_{PV} (km/s)','V_{PH} (km/s)','\rho (g/cm^3)','\eta','\xi = (V_{SH}/V_{SV})^2'};
z = card_final.z(I);

figure(100); clf;
set(gcf,'position',[50 100 1700 650]);
for k = 1:7
    subplot(1,7,k); box on; hold on;
    v_true = vals{k}(card_true); v_start = vals{k}(card_start); v_final = vals{k}(card_final);
    h1 = plot(v_true, z, '-k', 'linewidth', 2);
    h2 = plot(v_start, z, '-b', 'linewidth', 2);
    h3 = plot(v_final, z, '-r', 'linewidth', 2);
    plot(v_final - stds{k}, z, '--r', 'linewidth', 1);
    plot(v_final + stds{k}, z, '--r', 'linewidth', 1);
    xlabel(labels{k});
    if k==1; ylabel('Depth (km)'); title(['Data: ',data_types]); end
    ylim([0 z_dampbot+50]);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
    if k==1; legend([h1 h2 h3],{'true','start','final'},'Location','southwest'); end
end
