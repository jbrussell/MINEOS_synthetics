% Synthetic Rayleigh wave phase velocity inversion using MINEOS kernels
%
% The inversion solves for 6 parameters: Vsv, Vsh, Vpv, Vph, Rho, and Eta.
% Constraint equations include:
% - Norm damping
% - First derivative (flatness) smoothing
% - Second derivative (roughness) smoothing
% - Damping Vp(v,h)/Vs(v,h) ratios toward that of starting model
% - Damping Rho/Vsv ratio toward that of starting model
% - Strong Damping xi = (Vsh/Vsv)^2 toward that of starting model
% - Strong Damping crust (Vsv < 4.0 km/s) toward that of starting model
% - Strong Damping of water layer
%
% This inversion does not currently recalculate kernels or dispersion, 
% though the framework to do so is there... That will be included in the future.
%
% !! a1_run_mineos_check.m: Must comment out "clear" command at the top
% prior to running this!
%
% !! a2_mk_kernels_cv.m: Must comment out "clear" command at the top. Must
% also ensure the correct mode "branch" is specified at the top. it
% must match parameter "nmode" in this script.
%
% jbrussell 9/26/26
%

clear; close all;

parameter_FRECHET; % Set the base card file you want to use in the the parameter file
CARD = param.CARD;
CARDPATH = param.CARDPATH;
CARDID = param.CARDID;
CARDTABLE = CARDTABLE;
TYPE = param.TYPE;
SID = param.STYPEID;
TID = param.TTYPEID;

periods = param.periods; % round(logspace(log10(20),log10(150),10));

%% Inversion parameters
% Regularization
eps_data = 1000; %1; % data fit
eps_H = 0.05; % norm damping
eps_J = 0.001; % first derivative smoothing
eps_F = 1; % second derivative smoothing
eps_vpvs = 10; % enforce vp/vs ratio
eps_rhovs = 10; % enforce rho/vs ratio

% Break smoothing constraints at the following depths. If don't want to 
% allow discontinuities, then set to empty vector []. You could also set
% this below after reading in the input card file, if you want the
% discontinuties contained in the starting model to be respected.
discs = []; % [km] Discontinuties.

% The inversion technically solves for the whole planet. Define the
% effective maximum depth of the inverted region.
z_dampbot = 400; % [km] damp to starting model below this depth.

% The crust can sometimes dominate the inversion. This flag allows you to
% keep the crust fixed while solving for everything below it. Crust is 
% defined as Vsv < 4.0 km/s.
is_dampcrust = 1; % strongly damp crust to starting model?
vsv_thresh_dampcrust = 4; % [km/s] velocity below which is considered "crust"

% Mode of interest
nmode = 0; % mode branch number (0=fund, 1=1st overtone, 2=2nd overtone, etc.)

% Other inversion parameters
nit = 4; % total number of iterations
nit_recalc_c = 9999; % (NOT IMPLEMENTED YET!) number of iterations after which to recalculate phase velocity and kernels

%% %%%%%%%%%%%%%%%%%% SYNTHETIC SETUP BELOW THIS POINT %%%%%%%%%%%%%%%%%%
%% Generate synthetic dataset
% Save base card file to a temp file because we are going to need to perturb
% it to generate our synthetic dataset by running a1_run_mineos_check. 
% After we generate the synhtetic dataset, we will use the temp file to 
% restore the base card. That will then be our starting model.

% Copy base card file to temp
system(['cp ',CARDPATH,CARD,' temp.card']);

ncard = [CARDPATH,CARD];
card_base = read_model_card(ncard);
% Calculate dispersion for base model
a1_run_mineos_check % !! Comment out "clear" command at the top of a1
% Get dispersion
if strcmp(TYPE,'S') == 1
    CASC = [CARDTABLE,CARDID,'.',SID,'.asc'];
    qfile = [CARDTABLE,CARDID,'.',SID,'.q'];
elseif strcmp(TYPE,'T') == 1
    CASC = [CARDTABLE,CARDID,'.',TID,'.asc'];
    qfile = [CARDTABLE,CARDID,'.',TID,'.q'];
else
    disp('Type does not exist! Use T or S')
end
[phv_base,grv_base,phvq_base] = readMINEOS_qfile_per(qfile,periods,nmode);

% Load base card file and perturb it
ncard_base = [CARDPATH,CARD];
card_base = read_model_card(ncard_base);
card_pert = card_base;
% Perturb the model...
I_pert = find(card_pert.z>=30 & card_pert.z<=75);
card_pert.vsv(I_pert) = card_pert.vsv(I_pert)*1.1;
card_pert.vsh(I_pert) = card_pert.vsh(I_pert)*1.1;
card_pert.vpv(I_pert) = card_pert.vpv(I_pert)*1.1;
card_pert.vph(I_pert) = card_pert.vph(I_pert)*1.1;
card_pert.rho(I_pert) = card_pert.rho(I_pert)*1.1;
I_pert = find(card_pert.z>=80 & card_pert.z<=220);
card_pert.vsv(I_pert) = card_pert.vsv(I_pert)*0.95;
card_pert.vsh(I_pert) = card_pert.vsh(I_pert)*0.95;
card_pert.vpv(I_pert) = card_pert.vpv(I_pert)*0.95;
card_pert.vph(I_pert) = card_pert.vph(I_pert)*0.95;
card_pert.rho(I_pert) = card_pert.rho(I_pert)*0.95;

% Write perturbed card file to base location
write_MINEOS_mod(card_pert,[CARDPATH,CARD])

% Generate synthetic dataset
a1_run_mineos_check % !! Comment out "clear" command at the top of a1
% Get dispersion
if strcmp(TYPE,'S') == 1
    CASC = [CARDTABLE,CARDID,'.',SID,'.asc'];
    qfile = [CARDTABLE,CARDID,'.',SID,'.q'];
elseif strcmp(TYPE,'T') == 1
    CASC = [CARDTABLE,CARDID,'.',TID,'.asc'];
    qfile = [CARDTABLE,CARDID,'.',TID,'.q'];
else
    disp('Type does not exist! Use T or S')
end
[phv_pert,grv,phvq] = readMINEOS_qfile_per(qfile,periods,nmode);

% Copy temp card back to its base
system(['cp temp.card ',CARDPATH,CARD]);

% Plot base and perturbed
figure(1); clf; 
subplot(1,2,1); box on; hold on;
plot(card_base.vsv/1000,card_base.z,'-k','linewidth',1.5);
plot(card_pert.vsv/1000,card_pert.z,'--r','linewidth',1.5);
ylim([0 z_dampbot]);
xlabel('V_{SV} (km/s)');
ylabel('Depth (km)');
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
legend({'Base';'Perturbed'},'location','southwest')

subplot(1,2,2); box on; hold on;
plot(periods, phv_base,'-k','linewidth',1.5);
plot(periods, phv_pert,'--r','linewidth',1.5);
xlabel('Periods (s)');
ylabel('Phase Velocity (km/s)');
set(gca,'fontsize',15,'linewidth',1.5);
% legend({'Base';'Perturbed'},'location','southwest')

% Save "observed" data vectors
cobs = phv_pert(:); % [km/s] "observations"
cstd = cobs * 0.005; % [km/s] observation uncertainties

%% %%%%%%%%%%%%%%%%%% END SYNTHETIC SETUP %%%%%%%%%%%%%%%%%%

%% Starting model
% load starting model
% Get MINEOS model
ncard = [CARDPATH,CARD];
card = read_model_card(ncard);

figure(2); clf;
box on; hold on;
plot(card.vsv/1000,card.z,'-k','linewidth',1.5);
ylim([0 z_dampbot]);
xlabel('V_{SV} (km/s)');
ylabel('Depth (km)');
title('Starting Model');
set(gca,'FontSize',15,'linewidth',1.5,'ydir','reverse');

%% Do linearized inversion

% Calculate dispersion for starting model
a1_run_mineos_check % !! Comment out "clear" command at the top of a1
% Get dispersion
if strcmp(TYPE,'S') == 1
    CASC = [CARDTABLE,CARDID,'.',SID,'.asc'];
    qfile = [CARDTABLE,CARDID,'.',SID,'.q'];
elseif strcmp(TYPE,'T') == 1
    CASC = [CARDTABLE,CARDID,'.',TID,'.asc'];
    qfile = [CARDTABLE,CARDID,'.',TID,'.q'];
else
    disp('Type does not exist! Use T or S')
end
[cstart,~,cq_start] = readMINEOS_qfile_per(qfile,periods,nmode);
cstart = cstart(:);

%% Do the inversion

[card_inv, cpre, unc] = run_mineos_inv_Rayl_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,card,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,z_dampbot,is_dampcrust,vsv_thresh_dampcrust,nit,nit_recalc_c,nmode,qfile);

%%
% Plot
figure(3); clf; 
set(gcf,'position',[370   372   967   580]);

subplot(2,2,[1 3]); box on; hold on;
plot(card_pert.vsv/1000,card_pert.z,'-k','linewidth',2);
plot(card.vsv/1000,card.z,'-b','linewidth',2);
plot(card_inv.vsv/1000,card_inv.z,'-r','linewidth',2);
plot(card_inv.vsv/1000-unc.vsv_std/1000,card_inv.z,'--r','linewidth',2);
plot(card_inv.vsv/1000+unc.vsv_std/1000,card_inv.z,'--r','linewidth',2);
xlabel('V_{SV}');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
legend({'true','start','final'},'Location','southwest')
% legend({'start','final'},'Location','southwest')

subplot(2,2,2); box on; hold on;
% cpre = dispR_surf96(periods,finalmod,nmode);
plot(periods,cstart,'-ob','linewidth',2);
plot(periods,cpre,'-or','linewidth',2);
errorbar(periods,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
legend({'c start','c final','c obs'},'Location','best')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

%% Plot full model
figure(4); clf
set(gcf,'position',[54         145        1622         714]);

subplot(1,6,1); box on; hold on;
plot(card_pert.vsv/1000,card_pert.z,'-k','linewidth',2);
plot(card.vsv/1000,card.z,'-b','linewidth',2);
plot(card_inv.vsv/1000,card_inv.z,'-r','linewidth',2);
plot(card_inv.vsv/1000-unc.vsv_std/1000,card_inv.z,'--r','linewidth',2);
plot(card_inv.vsv/1000+unc.vsv_std/1000,card_inv.z,'--r','linewidth',2);
xlabel('V_{SV} (km/s)');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
legend({'true','start','final'},'Location','southwest')

subplot(1,6,2); box on; hold on;
plot(card_pert.vsh/1000,card_pert.z,'-k','linewidth',2);
plot(card.vsh/1000,card.z,'-b','linewidth',2);
plot(card_inv.vsh/1000,card_inv.z,'-r','linewidth',2);
plot(card_inv.vsh/1000-unc.vsh_std/1000,card_inv.z,'--r','linewidth',2);
plot(card_inv.vsh/1000+unc.vsh_std/1000,card_inv.z,'--r','linewidth',2);
xlabel('V_{SH} (km/s)');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');

subplot(1,6,3); box on; hold on;
plot(card_pert.vpv/1000,card_pert.z,'-k','linewidth',2);
plot(card.vpv/1000,card.z,'-b','linewidth',2);
plot(card_inv.vpv/1000,card_inv.z,'-r','linewidth',2);
plot(card_inv.vpv/1000-unc.vpv_std/1000,card_inv.z,'--r','linewidth',2);
plot(card_inv.vpv/1000+unc.vpv_std/1000,card_inv.z,'--r','linewidth',2);
xlabel('V_{PV} (km/s)');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');

subplot(1,6,4); box on; hold on;
plot(card_pert.vph/1000,card_pert.z,'-k','linewidth',2);
plot(card.vph/1000,card.z,'-b','linewidth',2);
plot(card_inv.vph/1000,card_inv.z,'-r','linewidth',2);
plot(card_inv.vph/1000-unc.vph_std/1000,card_inv.z,'--r','linewidth',2);
plot(card_inv.vph/1000+unc.vph_std/1000,card_inv.z,'--r','linewidth',2);
xlabel('V_{PH} (km/s)');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');

subplot(1,6,5); box on; hold on;
plot(card_pert.rho,card_pert.z,'-k','linewidth',2);
plot(card.rho,card.z,'-b','linewidth',2);
plot(card_inv.rho,card_inv.z,'-r','linewidth',2);
plot(card_inv.rho-unc.rho_std,card_inv.z,'--r','linewidth',2);
plot(card_inv.rho+unc.rho_std,card_inv.z,'--r','linewidth',2);
xlabel('\rho (kg/m^3)');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');

subplot(1,6,6); box on; hold on;
plot(card_pert.eta,card_pert.z,'-k','linewidth',2);
plot(card.eta,card.z,'-b','linewidth',2);
plot(card_inv.eta,card_inv.z,'-r','linewidth',2);
plot(card_inv.eta-unc.eta_std,card_inv.z,'--r','linewidth',2);
plot(card_inv.eta+unc.eta_std,card_inv.z,'--r','linewidth',2);
xlabel('\eta');
ylabel('Depth');
ylim([0 z_dampbot+100]);
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');

%% Plot Vsv Resolution kernel and Backus-Gilbert spread

nlayer = length(card_inv.z);
dr = gradient(card_inv.rad); % m
iv  = 1:nlayer; % indices for vsv
Rvv = unc.R(iv,iv);
z   = card_inv.z;

figure(5); clf; 
set(gcf,'position',[425   379   902   542])

sgtitle('V_{SV} Resolution','fontsize',20,'fontweight','bold')

subplot(1,2,1); box on; hold on;
box on; hold on
lgd5 = {};
ii = 0;
% averaging kernels: row k shows how true Vsv at all depths maps into estimated Vsv at depth z(k)
for k = find(z > 20 & z < z_dampbot, 1):10:find(z < z_dampbot, 1, 'last')
    ii = ii + 1;
    plot(Rvv(k,:)./dr(:)'*1000, z, 'linewidth', 1.5);   % divide by dr -> per-km, handles uneven layers
    lgd5{ii} = [num2str(z(k)),'km'];
end
legend(lgd5,'location','southeast');
ylim([0 z_dampbot+100]);
xlabel('Averaging kernel (1/km)');
ylabel('Depth (km)')
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');

% Calculate Backus spread
% The Backus–Gilbert spread boils the averaging kernel down to one number: 
% the depth width over which the estimate is averaging. In practice it is 
% the vertical resolution length at that depth.
subplot(1,2,2); box on; hold on;
% summary metrics
diagR = diag(Rvv); % ~1 = well resolved, ~0 = constraint-determined
dof_vsv = trace(Rvv); % effective # of independent Vsv parameters
dof_vsv_cumsum = flip(cumsum(flip(diag(Rvv))));
zlay_effective = interp1(dof_vsv_cumsum+[0:length(dof_vsv_cumsum)-1]'*1e-10,z+[0:length(z)-1]'*1e-10,[1:floor(dof_vsv)]);
width_bg = zeros(nlayer,1); % averaging width (Backus-Gilbert-style spread)
for k = 1:nlayer
    a = Rvv(k,:)'; if sum(a)==0, continue; end
    width_bg(k) = sqrt( sum(a.^2 .* (z - z(k)).^2) / sum(a.^2) ) ;   % km Half width (add factor of x2 for full-width)
end
h5(1) = plot(width_bg,z,'-b','linewidth',2);
for ii = 1:length(zlay_effective)
    h5(2) = yline(zlay_effective(ii),'-r','linewidth',2);
end
legend(h5,{'Spread';'Eff. Layer Depth'},'location','southwest');
ylim([0 z_dampbot+100]);
xlabel('Backus-Gilbert Spread (km)');
ylabel('Depth (km)')
set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
