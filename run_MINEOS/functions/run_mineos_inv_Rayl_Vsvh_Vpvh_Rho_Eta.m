function [card_final,cpre,unc] = run_mineos_inv_Rayl_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,card,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,z_dampbot,is_dampcrust,vsv_thresh_dampcrust,nit,nit_recalc_c,nmode,qfile)
% Do linearized inversion of Rayleigh wave phase velocities for Vs, Vp, 
% and Density (rho) using mineos to generate the kernels and to calculate 
% phase velocity.
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
% Thus, although we technically do solve for all 6 parameters,
% we solve primarily for Vsv and simply scale the other 4 parameters (other eta)
% accordingly. Smoothing constraints on Vsv will therefore apply also to 
% Vsh, Vpv, Vph, and rho.
%
% INPUTS (N=number of data; M=number of layers in model)
% cobs - phase velocity observed km/s [N x 1]
% cstd - phase velocity uncertainty km/s [N x 1]
% periods - seconds [N x 1]
% card - starting card: dz (km), vp (km/s), vs (km/s), rho (kg/m^3) [M x 4]
% discs - depth of desired discontinuities in smoothing (km) [any length x 1]
% eps_data - weight for data fit as fraction of G norm [scalar]
% eps_H - weight for norm damping as fraction of G norm [scalar]
% eps_J - weight for 1st derivative smoothing as fraction of G norm [scalar]
% eps_F - weight for 2nd derivative smoothing as fraction of G norm [scalar]
% eps_vpvs - weight for enforcing vp/vs as fraction of G norm [scalar]
% eps_rhovs - weight for enforcing rho/vs as fraction of G norm [scalar]
% z_dampbot - depth below which to damp to starting model (km) [scalar]
% is_dampcrust - strongly damp crust to starting model?
% vsv_thresh_dampcrust - velocity below which is considered "crust" (km/s) [scalar]
% nit - number of iterations [scalar]
% nit_recalc_c - number of iterations after which to recalculate phase velocity and kernels [scalar]
% nmode - mode number [scalar]
% qfile - path to base qfile containing dispersion calculated by a1_run_mineos_check
%
% OUTPUTS
% card_final - final MINEOS card structure
% cpre - phase velocity predicted for final model (km/s) [N x 1]
% unc - formal uncertainties on all model parameters (km/s) [struct: M x 1]
%
% jbrussell - 9/26/2026

cobs = cobs*1000; % km/s -> m/s

eps_large = 1e9; % large weight to force constraint equation

vpv_vsv = card.vpv ./ card.vsv; vpv_vsv(isinf(vpv_vsv))=0;
vph_vsh = card.vph ./ card.vsh; vph_vsh(isinf(vph_vsh))=0;
rho_vsv = card.rho ./ card.vsv; rho_vsv(isinf(rho_vsv))=0;
rho_vsh = card.rho ./ card.vsh; rho_vsh(isinf(rho_vsh))=0;

% % Calculate kernels for G matrix using SURF96
% ifnorm = 0; % for plotting only
% ifplot = 0;
% [dcdvs, dcdvp, dudvs, dudvp, zkern, dcdrho, dudrho] = calc_kernel96(startmod, periods, 'R', ifnorm, ifplot,nmode,fref);
% G = [dcdvs' dcdvp' dcdrho'];

a2_mk_kernels_cv; % !! Comment out "clear" command at the top of a1
disp(FRECH_S)

dr = gradient(FRECH_S(ip).rad); % m
G = [ [FRECH_S(:).vsv]'.*dr' [FRECH_S(:).vsh]'.*dr' [FRECH_S(:).vpv]'.*dr' [FRECH_S(:).vph]'.*dr' [FRECH_S(:).rho]'.*dr' [FRECH_S(:).eta]'.*dr' ];
G = G * 1000; % km/s -> m/s

% Data weighting
% min_pct = 0.005; % minimum error percentage of observed
% I_error_too_small = find(cstd./cobs < min_pct);
% cstd(I_error_too_small) = cobs(I_error_too_small)*min_pct;
W = diag(1./cstd);

% Set up smoothing and damping matrices (applies directly to Vs only)
nlayer = length(card.rad);
dr_mat = repmat(dr,1,length(dr));
% Damping matrix
H00 = eye(nlayer);
h0 = card.vsv;
% first derviative flatness
J00 = build_flatness(nlayer) ./ dr_mat;
j0 = zeros(size(J00,2),1);
% second derivative smoothing
F00 = build_smooth( nlayer ) ./ dr_mat.^2;
f0 = zeros(size(F00,2),1);

% Break constraints at discontinuities
% z_brks = [waterdepth, seddepth, mohodepth];
% z_brks = [ seddepth, mohodepth];
% z_brks = sort([discs; zdisc_Q(:)]);
z_brks = sort(discs);
z = card.z;
J00 = break_constraint(J00, z, z_brks);
F00 = break_constraint(F00, z, z_brks);

% Rescale the kernels
NA=norm(W*G,1);
NR=norm(H00,1);
eps_H0 = eps_H*NA/NR;
NR=norm(J00,1);
eps_J0 = eps_J*NA/NR;
NR=norm(F00,1);
eps_F0 = eps_F*NA/NR;

% Damp towards starting model
ind_dampstart = find(z > z_dampbot);
H00(ind_dampstart,ind_dampstart) = H00(ind_dampstart,ind_dampstart).*linspace(1,1000,length(ind_dampstart));
h0(ind_dampstart) = h0(ind_dampstart).*linspace(1,1000,length(ind_dampstart))';
% Kill water layers
ind_h2o = find(card.vsv==0);
H00(ind_h2o,ind_h2o) = H00(ind_h2o,ind_h2o)*eps_large;
h0(ind_h2o) = h0(ind_h2o)*eps_large;
% Damp crust to starting model
if is_dampcrust
    ind_cr = find(card.vsv < 1000*vsv_thresh_dampcrust);
    H00(ind_cr,ind_cr) = H00(ind_cr,ind_cr)*eps_large;
    h0(ind_cr) = h0(ind_cr)*eps_large;
end

% Add Vsh, Vpv, Vph, rho, and eta dummy zeros
H00=[H00 zeros(size(H00)) zeros(size(H00)) zeros(size(H00)) zeros(size(H00)) zeros(size(H00))];
J00=[J00 zeros(size(J00)) zeros(size(J00)) zeros(size(J00)) zeros(size(J00)) zeros(size(J00))];
F00=[F00 zeros(size(F00)) zeros(size(F00)) zeros(size(F00)) zeros(size(F00)) zeros(size(F00))];

% Add Vp/Vs and Rho/Vs constraints
%                  Vsv               Vsh           Vpv          Vph           Rho          Eta 
VPV_VSV_mat = [-vpv_vsv.*eye(nlayer) zeros(nlayer) eye(nlayer) zeros(nlayer) zeros(nlayer) zeros(nlayer)];
VPH_VSH_mat = [zeros(nlayer) -vph_vsh.*eye(nlayer) zeros(nlayer) eye(nlayer) zeros(nlayer) zeros(nlayer)];
RHO_VSV_mat = [-rho_vsv.*eye(nlayer) zeros(nlayer) zeros(nlayer) zeros(nlayer) eye(nlayer) zeros(nlayer)];
% RHO_VSH_mat = [zeros(nlayer) -rho_vsh.*eye(nlayer) zeros(nlayer) zeros(nlayer) eye(nlayer) zeros(nlayer)];
vpv_vsv_vec = zeros(nlayer,1);
vph_vsh_vec = zeros(nlayer,1);
rho_vsv_vec = zeros(nlayer,1);
% rho_vsh_vec = zeros(nlayer,1);
NR=norm(VPV_VSV_mat,1);
eps_vpvvsv0 = eps_vpvs*NA/NR;
NR=norm(VPH_VSH_mat,1);
eps_vphvsh0 = eps_vpvs*NA/NR;
NR=norm(RHO_VSV_mat,1);
eps_rhovsv0 = eps_rhovs*NA/NR;

% Set Vp of water layer
VPV_h2o_mat = zeros(length(ind_h2o),6*nlayer);
VPH_h2o_mat = zeros(length(ind_h2o),6*nlayer);
for ii = 1:length(ind_h2o)
    VPV_h2o_mat(ii, 2*nlayer + ind_h2o(ii)) = 1;
    VPH_h2o_mat(ii, 3*nlayer + ind_h2o(ii)) = 1;
end
vpv_h2o_vec = 1.5 * ones(length(ind_h2o),1);
vph_h2o_vec = 1.5 * ones(length(ind_h2o),1);
% NR=norm(VPV_h2o_mat,1);

% Set Rho of water layer
RHO_h2o_mat = zeros(length(ind_h2o),6*nlayer);
for ii = 1:length(ind_h2o)
    RHO_h2o_mat(ii, 4*nlayer + ind_h2o(ii)) = 1;
end
rho_h2o_vec = 1.03 * ones(length(ind_h2o),1);
% NR=norm(RHO_h2o_mat,1);

% % Set Vsh = Vsv, isotropic (Vsv - Vsh = 0)
% VSH_VSV_iso_mat = zeros(nlayer,6*nlayer);
% for ii = 1:nlayer
%     VSH_VSV_iso_mat(ii, ii) = 1; % Vsv
%     VSH_VSV_iso_mat(ii, ii + nlayer) = -1; % Vsh
% end
% vsh_vsv_iso_vec = zeros(nlayer,1);

% Damp radial anisotropy toward that of starting model
% Vsh/Vsv = sqrt(xi) --> Vsh = Vsv * sqrt(xi) --> Vsh - Vsv * sqrt(xi) = 0
VSH_VSV_start_mat = zeros(nlayer,6*nlayer);
xi = card.vsh.^2 ./ card.vsv.^2;
xi(isnan(xi)) = 1;
for ii = 1:nlayer
    VSH_VSV_start_mat(ii, ii) = -sqrt(xi(ii)); % Vsv
    VSH_VSV_start_mat(ii, ii + nlayer) = 1; % Vsh
end
vsh_vsv_start_vec = zeros(nlayer,1);
% NR=norm(VSH_VSV_start_mat,1);

% Strongly damp eta to starting value
ETA_damp_mat = zeros(nlayer,6*nlayer);
for ii = 1:nlayer
    ETA_damp_mat(ii, 5*nlayer + ii) = 1;
end
eta_damp_vec = card.eta(:);
% NR=norm(ETA_damp_mat,1);

% combine all constraints
H = [H00*eps_H0; J00*eps_J0; F00*eps_F0; VPV_VSV_mat*eps_vpvvsv0; VPH_VSH_mat*eps_vphvsh0; RHO_VSV_mat*eps_rhovsv0; VPV_h2o_mat*eps_large; VPH_h2o_mat*eps_large; RHO_h2o_mat*eps_large; VSH_VSV_start_mat*eps_large; ETA_damp_mat*eps_large];
h = [h0*eps_H0;  j0*eps_J0;  f0*eps_F0;  vpv_vsv_vec*eps_vpvvsv0; vph_vsh_vec*eps_vphvsh0; rho_vsv_vec*eps_rhovsv0; vpv_h2o_vec*eps_large; vph_h2o_vec*eps_large; rho_h2o_vec*eps_large; vsh_vsv_start_vec*eps_large; eta_damp_vec*eps_large];

% Data vector
% Calculate dispersion for starting model
a1_run_mineos_check % !! Comment out "clear" command at the top of a1
% Get dispersion
[cstart,~,cq_start] = readMINEOS_qfile_per(qfile,periods,nmode);
cstart = cstart(:)*1000; % km/s -> m/s
cpre = cstart(:);

% Least squares inversion
card_pre = card;
vsv_pre = card_pre.vsv;
vsh_pre = card_pre.vsh;
vpv_pre = card_pre.vpv;
vph_pre = card_pre.vph;
rho_pre = card_pre.rho;
eta_pre = card_pre.eta;
m_pre = [vsv_pre; vsh_pre; vpv_pre; vph_pre; rho_pre; eta_pre];
clrs = jet(nit);
isfigure = 1;
clear vs
for ii = 1:nit
    
    % reformulate inverse problem so that constraints apply directly to model
    % G(m-m0) = dobs-d0
    % Gm = G*m0 + dobs-d0
    %    = Ddobs
    dc = cobs - cpre;
    d = dc + G*m_pre;
    
    % least squares
    F = [W*G*eps_data; H];
    f = [W*d*eps_data; h];
    m = (F'*F)\F'*f;
    vsv = m(1:nlayer);
    vsh = m(  nlayer+1:2*nlayer);
    vpv = m(2*nlayer+1:3*nlayer);
    vph = m(3*nlayer+1:4*nlayer);
    rho = m(4*nlayer+1:5*nlayer);
    eta = m(5*nlayer+1:6*nlayer);
    vsv(vsv<eps) = 0;
    vsh(vsh<eps) = 0;
    m = [vsv; vsh; vpv; vph; rho; eta];
    
    % update data vector
    dm = m-m_pre;
    dc_update = G * dm;
    cpre = cpre + dc_update;
    
    % update model
    m_pre = m_pre + dm;
    vsv_pre = m_pre(1:nlayer);
    vsh_pre = m_pre(  nlayer+1:2*nlayer);
    vpv_pre = m_pre(2*nlayer+1:3*nlayer);
    vph_pre = m_pre(3*nlayer+1:4*nlayer);
    rho_pre = m_pre(4*nlayer+1:5*nlayer);
    eta_pre = m_pre(5*nlayer+1:6*nlayer);
    card_pre.vsv = vsv_pre;
    card_pre.vsh = vsh_pre;
    card_pre.vpv = vpv_pre;
    card_pre.vph = vph_pre;
    card_pre.rho = rho_pre;
    card_pre.eta = eta_pre;
    
    % Model uncertainties
    m_std = diag(inv(F'*F)).^(1/2);
    unc.vsv_std = m_std(1:nlayer);
    unc.vsh_std = m_std(  nlayer+1:2*nlayer);
    unc.vpv_std = m_std(2*nlayer+1:3*nlayer);
    unc.vph_std = m_std(3*nlayer+1:4*nlayer);
    unc.rho_std = m_std(4*nlayer+1:5*nlayer);
    unc.eta_std = m_std(5*nlayer+1:6*nlayer);
    
    if mod(ii,nit_recalc_c)==0
        error('HAVENT YET IMPLEMENTED THIS FOR MINEOS!')
        cpre = dispR_surf96(periods,premod,nmode,'C',fref);
        dc = cobs - cpre;
        
        ifnorm = 0; % for plotting only
        ifplot = 0;
        [dcdvs, dcdvp, dudvs, dudvp, zkern, dcdrho, dudrho] = calc_kernel96(startmod, periods, 'R', ifnorm, ifplot,nmode,fref);
        G = [dcdvs' dcdvp' dcdrho'];
    end
    
    if isfigure
        if ii == 1 
            figure(1000); clf; set(gcf,'position',[370   372   967   580]);
            
            subplot(2,2,[1 3]); box on; hold on;
            plot(card.vsv,card.z,'-k','linewidth',4);
            ylim([0 z_dampbot+100]);
            
            subplot(2,2,2); box on; hold on;
            plot(periods,cstart,'-ok','linewidth',4);
        end
        subplot(2,2,[1 3]);
        plot(card_pre.vsv,card_pre.z,'-k','color',clrs(ii,:),'linewidth',2);
        xlabel('Vs (km/s)');
        ylabel('Depth (km)');
        title('Starting Model');
        set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
        
        subplot(2,2,2);
        plot(periods,cpre,'-o','color',clrs(ii,:),'linewidth',2);
        errorbar(periods,cobs,2*cstd,'or','linewidth',4);
        xlabel('Period');
        ylabel('Phase Velocity');
        set(gca,'FontSize',18,'linewidth',1.5);
        % pause;
        drawnow
    end
end

card_final = card_pre;
cpre = cpre/1000; % m/s -> km/s
% cpre = dispR_surf96(periods,finalmod,nmode,'C',fref);

end

