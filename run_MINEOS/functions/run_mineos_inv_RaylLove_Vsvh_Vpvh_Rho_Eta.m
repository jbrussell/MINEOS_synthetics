function [card_final,cpre,unc,Kfinal] = run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,wavetype,card,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,zmax,nit,nit_recalc_c,nmode,cardid_work,varargin)
% [card_final,cpre,unc,Kfinal] = run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(..., nmode, cardid_work, isfigure, vsv_thresh_xi, nit_recalc_G, use_Q, fref)
%
% Linearized inversion of Rayleigh and/or Love wave phase velocities for a
% transversely isotropic (radially anisotropic) model: Vsv, Vsh, Vpv, Vph,
% density (rho) and eta, using MINEOS to calculate the phase velocities
% and kernels (run_mineos_phv_kernels.m). This is the MINEOS version of
% run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m in SURF96_Rayl_1Dinv
% (https://github.com/jbrussell/SURF96_Rayl_1Dinv), with the same
% constraints and weights. The model parameters are the values at the card
% nodes shallower than zmax; deeper nodes are held fixed.
%
% Constraint equations:
% - Norm damping of Vsv and Vsh toward the starting model (increasing
%   linearly from 1 at z_dampbot to 1000 at zmax)
% - First derivative (flatness) smoothing of Vsv and Vsh
% - Second derivative (roughness) smoothing of Vsv and Vsh
% - Damping Vpv/Vsv, Vph/Vsh and Rho/Vsv toward the starting-model ratios
% - Damping xi = (Vsh/Vsv)^2 toward the starting model (weight eps_xi):
%   use a large eps_xi for Rayleigh-only or Love-only data
% - Damping eta toward the starting model (weight eps_eta)
% - (optional) xi held at the starting model where 0 < Vsv < vsv_thresh_xi
% - Water layer(s) fixed to the starting model
% The smoothing is broken at the discontinuities (pairs of card nodes at
% the same radius) at the depths in discs.
%
% INPUTS (N = number of data)
% cobs - phase velocity observed km/s [N x 1]
% cstd - phase velocity uncertainty km/s [N x 1]
% periods - seconds [N x 1]
% wavetype - 'R' or 'L' for each datum [N x 1 char], or a single character
% card - starting MINEOS card structure (read_model_card.m; m/s, kg/m^3).
%        card.fname must point to the card file; the .qmod file with the
%        same name is used for the recalculations.
% discs - depths (km) of the discontinuities at which to break the
%         smoothing (any length; [] = smooth across all of them)
% eps_data - weight for data fit [scalar]
% eps_H - weight for norm damping as fraction of G norm [scalar]
% eps_J - weight for 1st derivative smoothing as fraction of G norm [scalar]
% eps_F - weight for 2nd derivative smoothing as fraction of G norm [scalar]
% eps_vpvs - weight for enforcing Vpv/Vsv and Vph/Vsh as fraction of G norm
% eps_rhovs - weight for enforcing rho/Vsv as fraction of G norm
% eps_xi - weight for damping xi to the starting model as fraction of G norm
% eps_eta - weight for damping eta to the starting model as fraction of G norm
% z_dampbot - depth below which to damp to the starting model (km)
% zmax - maximum depth of the inverted region (km)
% nit - number of iterations
% nit_recalc_c - number of iterations after which to recalculate the phase
%                velocities with MINEOS (otherwise the linearized prediction
%                G*dm is used)
% nmode - mode branch (0 = fundamental; the kernels are set up for 0 only)
% cardid_work - name of the card written to ./CARDS/ for the MINEOS
%               recalculations (e.g. 'c2_inv')
% isfigure - plot the progress of the inversion (OPTIONAL, default 1)
% vsv_thresh_xi - hold xi at the starting model where 0 < Vsv < vsv_thresh_xi
%                 (km/s) in the starting model (OPTIONAL, default 0 = off)
% nit_recalc_G - number of iterations after which to recalculate the
%                kernels (OPTIONAL, default nit_recalc_c)
% use_Q - 1 = include attenuation (see below) (OPTIONAL, default 0)
% fref - reference frequency (Hz) of the Q model (OPTIONAL, default 1)
%
% ATTENUATION: with use_Q = 1, the phase velocities and kernels include the
% causal (physical dispersion) correction for the card's Qmu and Qkappa,
% with the card velocities defined at fref (mineos_apply_Q.m; the same
% first-order correction as the CPS programs and the SURF96 version). Q is
% not inverted for: Qmu and Qkappa are held at their starting values. With
% use_Q = 0 the forward problem is elastic.
%
% OUTPUTS
% card_final - final MINEOS card structure (m/s, kg/m^3)
% cpre - phase velocity of the final model calculated with MINEOS (km/s)
% unc - formal uncertainties at the inverted nodes (km/s, g/cm^3):
%       vsv_std, vsh_std, vpv_std, vph_std, rho_std, eta_std, the node
%       indices I (into the card) and depths z, the model covariance Cm and
%       resolution R [6M x 6M] (order Vsv, Vsh, Vpv, Vph, Rho, Eta) and the
%       data resolution Ndat [N x N], and the phase velocities of the
%       starting model cstart (km/s) [N x 1]
% Kfinal - kernels (causal if use_Q) of the final model at the inverted nodes, in km/s per
%       km/s (g/cm^3 for rho) integrated over each node's depth interval:
%       Kfinal.(wavetype).dcdvsv etc. [M x Nper], Kfinal.z, Kfinal.periods
%
% UNITS: MINEOS kernels K satisfy dc (km/s) = sum K .* dr (m) .* dm, with dm
% in the card units (m/s, kg/m^3; eta dimensionless). Here the model is
% in km/s and g/cm^3 (as in the SURF96 version), so G = K.*dr*1000 for the
% velocities and density and G = K.*dr for eta.
%
% jbrussell - 11/2026

isfigure = 1;
if ~isempty(varargin) && ~isempty(varargin{1})
    isfigure = varargin{1};
end
vsv_thresh_xi = 0;
if length(varargin) >= 2 && ~isempty(varargin{2})
    vsv_thresh_xi = varargin{2};
end
nit_recalc_G = nit_recalc_c;
if length(varargin) >= 3 && ~isempty(varargin{3})
    nit_recalc_G = varargin{3};
end
use_Q = 0;
if length(varargin) >= 4 && ~isempty(varargin{4})
    use_Q = varargin{4};
end
fref = 1;
if length(varargin) >= 5 && ~isempty(varargin{5})
    fref = varargin{5};
end
if ~use_Q
    fref = []; % elastic
end

cobs = cobs(:); cstd = cstd(:); periods = periods(:);
N = length(cobs);
wavetype = upper(wavetype(:));
if length(wavetype) == 1
    wavetype = repmat(wavetype, N, 1);
end
if length(wavetype) ~= N
    error('wavetype must have one entry per datum (or be a single character)');
end

eps_large = 1e4; % large weight (relative to the G norm) to force constraint equations

% Files used for the MINEOS recalculations
qmod_start = regexprep(card.fname, '\.card$', '.qmod');
qmod_work = ['./CARDS/', cardid_work, '.qmod'];
if ~strcmp(qmod_start, qmod_work)
    copyfile(qmod_start, qmod_work);
end

% Inverted nodes: card nodes shallower than zmax
I = find(card.z <= zmax);
nlayer = length(I);
z = card.z(I);
dr = gradient(card.rad); % node weights (m) for the kernel integrals
dr = dr(I);

% Kernels and phase velocities of the starting model
[G, cstart] = build_G(card, I, dr, periods, wavetype, nmode, cardid_work, 1, fref);

% Data weighting
W = diag(1./cstd);

% Starting model (km/s, g/cm^3)
m0 = card2vec(card, I);
vsv0 = m0(1:nlayer); vsh0 = m0(nlayer+1:2*nlayer);
vpv0 = m0(2*nlayer+1:3*nlayer); vph0 = m0(3*nlayer+1:4*nlayer);
rho0 = m0(4*nlayer+1:5*nlayer); eta0 = m0(5*nlayer+1:6*nlayer);
vpv_vsv = vpv0 ./ vsv0;
vph_vsh = vph0 ./ vsh0;
rho_vsv = rho0 ./ vsv0;
xi0 = vsh0.^2 ./ vsv0.^2;
ind_h2o = find(vsv0 == 0);
ind_solid = find(vsv0 > 0);
xi0(ind_h2o) = 1;

% Node spacing for the smoothing (ignoring the zero spacing across
% discontinuities)
dz = node_spacing(z);
dz_mat = repmat(dz,1,nlayer);
% Damping matrix
H00 = eye(nlayer);
% first derivative flatness
J00 = build_flatness(nlayer) ./ dz_mat;
% second derivative smoothing
F00 = build_smooth(nlayer) ./ dz_mat.^2;

% Break the smoothing at discontinuities (and never smooth across a pair
% of nodes at the same depth)
J00 = break_constraint_nodes(J00, z, discs);
F00 = break_constraint_nodes(F00, z, discs);

% Rescale the constraint weights relative to the data kernels
NA = norm(W*G,1);
eps_H0 = eps_H*NA/norm(H00,1);
eps_J0 = eps_J*NA/norm(J00,1);
eps_F0 = eps_F*NA/norm(F00,1);

% Damp towards starting model below z_dampbot
hv = vsv0; hh = vsh0;
ind_dampstart = find(z > z_dampbot);
ramp = 1 + 999*(z(ind_dampstart) - z_dampbot)/(zmax - z_dampbot);
H00(ind_dampstart,ind_dampstart) = H00(ind_dampstart,ind_dampstart).*ramp(:)';
hv(ind_dampstart) = hv(ind_dampstart).*ramp;
hh(ind_dampstart) = hh(ind_dampstart).*ramp;

% Column blocks of the model vector m = [Vsv; Vsh; Vpv; Vph; Rho; Eta]
Z = zeros(nlayer);
Iden = eye(nlayer);
blk = @(k) [repmat(Z,1,k-1), Iden, repmat(Z,1,6-k)];

% Damping and smoothing on Vsv (block 1) and Vsh (block 2)
Hd = [H00*blk(1); H00*blk(2)];  hd = [hv; hh];
Jd = [J00*blk(1); J00*blk(2)];  jd = zeros(2*size(J00,1),1);
Fd = [F00*blk(1); F00*blk(2)];  fd = zeros(2*size(F00,1),1);

% Ratio constraints (solid nodes only)
%       Vpv - (Vpv/Vsv)0 Vsv = 0,  Vph - (Vph/Vsh)0 Vsh = 0,  Rho - (Rho/Vsv)0 Vsv = 0
S = Iden(ind_solid,:);
VPV_mat = S*blk(3) - S*diag(fixnan(vpv_vsv))*blk(1);
VPH_mat = S*blk(4) - S*diag(fixnan(vph_vsh))*blk(2);
RHO_mat = S*blk(5) - S*diag(fixnan(rho_vsv))*blk(1);
eps_vpvs0 = eps_vpvs*NA/norm(VPV_mat,1);
eps_rhovs0 = eps_rhovs*NA/norm(RHO_mat,1);

% xi constraint (solid nodes):  Vsh - sqrt(xi0) Vsv = 0
XI_mat = S*blk(2) - S*diag(sqrt(xi0))*blk(1);
eps_xi0 = eps_xi*NA/norm(XI_mat,1);
% hold xi at the starting model in the crust/sediments (Vsv < vsv_thresh_xi)
w_xi = ones(length(ind_solid),1);
w_xi(vsv0(ind_solid) < vsv_thresh_xi) = eps_large*NA/eps_xi0;
XI_mat = diag(w_xi)*XI_mat;

% eta damping toward the starting model
ETA_mat = blk(6);
eps_eta0 = eps_eta*NA/norm(ETA_mat,1);

% Water layer(s): all parameters fixed to the starting model
W_h2o = zeros(6*length(ind_h2o), 6*nlayer); w_h2o = zeros(6*length(ind_h2o),1);
irow = 0;
for k = 1:6
    for ii = ind_h2o(:)'
        irow = irow + 1;
        W_h2o(irow, (k-1)*nlayer + ii) = 1;
        w_h2o(irow) = m0((k-1)*nlayer + ii);
    end
end

% combine all constraints
H = [Hd*eps_H0; Jd*eps_J0; Fd*eps_F0; VPV_mat*eps_vpvs0; VPH_mat*eps_vpvs0; RHO_mat*eps_rhovs0; ...
     XI_mat*eps_xi0; ETA_mat*eps_eta0; W_h2o*eps_large*NA];
h = [hd*eps_H0; jd*eps_J0; fd*eps_F0; zeros(3*length(ind_solid),1); ...
     zeros(length(ind_solid),1); eta0*eps_eta0; w_h2o*eps_large*NA];

% Data vector
if any(isnan(cstart))
    error('run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta: %d phase velocities of the starting model could not be calculated', sum(isnan(cstart)));
end
cpre = cstart;

% Least squares inversion
card_pre = card;
m_pre = m0;
clrs = jet(nit);
for ii = 1:nit

    % reformulate inverse problem so that constraints apply directly to model
    % G(m-m0) = dobs-d0  -->  Gm = G*m0 + dobs-d0
    dc = cobs - cpre;
    d = dc + G*m_pre;

    % least squares
    F = [W*G*eps_data; H];
    f = [W*d*eps_data; h];
    m = (F'*F)\(F'*f);
    m(1:2*nlayer) = max(m(1:2*nlayer), 0); % no negative shear velocities

    % update data vector
    dm = m - m_pre;
    cpre = cpre + G*dm;

    % update model
    m_pre = m;
    card_pre = vec2card(m_pre, card_pre, I);

    % Model uncertainties and resolution
    FtF = F'*F;
    Ginv = eps_data^2 * (FtF \ (G'*(W'*W))); % generalized inverse [6M x N]
    unc.R = Ginv * G;                        % model resolution [6M x 6M]
    unc.Ndat = G * Ginv;                     % data resolution [N x N]
    unc.Cm = Ginv * diag(cstd.^2) * Ginv';   % model covariance from data errors
    m_std = sqrt(diag(unc.Cm));
    unc.vsv_std = m_std(1:nlayer);
    unc.vsh_std = m_std(nlayer+1:2*nlayer);
    unc.vpv_std = m_std(2*nlayer+1:3*nlayer);
    unc.vph_std = m_std(3*nlayer+1:4*nlayer);
    unc.rho_std = m_std(4*nlayer+1:5*nlayer);
    unc.eta_std = m_std(5*nlayer+1:6*nlayer);
    unc.I = I;
    unc.z = z;
    unc.cstart = cstart;

    do_c = mod(ii,nit_recalc_c)==0 && ii < nit; % recalculate phase velocities
    do_G = mod(ii,nit_recalc_G)==0 && ii < nit; % recalculate kernels
    if do_G
        % kernels and phase velocities (one MINEOS calculation)
        [G, cnew] = build_G(card_pre, I, dr, periods, wavetype, nmode, cardid_work, 1, fref);
    elseif do_c
        % phase velocities only
        [~, cnew] = build_G(card_pre, I, dr, periods, wavetype, nmode, cardid_work, 0, fref);
    end
    if do_c || do_G
        inan = isnan(cnew);
        if any(inan)
            warning('run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta: %d phase velocities could not be calculated; keeping the linearized prediction there', sum(inan));
            cnew(inan) = cpre(inan);
        end
        cpre = cnew;
    end

    if isfigure
        figure(2); % (the MINEOS kernel calculation draws its own figures)
        if ii == 1
            clf; set(gcf,'position',[370 372 1100 580]);
            subplot(2,3,[1 4]); box on; hold on;
            plot(vsv0, z, '-k', 'linewidth', 4);
            subplot(2,3,[2 5]); box on; hold on;
            plot(vsh0, z, '-k', 'linewidth', 4);
        end
        subplot(2,3,[1 4]);
        plot(card_pre.vsv(I)/1000, z, '-', 'color', clrs(ii,:), 'linewidth', 2);
        xlabel('Vsv (km/s)'); ylabel('Depth (km)');
        title('Iterations (black: start)');
        set(gca,'FontSize',15,'linewidth',1.5,'ydir','reverse');
        subplot(2,3,[2 5]);
        plot(card_pre.vsh(I)/1000, z, '-', 'color', clrs(ii,:), 'linewidth', 2);
        xlabel('Vsh (km/s)'); ylabel('Depth (km)');
        set(gca,'FontSize',15,'linewidth',1.5,'ydir','reverse');
        for iw = 1:2
            wt = 'RL'; I_w = wavetype == wt(iw);
            if ~any(I_w); continue; end
            subplot(2,3,3*iw); box on; hold on;
            if ii == 1
                plot(periods(I_w),cstart(I_w),'-ok','linewidth',3);
            end
            plot(periods(I_w),cpre(I_w),'-o','color',clrs(ii,:),'linewidth',2);
            he = errorbar(periods(I_w),cobs(I_w),2*cstd(I_w),'or'); set(he,'linewidth',2);
            xlabel('Period (s)'); ylabel([wt(iw),' phase velocity (km/s)']);
            set(gca,'FontSize',15,'linewidth',1.5);
        end
        drawnow
    end
end

card_final = card_pre;
% Phase velocities (and kernels) of the final model
[~, cpre, Kfinal] = build_G(card_final, I, dr, periods, wavetype, nmode, cardid_work, nargout >= 4, fref);
Kfinal.z = z;

end

%% ------------------------------------------------------------------------
function [G, c, Kout] = build_G(card, I, dr, periods, wavetype, nmode, cardid_work, do_kernels, fref)
% Write the card, run MINEOS and build the kernel matrix [N x 6M] for
% m = [Vsv; Vsh; Vpv; Vph; Rho; Eta] (km/s, g/cm^3) at the nodes I, and the
% phase velocities c [N x 1] (km/s). fref = [] for elastic, otherwise the
% reference frequency of the Q model (causal phase velocities and kernels)
write_MINEOS_mod(card, ['./CARDS/', cardid_work, '.card']);
N = length(periods);
nlayer = length(I);
nnode = length(card.rad);
G = zeros(N, 6*nlayer);
c = nan(N,1);
Kout = struct();
flds = {'vsv','vsh','vpv','vph','rho','eta'};
scale = [1000 1000 1000 1000 1000 1]; % card units -> km/s, g/cm^3
for wt = 'RL'
    I_w = find(wavetype == wt);
    if isempty(I_w); continue; end
    [c(I_w), K] = run_mineos_phv_kernels(cardid_work, wt, periods(I_w), nmode, do_kernels, [], fref);
    if ~do_kernels; continue; end
    % MINEOS kernels cover the top nK nodes of the card (all nodes for
    % spheroidal modes; the mantle and crust for toroidal modes)
    nK = length(K.rad);
    iK = I - (nnode - nK);
    if any(iK < 1)
        error('run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta: kernels do not cover the inverted nodes');
    end
    for k = 1:6
        Kk = K.(flds{k})(iK,:) .* dr(:) * scale(k); % [M x Nper]
        Kk(isnan(Kk)) = 0;
        G(I_w, (k-1)*nlayer+(1:nlayer)) = Kk';
        Kout.(wt).(['dcd',flds{k}]) = Kk;
    end
    Kout.(wt).periods = periods(I_w);
end
end

function m = card2vec(card, I)
% card (m/s, kg/m^3) -> [Vsv; Vsh; Vpv; Vph; Rho; Eta] (km/s, g/cm^3) at nodes I
m = [card.vsv(I); card.vsh(I); card.vpv(I); card.vph(I); card.rho(I); card.eta(I)];
n = length(I);
m(1:5*n) = m(1:5*n)/1000;
end

function card = vec2card(m, card, I)
n = length(I);
card.vsv(I) = m(1:n)*1000;
card.vsh(I) = m(n+1:2*n)*1000;
card.vpv(I) = m(2*n+1:3*n)*1000;
card.vph(I) = m(3*n+1:4*n)*1000;
card.rho(I) = m(4*n+1:5*n)*1000;
card.eta(I) = m(5*n+1:6*n);
end

function dz = node_spacing(z)
% node spacing for the smoothing operators, using the neighbor at a
% different depth where two nodes share a depth (discontinuity)
n = length(z);
dz = zeros(n,1);
for i = 1:n
    nb = [];
    if i > 1 && z(i-1) ~= z(i); nb(end+1) = abs(z(i)-z(i-1)); end %#ok<AGROW>
    if i < n && z(i+1) ~= z(i); nb(end+1) = abs(z(i+1)-z(i)); end %#ok<AGROW>
    if isempty(nb); nb = 1; end
    dz(i) = mean(nb);
end
end

function X = break_constraint_nodes(X, z, discs)
% Break the smoothing at the discontinuities in discs: remove the equations
% that involve nodes on both sides of the pair of nodes (same depth) at
% each of these depths. Pairs at other depths are smoothed across.
ipair = find(diff(z) == 0); % nodes ipair and ipair+1 share a depth
for ii = ipair(:)'
    if isempty(discs) || ~any(abs(discs - z(ii)) < 1e-3)
        continue
    end
    cols_a = 1:ii; cols_b = ii+1:size(X,2);
    cross = any(X(:,cols_a)~=0,2) & any(X(:,cols_b)~=0,2);
    X(cross,:) = 0;
end
end

function x = fixnan(x)
x(~isfinite(x)) = 0;
end
