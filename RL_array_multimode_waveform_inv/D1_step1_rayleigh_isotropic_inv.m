% D1_STEP1_RAYLEIGH_ISOTROPIC_INV  Step 1 of the two-step inversion (driver).
%
% Linearized MINEOS inversion of fundamental-mode Rayleigh wave phase
% velocities (20-100 s) for an ISOTROPIC Vs model (xi = 1), with the
% inversion engine of ../run_MINEOS/c2_RLphv_inv_vsvh_vpvh_rho_eta.m
% (run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m):
%   - data: MINEOS phase velocities of the true model (PA5), branch 0
%   - starting model (opt.s1_start_model):
%       'prem'        isotropic PREM (Voigt average) below the PA5 Moho, PA5
%                     water/sediments/crust above it; eta = 1
%       'true_smooth' isotropic PA5 (Voigt average) with its lid/LVZ
%                     structure smoothed away (functions/smooth_lid_lvz_card.m):
%                     no lid, no sharp LVZ, PA5 below opt.s1_smooth_blend(2)
%                     and above the Moho; Vpv/Vsv, Vph/Vsh, rho/Vsv, eta and Q
%                     of PA5 at every node
%   - xi held at 1 (eps_xi large, starting xi = 1), eta held at the starting
%     model
%   - Vpv/Vsv, Vph/Vsh and rho/Vsv held at the starting-model ratios (with
%     'true_smooth' these ratios and eta are imposed exactly on the final
%     model, so the step-2 reference has the Vp/Vs, rho/Vs and eta of PA5)
%   - MINEOS phase velocities and kernels recalculated after every iteration
% The result (CARDS/<opt.s1_final>.card) is the reference model of step 2
% (isotropic in Vs).
%
% Run from this folder. Needs the MINEOS binaries (../FORTRAN/bin).
% Output: opt.s1_file (step1_result.mat), figs_step1/.
%
% Next: D2_step2_prep_mineos_idagrn.m

close all;
parameters_RL;
mdir = opt.mdir;
here = pwd;
cleanupObj = onCleanup(@() cd(here));
if ~exist('figs_step1', 'dir'), mkdir('figs_step1'); end

%% models
card_prem = make_ref_card(mdir, opt.premref, opt.premid, opt.trueid, opt.zmoho);  % anisotropic PREM + PA5 crust
card_true = read_model_card(fullfile(mdir, 'CARDS', [opt.trueid '.card']));
smooth_ = strcmp(opt.s1_start_model, 'true_smooth');
if smooth_
    card_start = smooth_lid_lvz_card(card_true, opt.zmoho, opt.s1_smooth_sig, opt.s1_smooth_blend);
    qid = opt.ref_qmod; start_lbl = 'start (smoothed isotropic PA5)';
elseif strcmp(opt.s1_start_model, 'prem')
    card_start = isotropic_card(card_prem);
    qid = opt.ref_qmod; start_lbl = 'start (isotropic PREM)';
else
    error('opt.s1_start_model must be ''prem'' or ''true_smooth''');
end
fstart = fullfile(mdir, 'CARDS', [opt.s1_start '.card']);
write_MINEOS_mod(card_start, fstart);
copyfile(fullfile(mdir, 'CARDS', [qid '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_start '.qmod']));

%% data: fundamental-mode Rayleigh phase velocities of the true model
periods = opt.s1_periods(:);
cd(mdir);
cobs = run_mineos_phv_kernels(opt.trueid, 'R', periods, 0, 0);
cstd = opt.s1_cstd * cobs;
card_start = read_model_card(['./CARDS/' opt.s1_start '.card']);   % fname relative to run_MINEOS
% break the smoothing only at the crustal discontinuities (water bottom,
% sediments, Moho); the PREM discontinuities (80, 220 km) are smoothed over
zc = card_start.z;
discs = unique(zc(find(diff(zc) == 0 & zc(1:end-1) <= opt.zmoho + 1e-6)));

%% inversion
e = opt.s1_eps;
[card_final, cpre, unc] = run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs, cstd, periods, 'R', ...
    card_start, discs, e.data, e.H, e.J, e.F, e.vpvs, e.rhovs, e.xi, e.eta, ...
    opt.s1_z_dampbot, opt.s1_zmax, opt.s1_nit, 1, 0, opt.s1_work, 1, 4, 1);
cstart = unc.cstart;
cd(here);

%% step-1 model: the reference of step 2
if smooth_
    % impose the starting (= PA5) Vp/Vs, rho/Vs and eta exactly (the
    % inversion holds them with finite weights only); Vsv = Vsh (xi = 1)
    so_ = card_start.vsv > 0;
    card_final.vsh(so_) = card_final.vsv(so_);
    card_final.vpv(so_) = card_start.vpv(so_)./card_start.vsv(so_).*card_final.vsv(so_);
    card_final.vph(so_) = card_start.vph(so_)./card_start.vsh(so_).*card_final.vsh(so_);
    card_final.rho(so_) = card_start.rho(so_)./card_start.vsv(so_).*card_final.vsv(so_);
    card_final.eta = card_start.eta;
    cd(mdir);                                  % phase velocities of the final card
    write_MINEOS_mod(card_final, fullfile('CARDS', [opt.s1_work '.card']));
    copyfile(fullfile('CARDS', [qid '.qmod']), fullfile('CARDS', [opt.s1_work '.qmod']));
    cpre = run_mineos_phv_kernels(opt.s1_work, 'R', periods, 0, 0);
    cd(here);
end
write_MINEOS_mod(card_final, fullfile(mdir, 'CARDS', [opt.s1_final '.card']));
copyfile(fullfile(mdir, 'CARDS', [qid '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_final '.qmod']));
fprintf('Rayleigh fundamental mode, reduced chi^2: start %.1f, final %.2f\n', ...
    mean(((cobs - cstart)./cstd).^2), mean(((cobs - cpre)./cstd).^2));
save(opt.s1_file, 'card_prem', 'card_start', 'card_final', 'card_true', 'cobs', 'cstd', ...
     'cstart', 'cpre', 'periods', 'discs', 'opt');

%% figure
zw = max(card_true.z(card_true.vsv == 0 & card_true.z < 20));
zp = (0:1:400)';
cp = @(c, f) card_at_depth(c, f, zp + zw)/1000;
fig = figure('position', [50 50 1300 520]);
subplot(1, 3, 1); hold on; box on;
h1 = plot(cp(card_true, 'vsv'), zp, 'r-', 'linewidth', 2);
h2 = plot(cp(card_true, 'vsh'), zp, 'r--', 'linewidth', 1);
h3 = plot(cp(card_prem, 'vsv'), zp, '-', 'color', [0 0.6 0], 'linewidth', 1);
h4 = plot(cp(card_start, 'vsv'), zp, 'g-', 'linewidth', 1.5);
h5 = plot(cp(card_final, 'vsv'), zp, 'b-', 'linewidth', 2);
set(gca, 'ydir', 'reverse'); ylim([0 400]); xlim([3.9 5]); xlabel('V_S (km/s)'); ylabel('Depth below seafloor (km)');
legend([h1 h2 h3 h4 h5], {'PA5 V_{SV}', 'PA5 V_{SH}', 'PREM V_{SV}', start_lbl, 'step 1 (isotropic)'}, ...
    'location', 'southwest', 'fontsize', 8);
title('Step 1: isotropic V_S from Rayleigh n = 0');
subplot(1, 3, 2); hold on; box on;
plot(cp(card_true, 'vpv'), zp, 'r-', cp(card_start, 'vpv'), zp, 'g-', cp(card_final, 'vpv'), zp, 'b-');
set(gca, 'ydir', 'reverse'); ylim([0 400]); xlabel('V_{PV} (km/s)'); title('V_P (scaled with V_S)');
subplot(1, 3, 3); hold on; box on;
he = errorbar(periods, cobs, 2*cstd, 'ko'); set(he, 'markerfacecolor', 'k');
plot(periods, cstart, 'g-o', periods, cpre, 'b-o', 'linewidth', 1.5);
xlabel('Period (s)'); ylabel('Rayleigh n = 0 phase velocity (km/s)');
legend({'PA5 (data)', 'start', 'step 1'}, 'location', 'southeast');
print(fig, '-dpng', '-r130', fullfile('figs_step1', 'fig_step1.png'));
