% D1_STEP1_RAYLEIGH_ISOTROPIC_INV  Step 1 of the two-step inversion (driver).
%
% Linearized MINEOS inversion of fundamental-mode Rayleigh wave phase
% velocities (opt.s1_periods, 20-150 s) for an ISOTROPIC Vs model (xi = 1), with the
% inversion engine of ../run_MINEOS/c2_RLphv_inv_vsvh_vpvh_rho_eta.m
% (run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m):
%   - data: MINEOS phase velocities of the true model (PA5), branch 0
%   - starting model (opt.s1_start_model):
%       'pa5_iso' (the only one) isotropic PA5 without its lid/LVZ: between the
%                 Moho and opt.s1_zlin (250 km below sea level) Vs is a
%                 straight line between the Voigt Vs of PA5 at those depths;
%                 Vpv/Vsv, Vph/Vsh, rho/Vsv, eta and Q of PA5 at every node;
%                 PA5 elsewhere (functions/pa5_iso_ref_card.m)
%   - xi held at 1 (eps_xi large, starting xi = 1), eta held at the start
%   - Vpv/Vsv, Vph/Vsh and rho/Vsv held at the starting-model ratios; with
%     'pa5_iso' these ratios and eta are imposed exactly on the result, so
%     the step-2 reference has PA5's Vp/Vs, rho/Vs and eta
%   - MINEOS phase velocities and kernels recalculated after every iteration
% The result (CARDS/<opt.s1_final>.card) is the reference model of step 2:
% an isotropic model with a lid and a low-velocity zone.
%
% Run from this folder. Needs the MINEOS binaries (../FORTRAN/bin).
% Output: opt.s1_file (step1_result.mat), figs_step1/.
%
% Next: D2_step2_prep_mineos_idagrn.m

parameters_RL;
if fig_close_mode, close all; end
mdir = opt.mdir;
here = pwd;
cleanupObj = onCleanup(@() cd(here));
if ~exist('figs_step1', 'dir'), mkdir('figs_step1'); end

%% models
card_true = read_model_card(fullfile(mdir, 'CARDS', [opt.trueid '.card']));
card_start = pa5_iso_ref_card(card_true, opt.zmoho, opt.s1_zlin);
start_lbl = 'start (isotropic PA5, no lid/LVZ)';
fstart = fullfile(mdir, 'CARDS', [opt.s1_start '.card']);
write_MINEOS_mod(card_start, fstart);
copyfile(fullfile(mdir, 'CARDS', [opt.ref_qmod '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_start '.qmod']));

%% data: fundamental-mode Rayleigh phase velocities of the true model
periods = opt.s1_periods(:);
cd(mdir);
cobs = run_mineos_phv_kernels(opt.trueid, 'R', periods, 0, 0);
cstd = opt.s1_cstd * cobs;
if optget(opt, 's1b_noise', true)             % Gaussian data noise (opt.s1b_noise; as D1_step1_Bayesian_disp.m)
    rs_ = rng; rng(optget(opt, 's1b_noise_seed', 1));        % (rng, not the legacy randn('state'): MATLAB forbids rng after it)
    cobs = cobs + cstd .* randn(size(cobs)); rng(rs_);
end
card_start = read_model_card(['./CARDS/' opt.s1_start '.card']);   % fname relative to run_MINEOS
% break the smoothing only at the crustal discontinuities (water bottom,
% sediments, Moho); any deeper discontinuities of the start are smoothed over
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
if pa5iso
    % impose the starting (= PA5) Vp/Vs, rho/Vs and eta exactly (the
    % inversion holds them with finite weights only); Vsh = Vsv (xi = 1)
    so_ = card_start.vsv > 0;
    card_final.vsh(so_) = card_final.vsv(so_);
    card_final.vpv(so_) = card_start.vpv(so_)./card_start.vsv(so_).*card_final.vsv(so_);
    card_final.vph(so_) = card_start.vph(so_)./card_start.vsh(so_).*card_final.vsh(so_);
    card_final.rho(so_) = card_start.rho(so_)./card_start.vsv(so_).*card_final.vsv(so_);
    card_final.eta = card_start.eta;
    cd(mdir);                                  % phase velocities of the final card
    write_MINEOS_mod(card_final, fullfile('CARDS', [opt.s1_work '.card']));
    copyfile(fullfile('CARDS', [opt.ref_qmod '.qmod']), fullfile('CARDS', [opt.s1_work '.qmod']));
    cpre = run_mineos_phv_kernels(opt.s1_work, 'R', periods, 0, 0);
    cd(here);
end
write_MINEOS_mod(card_final, fullfile(mdir, 'CARDS', [opt.s1_final '.card']));
copyfile(fullfile(mdir, 'CARDS', [opt.ref_qmod '.qmod']), fullfile(mdir, 'CARDS', [opt.s1_final '.qmod']));
fprintf('Rayleigh fundamental mode, reduced chi^2: start %.1f, final %.2f\n', ...
    mean(((cobs - cstart)./cstd).^2), mean(((cobs - cpre)./cstd).^2));
% linearized uncertainty in the format of the Bayesian versions (s1_unc):
% the data errors propagated through the damped generalized inverse
% (unc.Cm of the last iteration; the damping/smoothing are NOT included as
% a prior covariance, so this is not a posterior), mean +/- 1 and 1.645 std
zw_ = max(card_true.z(card_true.vsv == 0 & card_true.z < 20));
[zs_, io_] = sort(unc.z - zw_); km_ = zs_ >= 0;
vs_ = card_final.vsv(unc.I(io_))/1000; sd_ = unc.vsv_std(io_);
zplot1 = (0:2:400)'; ip_ = @(x) interp1(zs_(km_), x(km_), zplot1, 'linear', NaN);
s1_unc = struct('z', zplot1, 'vsv_mean', ip_(vs_), 'vsv_std', ip_(sd_), 'vsv_p05', ip_(vs_ - 1.645*sd_), ...
    'vsv_p16', ip_(vs_ - sd_), 'vsv_p50', ip_(vs_), 'vsv_p84', ip_(vs_ + sd_), 'vsv_p95', ip_(vs_ + 1.645*sd_), ...
    'data', 'linearized (Ginv Cd Ginv'')');
method = 'linearized';
save(opt.s1_file, 'card_start', 'card_final', 'card_true', 'cobs', 'cstd', ...
     'cstart', 'cpre', 'periods', 'discs', 'opt', 'method', 's1_unc');

%% figure
zw = max(card_true.z(card_true.vsv == 0 & card_true.z < 20));
zp = (0:1:400)';
cp = @(c, f) card_at_depth(c, f, zp + zw)/1000;
fig = figure('position', [50 50 1300 520]);
subplot(1, 3, 1); hold on; box on;
h1 = plot(cp(card_true, 'vsv'), zp, 'r-', 'linewidth', 2);
h2 = plot(cp(card_true, 'vsh'), zp, 'r--', 'linewidth', 1);
h4 = plot(cp(card_start, 'vsv'), zp, 'g-', 'linewidth', 1.5);
h5 = plot(cp(card_final, 'vsv'), zp, 'b-', 'linewidth', 2);
set(gca, 'ydir', 'reverse'); ylim([0 400]); xlim([3.9 5]); xlabel('V_S (km/s)'); ylabel('Depth below seafloor (km)');
legend([h1 h2 h4 h5], {'PA5 V_{SV}', 'PA5 V_{SH}', start_lbl, 'step 1 (isotropic)'}, ...
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
