function [lin, cR, card] = s1_lin_at(mdir, card, cardid, qmodid, periods, maxF, zw, zmin, zmax)
% [lin, cR, card] = s1_lin_at(mdir, card, cardid, qmodid, periods, maxF, zw, zmin, zmax)
%
% Linearized MINEOS Rayleigh dispersion (fundamental mode) about a card,
% for the Vsv-only step-1 inversions (D1_step1_Bayesian_disp.m,
% D1_step1_Bayesian_waveform.m): writes the card to
% <mdir>/CARDS/<cardid>.card (with a copy of <qmodid>.qmod), runs MINEOS
% for the spheroidal branch-0 phase velocities and kernels at the periods
% (mineos_modes.m) and builds the linearization between zmin and zmax km
% below the seafloor (mineos_lin_setup.m; zw = water depth, km).
%
% lin  : linearization struct (mineos_lin_setup.m). Only the Rayleigh part
%        is meaningful: the Love fields are copies of the Rayleigh ones
%        (Love data are never used in step 1).
% cR   : [Nper x 1] branch-0 Rayleigh phase velocities of the card (km/s,
%        elastic; NaN where the branch does not exist)
% card : the card as read back from the file (read_model_card.m)
write_MINEOS_mod(card, fullfile(mdir, 'CARDS', [cardid '.card']));
qsrc = fullfile(mdir, 'CARDS', [qmodid '.qmod']); qdst = fullfile(mdir, 'CARDS', [cardid '.qmod']);
if ~strcmp(qsrc, qdst), copyfile(qsrc, qdst); end
card = read_model_card(fullfile(mdir, 'CARDS', [cardid '.card']));
mineos_lock(mdir, 'acquire');                    % (shared MINEOS scratch files: mineos_lock.m)
try
    o = mineos_modes(mdir, cardid, 'S', 1, maxF, periods(:), 0);
catch err
    mineos_lock(mdir, 'release'); rethrow(err);
end
mineos_lock(mdir, 'release');
cR = o.phv(:, 1);
lin = mineos_lin_setup(card, o.K, o.K, cR, cR, zw, zmin, zmax);
end
