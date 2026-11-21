function P = update_reference(tag, prepfile_new, refid_new)
% P = update_reference(tag, prepfile_new, refid_new)
%
% Iterated reference model: build a new reference card from the posterior
% mean of a finished run (setup<tag>.mat and chain<tag>_XX.mat) and rerun
% MINEOS and idagrn6 for it (prep_mineos_idagrn.m with opt.ref_card; the
% true-model data are reused). The next stage then uses
%   - reference spectra Uref (idagrn6_mask branches) and the path-
%     correction reference of the updated model, and
%   - MINEOS phase velocities and kernels linearized about it,
% which removes most of the linearization error and the branch mismatch
% (e.g. the near-crossing of the Love fundamental and first overtone in
% PA5-like lithosphere) of a reference far from the truth.
%
% The new card is the current reference card with Vsv and Vsh at the nodes
% between cfg.zmin and cfg.zmax (below the seafloor) replaced by the
% posterior mean; Vp, density, eta and everything outside that depth range
% are unchanged.
%
%   P = update_reference('_s1', 'prep_s2.mat', 'RLinv_ref_s2')
%
% jbrussell 11/2026
S = load(['setup' tag '.mat'], 'cfg', 'mc', 'zplot', 'opt');
cfg = S.cfg; mc = S.mc; zplot = S.zplot; opt = S.opt;
V = 0; Vh = 0; ns = 0;
for ic = 1:mc.nchains
    D = load(sprintf('chain%s_%02d.mat', tag, ic));
    V = V + sum(D.sv.vsv, 2); Vh = Vh + sum(D.sv.vsh, 2); ns = ns + size(D.sv.vsv, 2);
end
vsv = V/ns; vsh = Vh/ns;
Pold = prep_mineos_idagrn(opt);
c = Pold.card_ref;
zsf = c.z - cfg.zw;
in = find(zsf >= cfg.zmin - 1e-6 & zsf <= cfg.zmax + 1e-6 & c.vsv > 0);
% mantle side at the Moho: the crust-side node keeps its value
atm = in(abs(zsf(in) - cfg.zmin) < 1e-6);
if numel(atm) > 1, in = setdiff(in, atm(2:end)); end
im = zplot >= cfg.zmin;                          % mantle samples only
c.vsv(in) = 1000*interp1(zplot(im), vsv(im), zsf(in), 'linear', 'extrap');
c.vsh(in) = 1000*interp1(zplot(im), vsh(im), zsf(in), 'linear', 'extrap');
o = opt; o.prepfile = prepfile_new; o.refid = refid_new; o.ref_card = c; o.P_true = Pold;
P = prep_mineos_idagrn(o);
end
