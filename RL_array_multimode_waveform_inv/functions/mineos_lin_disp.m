function [CL, CR, X] = mineos_lin_disp(L, vsv, vsh, eta, which)
% [CL, CR, X] = mineos_lin_disp(L, vsv, vsh, eta, which)
%
% Linearized MINEOS phase velocities [Nper x nm] (Love CL, Rayleigh CR)
% for Vsv and Vsh (km/s) and eta at the model nodes of L
% (mineos_lin_setup.m):
%
%   c_n = c_ref,n + K_vsv dVsv + K_vsh dVsh + K_vpv dVpv + K_vph dVph
%                 + K_rho drho + K_eta deta,
%
% with every parameter's own kernel and perturbation dX = X - X_ref:
%   Vpv = (Vpv/Vsv)_ref Vsv,  Vph = (Vph/Vsh)_ref Vsh,  rho = (rho/Vsv)_ref Vsv,
%   eta = the reference eta if eta is empty, else the given profile.
% which = 'both' (default), 'love' or 'rayleigh' ([] for the other).
% X: struct of the node values (vsv, vsh, vpv, vph, rho, eta).
if nargin < 5 || isempty(which), which = 'both'; end
if isempty(eta), eta = L.eta0; end
X.vsv = vsv(:); X.vsh = vsh(:);
X.vpv = L.rpv .* X.vsv; X.vph = L.rph .* X.vsh; X.rho = L.rrho .* X.vsv;
X.eta = eta(:);
d.vsv = X.vsv - L.vsv0; d.vsh = X.vsh - L.vsh0;
d.vpv = X.vpv - L.vpv0; d.vph = X.vph - L.vph0;
d.rho = X.rho - L.rho0; d.eta = X.eta - L.eta0;
f = {'vsv', 'vsh', 'vpv', 'vph', 'rho', 'eta'};
nm = numel(L.KL);
CL = []; CR = [];
if ~strcmp(which, 'rayleigh')
    CL = L.cL;
    for n = 1:nm
        for i = 1:numel(f), CL(:, n) = CL(:, n) + L.KL{n}.(f{i}) * d.(f{i}); end
    end
end
if ~strcmp(which, 'love')
    CR = L.cR;
    for n = 1:nm
        for i = 1:numel(f), CR(:, n) = CR(:, n) + L.KR{n}.(f{i}) * d.(f{i}); end
    end
end
end
