function [c, K, gamma] = mineos_apply_Q(card, c0, K, periods, fref)
% [c, K, gamma] = mineos_apply_Q(card, c0, K, periods, fref)
%
% Apply the causal (physical dispersion) correction of a constant-Q model
% to elastic MINEOS phase velocities and phase-velocity kernels, using the
% first-order formulas of the CPS programs (subroutine gammap of
% tregn96/tlegn96 and cslmod of srfker96), as in calc_der96_TI.m and
% calc_anelastic_kernel96_TI.m of SURF96_Rayl_1Dinv:
%
%   x     = sum_r [ (dc/dVsv Vsv + dc/dVsh Vsh) Qs^-1
%                 + (dc/dVpv Vpv + dc/dVph Vph) Qp^-1 ] dr
%   c     = c0 + x ln(f/fref) / pi
%   gamma = w x / (2 c0^2)                       (attenuation coefficient)
%   dc/dVs -> dc/dVs (1 + ln(f/fref)/pi Qs^-1),  dc/dVp -> dc/dVp (1 + ln(f/fref)/pi Qp^-1)
%
% The velocities of the card are the values at the reference frequency
% fref. The eta and density kernels are left as their elastic values (as
% in srfker96). Q is taken from the card (Qs = Qmu, 1/Qp = L/Qmu +
% (1-L)/Qkappa with L = 4/3 (Vs/Vp)^2 and Voigt-average Vs and Vp; Qs^-1 = 0
% in fluids), as in card2mod_TI.m.
%
% INPUTS
% card    - MINEOS card structure (read_model_card.m; m/s, kg/m^3)
% c0      - elastic phase velocities (km/s) [N x 1]
% K       - elastic MINEOS kernels from run_mineos_phv_kernels.m (they
%           cover the top length(K.rad) nodes of the card) [Nnode x N]
% periods - periods (s) [N x 1]
% fref    - reference frequency (Hz) of the card velocities
%
% OUTPUTS
% c     - causal phase velocities (km/s) [N x 1]
% K     - causal kernels (same fields and units as the input)
% gamma - attenuation coefficients (1/km) [N x 1]
%
% jbrussell 11/2026

c0 = c0(:); periods = periods(:);
nnode = length(card.rad);
nK = length(K.rad);
iK = nnode - nK + (1:nK);

vsv = card.vsv(iK); vsh = card.vsh(iK); vpv = card.vpv(iK); vph = card.vph(iK);
qmu = card.qmu(iK); qkap = card.qkap(iK);
dr = gradient(card.rad);
dr = dr(iK);

% Q^-1 (Q <= 0 means elastic)
vs = sqrt((2*vsv.^2 + vsh.^2)/3);
vp = sqrt((vpv.^2 + 4*vph.^2)/5);
iqm = zeros(nK,1); iqm(qmu > 0) = 1./qmu(qmu > 0);
iqk = zeros(nK,1); iqk(qkap > 0) = 1./qkap(qkap > 0);
L = 4/3 * (vs./vp).^2;
qa = L.*iqm + (1-L).*iqk; % Qp^-1
qb = iqm;                 % Qs^-1
qb(vsv <= 10) = 0;        % fluid (Vsv <= 0.01 km/s)

N = length(periods);
c = nan(N,1);
gamma = nan(N,1);
for ip = 1:N
    if isnan(c0(ip)); continue; end
    f = 1/periods(ip); w = 2*pi*f;
    faclog = log(f/fref)/pi;
    sumS = (K.vsv(:,ip).*vsv + K.vsh(:,ip).*vsh) .* dr;
    sumP = (K.vpv(:,ip).*vpv + K.vph(:,ip).*vph) .* dr;
    sumS(isnan(sumS)) = 0; sumP(isnan(sumP)) = 0;
    x = sum(sumS.*qb + sumP.*qa); % km/s
    c(ip) = c0(ip) + x*faclog;
    gamma(ip) = 0.5*w*x/(c0(ip)^2);
    K.vsv(:,ip) = K.vsv(:,ip).*(1 + faclog*qb);
    K.vsh(:,ip) = K.vsh(:,ip).*(1 + faclog*qb);
    K.vpv(:,ip) = K.vpv(:,ip).*(1 + faclog*qa);
    K.vph(:,ip) = K.vph(:,ip).*(1 + faclog*qa);
end
end
