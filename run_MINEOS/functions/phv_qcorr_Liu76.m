function [phvq] = phv_qcorr_Liu76(phv,q,periods,fref_hz)
    %
    % Correct phase velocities for physical dispersion following:
    % Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
    % anelasticity: implications for seismology and mantle composition,
    % Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
    %
    % jbrussell 5/3/2022

    w = 2*pi./periods(:);
    w_ref = 2*pi*fref_hz;
    phvq = phv(:) .* ( 1 + 1./(pi.*q(:)) .* log(w(:)./w_ref) );

end
