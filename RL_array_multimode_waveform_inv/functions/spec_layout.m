function [mu, sb2, it, ith] = spec_layout(lin_r, nm, K, Ka, sb, tau_prior, dth_prior)
% SPEC_LAYOUT  Prior mean/variance and parameter indices for event_G_spec.
%   lin_r     : struct with tau0 (nm x Kt, or 1 x nm for a constant tau_n)
%               and dth0 (nm x Ka) for one record
%   tau_prior : prior std of the time-shift coefficients (s), scalar or
%               1 x Kt, e.g. [15 10 10] = 15 s band mean, 10 s for the
%               linear and quadratic terms in normalized log(w)
%   dth_prior : 1 x Ka prior std of the arrival-angle coefficients (deg)
% The priors are centred on zero TOTAL tau and dth (hence the -lin terms).
% it  : nm x Kt indices of the dtau coefficients in the solution (index the
%       solution with reshape(xh(it), size(it)): with nm = 1, xh(it) of a column
%       xh would be a column, not 1 x Kt)
% ith : nm x Ka indices of the dth coefficients
nb  = 2*nm*K;
T0 = lin_r.tau0;
if size(T0, 1) == 1 && size(T0, 2) == nm && nm > 1, T0 = T0(:); end
Kt = size(T0, 2);
tp = tau_prior(:)'; if numel(tp) < Kt, tp = [tp, tp(end)*ones(1, Kt - numel(tp))]; end
t0 = T0.';                                            % Kt x nm -> mode-major
th0 = lin_r.dth0.';                                   % Ka x nm -> mode-major
mu  = [repmat([1; 0], nm*K, 1); -t0(:); -th0(:)];
sb2 = [sb^2*ones(nb, 1); repmat(tp(1:Kt)'.^2, nm, 1); repmat(dth_prior(1:Ka)'.^2, nm, 1)];
it  = reshape(nb + (1:nm*Kt), Kt, nm).';
ith = reshape(nb + nm*Kt + (1:nm*Ka), Ka, nm).';
end
