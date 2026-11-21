function [L, bhat] = event_loglike(d, G, wtype, sig, mu_b, sb2)
% EVENT_LOGLIKE  Log marginal likelihood for one event with the modal
% amplitudes integrated out analytically.
%
%   d = G*b + n,  b ~ N(mu_b, diag(sb2)),  n ~ N(0, diag(sig(wtype).^2))
%   => d ~ N(G*mu_b, Cd + G*Cb*G')
%
% Evaluated with the Woodbury identity so only a (2*nmodes)^2 system is
% solved. bhat is the posterior mean of b (useful for plotting fits).

Wd  = 1 ./ reshape(sig(wtype), [], 1).^2;
GW  = G .* Wd;                         % = Cd^-1 * G
A   = diag(1./sb2) + G' * GW;          % = Cb^-1 + G' Cd^-1 G
r0  = d - G*mu_b;
b   = GW' * r0;
Rc  = chol(A);
y   = Rc' \ b;
quad   = sum(r0.^2 .* Wd) - y'*y;
logdet = -sum(log(Wd)) + sum(log(sb2)) + 2*sum(log(diag(Rc)));
L = -0.5*(quad + logdet + numel(d)*log(2*pi));
if nargout > 1
    bhat = Rc \ (Rc' \ (GW'*d + mu_b./sb2));
end
end
