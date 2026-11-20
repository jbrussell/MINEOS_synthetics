function [mu, sb2, it, ith] = spec_layout(lin_r, nm, K, Ka, sb, tau_prior, dth_prior)
% SPEC_LAYOUT  Prior mean/variance and parameter indices for event_G_spec.
%   lin_r     : struct with tau0 (1 x nm) and dth0 (nm x Ka) for one record
%   dth_prior : 1 x Ka prior std of the arrival-angle coefficients (deg),
%               e.g. [5 3] = 5 deg constant, 3 deg linear-in-log(w) term
% The priors are centred on zero TOTAL tau and dth (hence the -lin terms),
% as in run_pt_LR. it, ith index dtau and dth (nm x Ka) in the solution.
nb  = 2*nm*K;
th0 = lin_r.dth0.';                                   % Ka x nm -> mode-major
mu  = [repmat([1; 0], nm*K, 1); -lin_r.tau0(:); -th0(:)];
sb2 = [sb^2*ones(nb, 1); tau_prior^2*ones(nm, 1); repmat(dth_prior(1:Ka)'.^2, nm, 1)];
it  = nb + (1:nm);
ith = reshape(nb + nm + (1:nm*Ka), Ka, nm).';
end
