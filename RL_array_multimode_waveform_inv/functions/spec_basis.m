function [Phi, Psi, Psit] = spec_basis(w, K, Ka, Kt)
% SPEC_BASIS  Frequency bases for event_G_spec.
%   Phi  : Nf x K hat functions, uniform in log(w), summing to 1, for the
%          spectral correction B_n(w). K = 1 gives the scalar b_n.
%   Psi  : Nf x Ka Legendre polynomials in xi = normalized log(w) on [-1, 1]
%          for the arrival angle dth_n(w): Ka = 1 constant, Ka = 2 adds a
%          term linear in log frequency, Ka = 3 a quadratic.
%   Psit : Nf x Kt Legendre polynomials (as Psi) for the time shift
%          tau_n(w) (OPTIONAL, default Kt = 1: a constant time shift).
if nargin < 4, Kt = 1; end
x = log(w(:)); xn = linspace(min(x), max(x), max(K, 2));
Phi = ones(numel(x), 1);
if K > 1
    Phi = zeros(numel(x), K);
    for k = 1:K, Phi(:, k) = max(0, 1 - abs(x - xn(k)) / (xn(2) - xn(1))); end
end
xi = 2*(x - min(x)) / (max(x) - min(x)) - 1;
P = [ones(size(xi)), xi, 1.5*xi.^2 - 0.5, 2.5*xi.^3 - 1.5*xi];
Psi = P(:, 1:Ka);
Psit = P(:, 1:Kt);
end
