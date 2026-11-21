function [Phi, Psi] = spec_basis(w, K, Ka)
% SPEC_BASIS  Frequency bases for event_G_spec.
%   Phi : Nf x K hat functions, uniform in log(w), summing to 1, for the
%         spectral correction B_n(w). K = 1 gives the scalar b_n.
%   Psi : Nf x Ka Legendre polynomials in xi = normalized log(w) on [-1, 1]
%         for the arrival angle: Ka = 1 constant dth_n, Ka = 2 adds a term
%         linear in log frequency, Ka = 3 a quadratic.
x = log(w(:)); xn = linspace(min(x), max(x), max(K, 2));
Phi = ones(numel(x), 1);
if K > 1
    Phi = zeros(numel(x), K);
    for k = 1:K, Phi(:, k) = max(0, 1 - abs(x - xn(k)) / (xn(2) - xn(1))); end
end
xi = 2*(x - min(x)) / (max(x) - min(x)) - 1;
P = [ones(size(xi)), xi, 1.5*xi.^2 - 0.5];
Psi = P(:, 1:Ka);
end
