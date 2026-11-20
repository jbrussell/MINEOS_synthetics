function v = bspline_profile(zc, vc, zq, nu)
% BSPLINE_PROFILE  Evaluate a 1-D profile defined by B-spline control points.
%
%   v = bspline_profile(zc, vc, zq)
%
%   zc, vc : control-point depths and values (any order; sorted by depth)
%   zq     : query depths (inside [min(zc), max(zc)])
%
% The profile is the parametric B-spline curve (z(u), v(u)), u in [0, 1],
% with control polygon (zc, vc), clamped uniform knots, and degree
% p = min(3, n-1) for n control points (cubic once n >= 4). Clamping makes
% the curve start and end at the first and last control points; because
% the control depths are sorted, z(u) is monotonic, so the curve is a
% single-valued, C2-smooth profile v(z). By the convex-hull property it
% never leaves the range of the control values. Adding or removing a
% control point changes the curve only locally (within p+1 neighbours),
% which keeps transdimensional proposals well behaved.
if nargin < 4, nu = 400; end
[zc, i] = sort(zc(:)); vc = vc(:); vc = vc(i);
n = numel(zc); p = min(3, n - 1);
u = linspace(0, 1, nu)';
N = bspline_basis(n, p, u);
zu = N*zc; vu = N*vc;
zu = zu + (0:nu-1)'*1e-10;           % strictly increasing for interp1
v = interp1(zu, vu, zq(:), 'linear', 'extrap');
v = reshape(v, size(zq));
end

function N = bspline_basis(n, p, u)
% Cox-de Boor basis (numel(u) x n) on a clamped uniform knot vector.
t = [zeros(1, p), linspace(0, 1, n - p + 1), ones(1, p)];
m = numel(u);
N = zeros(m, numel(t) - 1);
for j = 1:numel(t) - 1
    N(:, j) = (u >= t(j) & u < t(j+1));
end
N(u == 1, n) = 1; N(u == 1, n+1:end) = 0;   % right end of the last span
for d = 1:p
    Nn = zeros(m, numel(t) - 1 - d);
    for j = 1:numel(t) - 1 - d
        a = 0; b = 0;
        if t(j+d) > t(j),     a = (u - t(j)) / (t(j+d) - t(j)) .* N(:, j); end
        if t(j+d+1) > t(j+1), b = (t(j+d+1) - u) / (t(j+d+1) - t(j+1)) .* N(:, j+1); end
        Nn(:, j) = a + b;
    end
    N = Nn;
end
N = N(:, 1:n);
end
