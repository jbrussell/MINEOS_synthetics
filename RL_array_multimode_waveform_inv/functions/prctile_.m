function q = prctile_(X, p)
% PRCTILE_  Percentile of each row of X across columns (no toolbox needed).
X = sort(X, 2); n = size(X, 2);
pos = 1 + (n - 1)*p/100;
lo = floor(pos); hi = min(lo + 1, n); fr = pos - lo;
q = (1 - fr)*X(:, lo) + fr*X(:, hi);
end
