function [winpara, W] = auto_win_select(u, t, dist, periods, o)
% [winpara, W] = auto_win_select(u, t, dist, periods, o)
%
% Automatic surface-wave window of one event recorded by an array: the
% data pick the window. Port of auto_win_select.m of ASWMS-Q
% (github.com/jbrussell/ASWMS-Q, a3_a_run_autowinpick.m; after Jin & Gaherty
% 2015), with the same algorithm and defaults:
%   1. per station, narrow-band Gaussian filters at o.bandnum frequencies
%      spanning 1/periods; envelopes |analytic signal|, zero outside the
%      group-velocity range dist/o.max_groupv .. dist/o.min_groupv, each
%      normalized by its maximum;
%   2. local envelope maxima above o.peakamptol of the largest are
%      candidate arrivals; starting from the largest peak at the centre
%      frequency, the best peak is tracked to lower and higher frequencies,
%      allowing a group-velocity change of o.groupv_diff_tol between
%      neighbouring bands (higher frequencies: a score of amplitude +
%      group-velocity continuity, favouring normal dispersion);
%   3. per frequency, a weighted (peak amplitude / sum of the peak
%      amplitudes), outlier-trimmed straight-line fit of the group delays
%      against distance gives the group velocity and an offset;
%   4. per station the window runs from the earliest (over frequency) of
%      peak time - o.cycle_before periods to the latest of peak time +
%      o.cycle_after periods; straight lines fitted to these against
%      distance give the window dist/v1 + t1 .. dist/v2 + t2.
%
% u       - [Nt x Nsta] displacement (any processing; the narrow-band
%           filters define the band)
% t       - [Nt x 1] time after origin (s), uniform sampling
% dist    - [Nsta x 1] epicentral distances (km)
% periods - periods spanning the band (s); the frequencies are
%           linspace(1/max(periods), 1/min(periods), o.bandnum)
% o       - optional struct, fields (defaults as in ASWMS-Q):
%           min_groupv (2), max_groupv (5) km/s, bandnum (20), cent_freq
%           (mean(1./periods)), min_width/max_width (0.1/0.1, Gaussian
%           width relative to the centre frequency), cycle_before (2),
%           cycle_after (5), peakamptol (0.5), peak_search_range (1),
%           groupv_diff_tol (0.2), positive_disp_weight (5), tmax_cap (Inf)
%
% winpara - [v1 t1 v2 t2]: window dist/v1 + t1 .. dist/v2 + t2 (s)
% W       - diagnostics: freqs, groupdelay [nf x Nsta], snr, groupv,
%           offset [nf], bgtime, endtime [Nsta], env {Nsta} [Nt x nf]
%           normalized envelopes
%
% jbrussell 11/2026 (port)
if nargin < 5, o = struct(); end
d_ = struct('min_groupv', 2, 'max_groupv', 5, 'bandnum', 20, 'cent_freq', mean(1./periods), ...
    'min_width', 0.1, 'max_width', 0.1, 'cycle_before', 2, 'cycle_after', 5, 'peakamptol', 0.5, ...
    'peak_search_range', 1, 'groupv_diff_tol', 0.2, 'positive_disp_weight', 5, 'tmax_cap', Inf);
fn = fieldnames(d_);
for k = 1:numel(fn), if ~isfield(o, fn{k}) || isempty(o.(fn{k})), o.(fn{k}) = d_.(fn{k}); end, end
t = t(:); dist = dist(:); [Nt, Nsta] = size(u); dt = t(2) - t(1);
freqs = linspace(1/max(periods), 1/min(periods), o.bandnum); nf = numel(freqs);
[~, ic] = min(abs(freqs - o.cent_freq));
% Gaussian filters on the positive frequencies (analytic signal)
faxis = (0:floor(Nt/2))'/dt/Nt;
if nf > 1
    width = o.max_width - (o.max_width - o.min_width)/(max(freqs) - min(freqs))*(freqs - min(freqs));
else
    width = o.max_width;
end
G = zeros(Nt, nf);
for ip = 1:nf
    G(1:numel(faxis), ip) = exp(-(faxis - freqs(ip)).^2/2./(width(ip)*freqs(ip)).^2);
end
groupdelay = zeros(nf, Nsta); snr = zeros(nf, Nsta); env = cell(Nsta, 1);
tol = o.groupv_diff_tol; pw = o.positive_disp_weight; sr = o.peak_search_range;
for is = 1:Nsta
    U = fft(u(:, is));
    E = abs(ifft(U .* G));                         % [Nt x nf] envelopes
    tmin = dist(is)/o.max_groupv; tmax = min(dist(is)/o.min_groupv, o.tmax_cap);
    inw = find(t > tmin & t < tmax);
    inw = inw(inw > 1 & inw < Nt);
    E(t <= tmin | t >= tmax, :) = 0;
    E = E ./ max(max(E, [], 1), realmin);
    env{is} = E;
    pk = struct('t', {}, 'a', {}, 'best', {}, 'snr', {});
    for ip = 1:nf
        e = E(:, ip); de = diff(e);
        k = inw(de(inw - 1) > 0 & de(inw) < 0);
        k = k(e(k) > max(e)*o.peakamptol);
        pk(ip).t = t(k); pk(ip).a = e(k); pk(ip).best = []; pk(ip).snr = 0;
    end
    if ~isempty(pk(ic).t)
        [~, k] = max(pk(ic).a);
        pk(ic).best = pk(ic).t(k); pk(ic).snr = pk(ic).a(k)/sum(pk(ic).a);
    end
    % to lower frequencies: largest peak within the group-velocity tolerance
    for ip = ic-1:-1:1
        if isempty(pk(ip).t) || isempty(pk(ip+1).best)
            pk(ip).best = pk(ip+1).best; pk(ip).snr = 0; continue;
        end
        [~, kc] = min(abs(pk(ip+1).best - pk(ip).t)); kb = kc; ab = pk(ip).a(kb);
        V0 = dist(is)/pk(ip+1).best;
        for k = kc-sr:kc+sr
            if k >= 1 && k <= numel(pk(ip).t)
                v = dist(is)/pk(ip).t(k);
                if v > V0*(1 - tol) && v < V0*(1 + tol) && pk(ip).a(k) > ab, kb = k; ab = pk(ip).a(k); end
            end
        end
        v = dist(is)/pk(ip).t(kb);
        if v < V0*(1 - tol) || v > V0*(1 + tol)
            pk(ip).best = pk(ip+1).best; pk(ip).snr = 0;
        else
            pk(ip).best = pk(ip).t(kb); pk(ip).snr = pk(ip).a(kb)/sum(pk(ip).a);
        end
    end
    % to higher frequencies: amplitude + continuity score (normal dispersion favoured)
    for ip = ic+1:nf
        if isempty(pk(ip).t) || isempty(pk(ip-1).best)
            pk(ip).best = pk(ip-1).best; pk(ip).snr = 0; continue;
        end
        [~, kc] = min(abs(pk(ip-1).best - pk(ip).t)); kb = kc;
        V0 = dist(is)/pk(ip-1).best; v = dist(is)/pk(ip).t(kb);
        if v/V0 > 1, vp = 1 - (v/V0 - 1)/tol; else, vp = 1 - (1 - v/V0)/tol/pw; end
        best = pk(ip).a(kb) + vp;
        for k = kc-sr:kc+sr
            if k >= 1 && k <= numel(pk(ip).t)
                v = dist(is)/pk(ip).t(k);
                if v > V0*(1 - tol) && v < V0*(1 + tol)
                    if v/V0 > 1, vp = 1 - (v/V0 - 1)*10; else, vp = 1 - (1 - v/V0)*10/pw; end
                    if pk(ip).a(k) + vp > best, kb = k; best = pk(ip).a(k) + vp; end
                end
            end
        end
        v = dist(is)/pk(ip).t(kb);
        if v < V0*(1 - tol) || v > V0*(1 + tol)
            pk(ip).best = pk(ip-1).best; pk(ip).snr = 0;
        else
            pk(ip).best = pk(ip).t(kb); pk(ip).snr = pk(ip).a(kb)/sum(pk(ip).a);
        end
    end
    for ip = 1:nf
        if isempty(pk(ip).best), pk(ip).best = 0; pk(ip).snr = 0; end
    end
    groupdelay(:, is) = [pk.best]'; snr(:, is) = [pk.snr]';
end
groupv = zeros(1, nf); offset = zeros(1, nf);
for ip = 1:nf
    [groupv(ip), offset(ip)] = groupv_fit(dist, groupdelay(ip, :)', snr(ip, :)', o.min_groupv, o.max_groupv);
end
bgtime = zeros(Nsta, 1); endtime = bgtime;
for is = 1:Nsta
    pt = dist(is)./groupv + offset;
    bgtime(is) = min(pt - o.cycle_before./freqs);
    endtime(is) = max(pt + o.cycle_after./freqs);
end
p1 = polyfit(dist, bgtime, 1); p2 = polyfit(dist, endtime, 1);
winpara = [1/p1(1), p1(2), 1/p2(1), p2(2)];
W = struct('freqs', freqs, 'groupdelay', groupdelay, 'snr', snr, 'groupv', groupv, ...
    'offset', offset, 'bgtime', bgtime, 'endtime', endtime);
W.env = env;
end

function [groupv, offset] = groupv_fit(x, y, w, vmin, vmax)
% weighted straight-line fit y = x/groupv + offset, two passes of outlier
% removal (2 standard deviations), groupv clipped to [vmin vmax]
A = [x(:), ones(numel(x), 1)]; w = w(:); y = y(:);
if sum(w > 0) < 2, groupv = vmin; offset = 0; return; end
p = lsq_(A, y, w);
for pass = 1:2
    e = A*p - y; e(w == 0) = 0;
    w(abs(e) > 2*std(e)) = 0;
    if sum(w > 0) < 2, break; end
    p = lsq_(A, y, w);
end
groupv = 1/p(1); offset = p(2);
if ~(groupv <= vmax), groupv = vmax; offset = sum((y - x(:)/groupv).*w)/sum(w); end
if groupv < vmin, groupv = vmin; offset = sum((y - x(:)/groupv).*w)/sum(w); end
end

function p = lsq_(A, y, w)
W = diag(w); p = (A'*W*A) \ (A'*W*y);
end
