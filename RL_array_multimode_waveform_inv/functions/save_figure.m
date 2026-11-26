function f = save_figure(fig, file, res)
% f = save_figure(fig, file, res)
%
% Save figure fig as a PNG (resolution res dpi, default 130) under file,
% with the script that made it in front of the file name: <dir>/<S>_<name>,
% where S is the driver (A0..A4, D1..D99: the first part of its name, e.g.
% D4 for D4_step2_plot_results.m) that is running - the innermost one when
% one driver calls another (D3 calls D4: the D4 figures get D4_). Without
% a driver on the call stack (a function called by hand) the name is kept.
% Returns the file written.
if nargin < 3 || isempty(res), res = 130; end
[d, n, e] = fileparts(file);
if isempty(e), e = '.png'; end
p = driver_prefix_();
if ~isempty(p) && ~strncmp(n, [p '_'], numel(p) + 1), n = [p '_' n]; end
f = fullfile(d, [n e]);
print(fig, '-dpng', sprintf('-r%d', res), f);
end

function p = driver_prefix_()
p = '';
s = dbstack;
for k = 1:numel(s)
    t = regexp(s(k).name, '^([AD]\d+)_', 'tokens', 'once');
    if ~isempty(t), p = t{1}; return; end
end
end
