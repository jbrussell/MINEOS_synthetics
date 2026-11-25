function plot_font_defaults(fs)
% plot_font_defaults(fs)
%
% Font size (pt) of all plot axes: sets the default font size of new axes
% (tick labels, axis labels, titles) and text objects to fs (default
% 15; opt.plot_fontsize). Called by parameters_RL.m, so it applies to
% every figure of the package scripts; call it before using a plotting
% function on its own (e.g. plot_noise_diagnostics.m).
if nargin < 1 || isempty(fs), fs = 15; end
set(0, 'defaultaxesfontsize', fs);
set(0, 'defaulttextfontsize', fs);
end
