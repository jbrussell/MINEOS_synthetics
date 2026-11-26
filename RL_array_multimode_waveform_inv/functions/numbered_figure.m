function fig = numbered_figure(n, varargin)
% fig = numbered_figure(n, 'property', value, ...)
%
% Figure number n, cleared, with the given properties (e.g. 'position').
% Every figure of the package has its own fixed number, so that rerunning a
% script redraws its figures in the same windows instead of opening new
% ones. Numbers (functions and drivers that use them):
%   1010-1095  D4 / sub_plot_results.m: 1010 fig1, 1020 fig2, 1030 fig3,
%              1040+e fig4 (event e), 1050 fig5, 1060 fig6, 1070 fig7,
%              1080 fig8, 1090 fig9 (plot_nuisance_truth.m), 1095 fig9_eta
%   1100       fig10 (plot_posterior_pdf.m)
%   1110, 1111 fig11, fig11b (plot_posterior_disp.m)
%   1200+e     fig12 (event e), 1250+e fig12b (plot_posterior_waveforms.m)
%   1301, 1302 fig13 Love, Rayleigh (plot_posterior_B.m)
%   1400       fig14 (plot_posterior_dth.m)
%   1501-1504  noise diagnostics (plot_noise_diagnostics.m)
%   1600       figs_compare.png (compare_runs.m)
%   2001-2004  D2 record sections (plot_record_sections.m): T all, T events,
%              Z all, Z events
%   2100+c     A4 event figure (component c), 2110+c A4 summary
%   3001       D1_step1_rayleigh_isotropic_inv, 3002 D1_step1_Bayesian_disp,
%   3003       D1_step1_Bayesian_waveform
%   3100       D5, 3200+k D7 (case k)
% The ASWMS-Q functions (functions/aswms/) keep their own numbers (1, 3, 36,
% 38, 39, 99, 101).
fig = figure(n);
clf(fig);
if ~isempty(varargin), set(fig, varargin{:}); end
end
