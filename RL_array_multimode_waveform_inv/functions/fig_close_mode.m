function c = fig_close_mode()
% c = fig_close_mode()
%
% Whether the scripts close their figures after saving them (and run
% "close all" at the start): opt.close_figs of parameters_RL.m (stored as
% application data of the root, 'RL_close_figs'); empty (default) = auto:
% true when running headless - Octave (the scripts make invisible figures
% there) or MATLAB without its desktop (-batch, -nodisplay) - so that the
% open figures do not pile up in long runs; false in interactive MATLAB,
% where the figures stay open on screen.
c = getappdata(0, 'RL_close_figs');
if isempty(c), c = exist('OCTAVE_VERSION', 'builtin') ~= 0 || ~usejava('desktop'); end
c = logical(c);
end
