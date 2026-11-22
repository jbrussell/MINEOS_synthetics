function f = tmp_name(prefix)
% f = tmp_name(prefix): a temporary file name in the current folder that is
% unique between processes (process id + time), for write-then-rename saves.
% (Not from rand: sub_setup_records.m seeds the generator, so concurrent
% processes would draw the same number.)
if exist('OCTAVE_VERSION', 'builtin'), pid = getpid(); else, pid = feature('getpid'); end
f = sprintf('%s_tmp%d_%s.mat', prefix, pid, datestr(now, 'HHMMSSFFF'));
end
