function mineos_lock(mdir, action)
% mineos_lock(mdir, 'acquire' | 'release')
% Serialize MINEOS/idagrn6 runs between processes: they share scratch
% files with fixed names in mdir (run_MINEOS). The lock is the directory
% <mdir>/.prep_lock (mkdir is atomic); its file owner.txt records the
% holder (process ID, host, start time).
%
% Stale locks are removed automatically before waiting:
%   - the holder ran on this host and its process no longer exists
%     (a crashed or killed MATLAB/Octave);
%   - the holder is this process itself: the lock is not re-entrant, so it
%     is left from an interrupted call (Ctrl-C in interactive MATLAB, after
%     which the process lives on);
%   - the lock has no owner.txt (made by an older version of this
%     function, or a holder killed between mkdir and writing it), is older
%     than 2 min, and no MINEOS or idagrn6 program is running.
% A lock held by a live process, or by a process on another host (shared
% file system), is waited for; the message names the holder. If a live
% holder is stuck, stop it, or remove <mdir>/.prep_lock by hand.
d = fullfile(mdir, '.prep_lock');
of = fullfile(d, 'owner.txt');
me = getpid_();
switch action
    case 'acquire'
        t_msg = -Inf;
        while true
            if exist(d, 'dir') ~= 7 && system(sprintf('mkdir "%s" 2>/dev/null', d)) == 0
                fid = fopen(of, 'w');
                if fid > 0
                    fprintf(fid, '%d\n%s\n%s\n', me, host_(), datestr(now, 'yyyy-mm-dd HH:MM:SS')); fclose(fid);
                end
                return
            end
            [stale, who] = is_stale_(d, of, me);
            if stale
                fprintf('removing stale MINEOS/idagrn6 lock (%s): %s\n', d, who);
                remove_(d, of);
                continue
            end
            if toc_(t_msg) >= 600                         % message every 10 min
                fprintf('waiting for another MINEOS/idagrn6 run (%s; %s)\n', d, who);
                t_msg = now;
            end
            pause(30);
        end
    case 'release'
        if exist(d, 'dir') ~= 7, return; end
        [pid, h] = owner_(of);
        if isempty(pid) || (pid == me && strcmp(h, host_()))
            remove_(d, of);
        else
            warning('mineos_lock: %s is held by process %d on %s, not by this process (%d): left in place', ...
                    d, pid, h, me);
        end
end
end

function [stale, who] = is_stale_(d, of, me)
stale = false;
[pid, h, t0] = owner_(of);
if isempty(pid)
    age = lock_age_min_(d);
    who = sprintf('no owner record, made %.0f min ago', age);
    if age > 2 && ~mineos_running_()
        stale = true; who = [who ', no MINEOS/idagrn6 program running'];
    end
    return
end
who = sprintf('held by process %d on %s since %s', pid, h, t0);
if ~strcmp(h, host_()), return; end                 % another host: cannot check
if pid == me
    stale = true; who = [who ' (this process: left by an interrupted run)'];
elseif ~alive_(pid)
    stale = true; who = [who ' (that process no longer exists)'];
end
end

function [pid, h, t0] = owner_(of)
pid = []; h = ''; t0 = '';
fid = fopen(of, 'r');
if fid < 0, return; end
c = textscan(fid, '%s', 'delimiter', '\n'); fclose(fid); c = c{1};
if numel(c) >= 1, pid = str2double(c{1}); if ~isfinite(pid), pid = []; end, end
if numel(c) >= 2, h = strtrim(c{2}); end
if numel(c) >= 3, t0 = strtrim(c{3}); end
end

function remove_(d, of)
if exist(of, 'file') == 2, delete(of); end
if exist(d, 'dir') == 7
    try
        rmdir(d);
    catch
        warning('mineos_lock: could not remove %s (not empty?); remove it by hand', d);
    end
end
end

function a = lock_age_min_(d)
a = 0;
e = dir(d);
i = find(strcmp({e.name}, '.'), 1);
if ~isempty(i), a = (now - e(i).datenum)*24*60; end
end

function tf = alive_(pid)
if ispc
    [~, out] = system(sprintf('tasklist /FI "PID eq %d" /NH', pid));
    tf = ~isempty(strfind(out, sprintf(' %d ', pid)));
else
    tf = system(sprintf('ps -p %d >/dev/null 2>&1', pid)) == 0;
end
end

function tf = mineos_running_()
% any MINEOS program (mineos_*, eig_recover, frechet*) or idagrn6 running on
% this host (by process name)
tf = false;
if ispc, return; end
[st, ~] = system('pgrep "^(mineos_.*|idagrn.*|eig_recover|frechet.*)$" >/dev/null 2>&1');
tf = st == 0;
end

function p = getpid_()
if exist('OCTAVE_VERSION', 'builtin'), p = getpid();
else, p = feature('getpid'); end
end

function h = host_()
persistent h_
if isempty(h_)
    [st, out] = system('hostname');
    if st == 0, h_ = strtrim(out); else, h_ = 'unknown'; end
end
h = h_;
end

function s = toc_(t)
s = (now - t)*86400;
end
