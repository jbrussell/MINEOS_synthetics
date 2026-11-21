function mineos_lock(mdir, action)
% mineos_lock(mdir, 'acquire' | 'release')
% Serialize MINEOS/idagrn6 runs between processes: they share scratch
% files with fixed names in mdir (run_MINEOS). The lock is the directory
% <mdir>/.prep_lock (mkdir is atomic). If a process dies while holding it,
% remove that directory by hand.
d = fullfile(mdir, '.prep_lock');
switch action
    case 'acquire'
        while true
            if exist(d, 'dir') ~= 7 && system(sprintf('mkdir "%s" 2>/dev/null', d)) == 0, return; end
            fprintf('waiting for another MINEOS/idagrn6 run (%s)\n', d); pause(30);
        end
    case 'release'
        if exist(d, 'dir') == 7, rmdir(d); end
end
end
