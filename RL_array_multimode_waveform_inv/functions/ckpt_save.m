function ckpt_save(f, S, sv, info, isv, it, rngst)
% ckpt_save(f, S, sv, info, isv, it, rngst)
% Checkpoint of run_pt_LR.m, written atomically (temporary file, then
% rename), so that a run killed while writing keeps the previous one.
tmpf = [f '.tmp.mat'];
save(tmpf, 'S', 'sv', 'info', 'isv', 'it', 'rngst', '-v7');
movefile(tmpf, f);
end
