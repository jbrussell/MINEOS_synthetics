function rm_if_exists(f)
% rm_if_exists(f)  delete file f if it exists (e.g. a finished checkpoint)
if exist(f, 'file') == 2, delete(f); end
end
