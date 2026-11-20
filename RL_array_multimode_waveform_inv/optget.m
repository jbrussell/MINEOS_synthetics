function v = optget(opt, f, d)
% OPTGET  Field F of struct OPT if present, otherwise default D.
if isfield(opt, f), v = opt.(f); else, v = d; end
end
