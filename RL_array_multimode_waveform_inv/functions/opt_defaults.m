function opt = opt_defaults(opt, def)
% opt = opt_defaults(opt, def)
% Fill the fields of struct opt that are missing with those of def.
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt, f{i}), opt.(f{i}) = def.(f{i}); end
end
end
