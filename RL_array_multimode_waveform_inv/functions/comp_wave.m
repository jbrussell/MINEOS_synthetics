function wv = comp_wave(opt, comp)
% wv = comp_wave(opt, comp)
% Wave type of a component: 1 (Love) or 2 (Rayleigh), from its place in
% opt.aswms_components ({Love Rayleigh}, default {'BHT', 'BHZ'}); a name
% not listed there: Love if it ends in T, else Rayleigh.
c = optget(opt, 'aswms_components', {'BHT', 'BHZ'});
wv = find(strcmp(c, comp), 1);
if isempty(wv), if ~isempty(comp) && upper(comp(end)) == 'T', wv = 1; else, wv = 2; end, end
end
