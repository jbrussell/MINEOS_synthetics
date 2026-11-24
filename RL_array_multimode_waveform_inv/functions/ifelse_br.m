function b = ifelse_br(br, nm)
% b = ifelse_br(br, nm)  branches used by D7: br, or 0..nm-1 if br is empty
if isempty(br), b = 0:nm-1; else, b = br; end
end
