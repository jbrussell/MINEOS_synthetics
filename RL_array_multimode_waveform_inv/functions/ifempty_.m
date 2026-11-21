function y = ifempty_(a, b)
% y = ifempty_(a, b): a, or b if a is empty
if isempty(a), y = b; else, y = a; end
end
