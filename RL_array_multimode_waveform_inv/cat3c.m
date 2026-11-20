function X = cat3c(c)
% CAT3C  Concatenate the arrays in cell array C along dimension 3.
X = cat(3, c{:});
end
