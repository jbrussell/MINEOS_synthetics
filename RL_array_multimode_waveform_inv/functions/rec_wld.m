function w = rec_wld(rec)
% w = rec_wld(rec)  log det of the record's whitening operator (rec.wld), or 0
w = 0;
if isfield(rec, 'wld') && ~isempty(rec.wld), w = rec.wld; end
end
