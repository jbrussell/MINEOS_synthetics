function W = rec_W(rec)
% W = rec_W(rec)  whitening operator of a record (rec.W, AR(p) noise models;
% sub_setup_records.m), or [] (AR(1) with the sampled correlation)
W = [];
if isfield(rec, 'W') && ~isempty(rec.W), W = rec.W; end
end
