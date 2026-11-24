function nl = rec_nlev(rec)
% nl = rec_nlev(rec)  per-sample noise levels of a record (rec.nlev, from the
% measured pre-event noise; sub_setup_records.m), or [] if it has none
nl = [];
if isfield(rec, 'nlev') && ~isempty(rec.nlev), nl = rec.nlev; end
end
