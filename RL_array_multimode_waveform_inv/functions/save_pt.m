function save_pt(fname, sv, st, ic, info) %#ok<INUSD>
% SAVE_PT  Save one parallel-tempering ensemble (callable from parfor).
save(fname, 'sv', 'st', 'ic', 'info');
end
