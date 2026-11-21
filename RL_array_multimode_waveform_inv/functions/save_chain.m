function save_chain(fname, sv, st, ic) %#ok<INUSD>
% SAVE_CHAIN  Save one chain's output (callable from inside parfor).
save(fname, 'sv', 'st', 'ic');
end
