function save_burn(fname, sv1, info1, Sf) %#ok<INUSD>
% SAVE_BURN  Save the burn-in phase of one chain/ensemble (callable from
% parfor): late burn-in samples sv1, info1 and the final states Sf of all
% temperatures (D3_step2_waveform_inv.m, opt.update_ref = 'burnin').
tmpf = tmp_name(fname);
save(tmpf, 'sv1', 'info1', 'Sf'); movefile(tmpf, fname);
end
