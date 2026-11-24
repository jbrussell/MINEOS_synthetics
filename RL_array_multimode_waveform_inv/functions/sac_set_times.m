function sac_set_times(fn, B, delta, npts)
% sac_set_times(fn, B, delta, npts)
% Set the begin time B (and E = B + (npts-1)*delta, O = 0) of a
% little-endian SAC file written by mksac.m (which forces B = 0): header
% words 5 (B), 6 (E) and 7 (O), 4-byte floats.
fid = fopen(fn, 'r+', 'ieee-le');
fseek(fid, 5*4, 'bof');
fwrite(fid, [B, B + (npts-1)*delta, 0], 'float32');
fclose(fid);
end
