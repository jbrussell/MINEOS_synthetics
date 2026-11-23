# Fallbacks for missing toolbox functions

`butter`, `filtfilt`, `tukeywin` (signal toolbox/package) and `deg2km`,
`km2deg`, `distance` (mapping toolbox/package) are used by the ASWMS-Q
functions in `../aswms`. `parameters_RL.m` adds this folder at the END of
the path, so a toolbox's own function takes precedence when it is
installed; these only fill in when it is not (e.g. Octave without the
signal and mapping packages).
