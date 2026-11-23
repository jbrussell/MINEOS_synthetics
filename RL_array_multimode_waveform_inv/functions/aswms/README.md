# ASWMS-Q functions

Copied from [ASWMS-Q](https://github.com/jbrussell/ASWMS-Q) at commit
`7255f9d` (2026-03-09). They are used by the data-preparation scripts
`A1_sac2eventmat.m`, `A2_cleanup_events.m`, `A3_run_autowinpick.m` and
`A3b_run_exam_winpara.m`. Real data and the synthetic "data" from
`D0_synth_make_sac_data.m` go through the same steps.

| File | Used by |
|---|---|
| `calc_SNR.m`, `readsac.m`, `vdist.m` | A1 (SAC files to `eventmat/*.mat`) |
| `auto_win_select.m`, `groupv_fit.m`, `build_gaus_filter.m`, `Is_inrange.m`, `plot_win_select.m` | A3 (automatic window picking) |
| `exam_winpara_phase.m`, `phasedb.mat` | A3b (interactive window review; MATLAB) |
| `setup_ErrorCode.m` | A2, A3 |

Changes from ASWMS-Q:
- `setup_parameters.m` is this package's own. It returns the
  `parameters` structure installed by `functions/aswms_setup.m`, which is
  built from `opt` (`parameters_RL.m`) by `functions/aswms_parameters.m`.
- `exam_winpara_phase.m` loads `phasedb.mat` from this folder instead
  of `./functions/`.

The other files are unchanged. Under Octave without the signal and
mapping packages, `functions/octave_shims` supplies `butter`, `filtfilt`,
`tukeywin`, `deg2km`, `km2deg` and `distance`.
