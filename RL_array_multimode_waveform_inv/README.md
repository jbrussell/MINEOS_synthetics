# Joint Love + Rayleigh multimode array waveform inversion with MINEOS and idagrn6 (NoMelt, PA5)

A recovery test of a transdimensional Bayesian (rjMcMC) inversion for local
**Vsv(z), Vsh(z) and ξ(z) = (Vsh/Vsv)²** beneath the 16 NoMelt broadband
OBSs, adapted from the `nomelt_LR4` demo package so that

- the **synthetic seismograms come from idagrn6** (full normal-mode
  summation from MINEOS mode tables, real moment tensors, real geometry),
- the **dispersion comes from MINEOS** (phase velocities and Vsv/Vsh
  kernels of branches n = 0–2),
- the data are **15 real Pacific-rim earthquakes** recorded while NoMelt was
  deployed (2012), and
- the **input (reference) model is PREM**; the true model is PA5
  (Gaherty et al., 1996).

It runs in MATLAB (R2016b+) and Octave (7+); no toolboxes are needed. The
MINEOS and idagrn6 programs (`../FORTRAN`, see `../FORTRAN/README`) must be
compiled, and idagrn6 needs SAC's `sacio` library (SAC Community Edition,
https://github.com/EarthScope/sac-community, installed in `/opt/local/sac`).

```
>> demo_LR                                         % stage 1: PREM reference
>> update_reference('_plain', 'prep_s2.mat', 'RLinv_ref_s2');
>> opt = struct('tag', '_s2', 'outdir', 'figs_s2', 'prepfile', 'prep_s2.mat'); demo_LR
$ CHAINS="1 2" octave-cli --eval demo_LR &         % chains in parallel processes
>> demo_LR_pt                                      % parallel tempering
```

The first call runs MINEOS and idagrn6 (`prep_mineos_idagrn.m`, about 2
minutes) and caches the results in `prep_mineos.mat`.

## What changed relative to `nomelt_LR4`

| | nomelt_LR4 | this package |
|---|---|---|
| events | 16 synthetic back-azimuths, random distances | 15 real 2012 events, Global CMT (`events_nomelt_2012.cmt`) |
| data | analytic plane-wave modes n = 0–2 with PA5 dispersion | idagrn6 seismograms of PA5, **branches 0–7**, T and Z displacement |
| reference spectra Uref | analytic source term | idagrn6_mask single-branch seismograms (n = 0–2) of the reference model |
| path to the array | reference model in data and Uref | **PA5 in the data, PREM in Uref** (corrected from the data, see below) |
| dispersion | Haskell/propagator solvers; Rayleigh isotropic | MINEOS: c = c_ref + K_vsv·δVsv + K_vsh·δVsh (radially anisotropic Rayleigh) |
| reference model | PA5 background | PREM below the Moho, PA5 water/sediments/crust above it |
| fixed Vp, ρ, η, deep structure | PA5 (the truth) | reference (PREM) |
| noise | 3% of the record peak | 3% of each record's peak (records normalized) |

### Files

| File | Purpose |
|---|---|
| `demo_LR.m`, `demo_LR_pt.m` | drivers (plain rjMcMC, parallel tempering) |
| `setup_LR.m` | configuration, records, data, noise, path correction |
| `prep_mineos_idagrn.m` | runs MINEOS and idagrn6 once (cached) |
| `mineos_tools/make_ref_card.m` | reference card: PREM mantle under the PA5 water/sediments/crust |
| `mineos_tools/mineos_modes.m` | MINEOS mode tables, phase/group velocities and kernels of any branch |
| `mineos_tools/idagrn_synth.m` | idagrn6 full or single-branch synthetics for a list of events |
| `mineos_lin_setup.m`, `mineos_dispersion.m` | MINEOS dispersion of a model (reference + kernels) |
| `update_reference.m` | iterated reference: new reference card from a posterior mean, MINEOS + idagrn6 rerun |
| `path_precorrect.m` | data-derived source-to-array path correction of Uref |
| `fit_nuisance_LR.m` | best nuisance fit for a fixed structure (log L of the truth, misfit checks) |
| `read_cmtsolution.m`, `gc_dist_az.m`, `card_at_depth.m`, `tukeywin_.m` | helpers |
| `events_nomelt_2012.cmt`, `stations_nomelt.txt` | events (Global CMT) and NoMelt stations |
| others | samplers, waveform operator and likelihood from `nomelt_LR4` |

`../run_MINEOS/parameter_FRECHET.m` gained an optional `MINEOS_NMODES`
override (number of branches in the mode table), alongside the existing
environment overrides used by `mineos_modes.m`.

RESULTS_PLACEHOLDER
