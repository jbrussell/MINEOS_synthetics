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
| `run_recovery_tests.m`, `compare_runs.m` | the six recovery tests below and their comparison figure |
| `results/` | comparison figure, summaries of all runs, figures of runs 1 and 5 |
| `read_cmtsolution.m`, `gc_dist_az.m`, `card_at_depth.m`, `tukeywin_.m` | helpers |
| `events_nomelt_2012.cmt`, `stations_nomelt.txt` | events (Global CMT) and NoMelt stations |
| others | samplers, waveform operator and likelihood from `nomelt_LR4` |

`../run_MINEOS/parameter_FRECHET.m` gained an optional `MINEOS_NMODES`
override (number of branches in the mode table), alongside the existing
environment overrides used by `mineos_modes.m`.

## Recovery tests

Data in all runs: idagrn6 seismograms of PA5 (branches 0–7), 15 events
(14 Rayleigh + 15 Love records after quality control), 3% band-limited
noise. Plain rjMcMC, 4 chains (2 for runs 2 and 3), 4000 iterations (3000
for the 8-branch runs). `run_recovery_tests.m` repeats them;
`results/compare_runs.png` shows all six.

| run | reference | branches modeled | mean \|err\| Vsv / Vsh / ξ | truth in 5–95% (Vsv/Vsh/ξ) | σ/σ_true (L ovt, L fund, R ovt, R fund) | ξ at 7–60 / 60–160 / 160–250 / 250–350 km |
|---|---|---|---|---|---|---|
| true (PA5) | | | | | | 1.077 / 1.064 / 1.001 / 1.001 |
| **1** | **PREM** | **0–2** | 0.36 / 0.51 / 0.092 | 26 / 34 / 37 % | 1.50 1.51 1.51 1.23 | 1.115 / 1.077 / 1.101 / 0.908 |
| 2 | iterated (posterior of 1) | 0–2 | 0.25 / 0.23 / 0.055 | 14 / 13 / 41 % | 1.54 2.77 1.57 1.48 | 1.107 / 1.106 / 0.968 / 0.942 |
| 3 | PA5 | 0–2 | 0.16 / 0.22 / 0.055 | 32 / 28 / 9 % | 1.51 1.55 1.46 1.00 | 1.011 / 1.083 / 1.063 / 1.069 |
| 4 | PA5, *data with branches 0–2 only* | 0–2 | 0.065 / 0.087 / 0.026 | 89 / 71 / 55 % | 1.02 0.99 1.04 1.00 | 1.044 / 1.063 / 1.033 / 1.033 |
| **5** | **PA5** | **0–7** | **0.064 / 0.066 / 0.019** | **85 / 78 / 89 %** | 1.00 0.98 1.00 1.00 | 1.060 / 1.058 / 1.022 / 1.030 |
| 6 | PREM | 0–7 | 0.25 / 0.32 / 0.055 | 49 / 40 / 42 % | 1.05 1.15 1.12 1.66 | 1.123 / 1.079 / 0.991 / 0.915 |

**With the PREM reference (run 1), Vsv, Vsh and ξ are not recovered.** The
posterior is narrow and wrong: ξ ≈ 1.12 in the lid and 1.10 at 160–250 km,
Vsv runs into the prior bounds below 150 km, all chains use the maximum
number of control points, and the noise level is overestimated by 50%: the
model is absorbing errors of the forward model, not resolving structure.
Only run 5 — a reference close to the truth **and** every branch that is in
the data in the model — recovers PA5 within its uncertainties, and even
there ξ is only resolved as a smooth profile (ξ ≈ 1.06 from the Moho to
~150 km, ~1.02–1.03 below); the step from 1.07 to 1.00 at the L
discontinuity (~160 km) is not resolved.

## Why the radial anisotropy is not resolved

In order of importance (each was isolated with the runs above or a direct
MINEOS/idagrn6 check):

1. **Mode branches in the data that are not in the model.** With real
   moment tensors, idagrn6 puts substantial energy in branches 3–7: 52–62%
   of the overtone windows and 16% of the Love fundamental window
   (Guerrero event, PA5 reference; 0% of the Rayleigh fundamental window).
   The Love and Rayleigh overtone group velocities overlap, and the Love
   fundamental (4.38–4.51 km/s) overlaps the Love overtones (4.40–5.60
   km/s), so windowing cannot separate them. The sampler then distorts Vsv,
   Vsh and ξ to fit that energy: with a *perfect* reference (run 3) the
   chains fit the data better than the true structure (log L 295,900 vs
   281,200), σ is 1.5× too high in exactly the windows with higher-mode
   energy, and ξ is 1.01 in the lid and 1.07 below 160 km. Removing branches
   3–7 from the data (run 4) or modeling them (run 5) fixes this. The
   original package simulated only branches 0–2.

2. **Love-wave branch osculation in PA5.** PA5's fast lid (Vsh 4.83 km/s)
   over a strong LVZ makes the MINEOS Love branches 0 and 1 swap character
   at 15–30 s: in PA5 branch 0 is an LVZ-channel mode (Vsh sensitivity
   centred at 150–163 km) and branch 1 is the lid-guided "fundamental"
   (92–105 km) and carries 54% of the Love fundamental-window energy, while
   in PREM branch 0 is the ordinary fundamental (35–80 km). The PREM-based
   reference spectra and dispersion of each branch therefore describe the
   wrong wave, and the inversion makes the Love branch-0 phase velocity
   0.15–0.20 km/s too fast (`results/fig3_dispersion_fits_s1.png`), i.e.
   Vsh, and ξ, far too high in the lid. The same osculation makes the Love
   dispersion strongly nonlinear: linearizing MINEOS about PREM
   overpredicts the PA5 Love fundamental by 0.14 km/s at 16 s even though
   the kernels pass a finite-difference test (a hybrid card, PREM with PA5
   Vsv/Vsh at 7–350 km, reproduces the PA5 Love dispersion to <0.004 km/s).
   One iteration of the reference (run 2) and modeling all branches (run 6)
   both leave this bias.

3. **Little information below ~150 km.** The overtone windows are at or
   below the noise (rms 0.017–0.029 of the record peak vs 0.03 noise) with
   real excitation (the original analytic sources had overtone amplitudes
   of 45–60% of the fundamental), and the array aperture (~500 km) gives
   weak constraints on long-period phase velocities: a 0.13 km/s error of
   the 100-s Rayleigh wave changes the phase across the array by only ~0.1
   rad. ξ below ~150 km is therefore set largely by the prior and the
   modeling errors (Vsv hits the 5.0 km/s prior bound in runs 1–3 and 6).

4. **The path to the array.** With PREM as reference, the fundamental
   Rayleigh wave arrives up to 107 s (several cycles) earlier in the data
   than in the reference spectra — far beyond what the time shifts and the
   smooth spectral correction of the nuisance model can absorb.
   `path_precorrect.m` measures the array-average delay from the data
   (unwrapped stacked cross-spectrum) and folds it into Uref: 1–5 s error
   for the records that pass QC (14 of 15). Love and overtone path delays
   (up to ~30 s) could not be measured reliably (overlapping windows,
   higher branches) and are left to the per-mode time shifts, which is one
   reason runs 1, 2 and 6 misfit (σ 1.2–1.7× too high).

5. **Structure held fixed at the reference.** Vp, density, η and the
   structure below 350 km are PREM in runs 1, 2 and 6: this alone shifts
   the long-period Rayleigh phase velocities by ~0.026 km/s at 114 s.

What would make ξ resolvable here: model every branch that carries energy
in the windows (run 5 vs 3), start from — or iterate to — a reference
whose branch structure matches the data (e.g. a lid/LVZ oceanic model
rather than PREM; relinearize MINEOS and recompute the idagrn6 reference
spectra at each stage), and recompute the dispersion with MINEOS rather
than with fixed kernels where the branches osculate. Better convergence
would need parallel tempering (`demo_LR_pt.m`, not run here): the median
R-hat of the plain chains is 1.2–4.7 across the runs.

## Caveats

- No site terms, azimuthal anisotropy or lateral heterogeneity; Q is the
  reference `.qmod` in both data and reference spectra.
- The source time function and the CMT centroid time are not applied
  (identical in data and reference spectra).
- Only tested in Octave 8.4 (not MATLAB), with MINEOS and idagrn6
  compiled with gfortran on Linux and SAC Community Edition 103.0.
