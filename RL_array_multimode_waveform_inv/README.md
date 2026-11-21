# Two-step Love + Rayleigh multimode array waveform inversion with MINEOS and idagrn6 (NoMelt, PA5)

A recovery test of a transdimensional Bayesian (rjMcMC) inversion for local
**Vsv(z), Vsh(z) and ξ(z) = (Vsh/Vsv)²** beneath the 16 NoMelt broadband
OBSs, from 15 real Pacific-rim earthquakes recorded in 2012 (Global CMT),
with idagrn6 synthetic seismograms of PA5 (Gaherty et al., 1996) as data
and MINEOS for the dispersion. Adapted from the `nomelt_LR4` demo package.

The inversion has two steps, so that the waveform inversion starts from a
reference model whose mode branches resemble the data's:

1. **Isotropic Rayleigh inversion** (`D1`): fundamental-mode Rayleigh
   phase velocities (20–100 s, MINEOS) are inverted for an isotropic Vs
   model (ξ = 1) with the linearized MINEOS inversion of
   `../run_MINEOS/c2_RLphv_inv_vsvh_vpvh_rho_eta.m`, starting from an
   isotropic PREM. The result has a lid and a low-velocity zone.
2. **Multimode waveform inversion** (`D2`–`D4`): Love + Rayleigh
   waveforms (17–112 s) are inverted for Vsv, Vsh and ξ with rjMcMC,
   using the step-1 model as the reference for the per-branch spectra,
   the MINEOS dispersion and kernels, and the time windows.

In both steps Vpv, Vph and density follow Vsv and Vsh: Vpv/Vsv, Vph/Vsh and
ρ/Vsv are held at the starting/reference model.

Runs in MATLAB (R2016b+) and Octave (7+); no toolboxes. Needs the MINEOS
and idagrn6 binaries (`../FORTRAN`, see `../FORTRAN/README`) and SAC's
`sacio` library for idagrn6 (SAC Community Edition,
https://github.com/EarthScope/sac-community, in `/opt/local/sac`).

## How to run

Run the drivers in order from this folder. Settings are in
`parameters_RL.m` (or define a struct `opt` with any of its fields before
calling a driver).

| Driver | What it does | Time (4 cores) |
|---|---|---|
| `D1_step1_rayleigh_isotropic_inv.m` | step 1: Rayleigh n = 0 phase velocities → isotropic Vs; writes the step-2 reference card `../run_MINEOS/CARDS/RLinv_step1.card` | ~5 min |
| `D2_step2_prep_mineos_idagrn.m` | MINEOS tables, dispersion and kernels; idagrn6 data (PA5, branches 0–15), reference spectra (step-1 model, branches 0–7) and PREM seismograms; cached in `prep_step2.mat` | ~10 min |
| `D3_step2_waveform_inv.m` | the waveform inversion with the plain (`opt.sampler = 'plain'`) or parallel-tempering (`'pt'`) sampler | hours |
| `D4_step2_plot_results.m` | figures and `summary.txt` in `figs_<sampler>/`; compares the two samplers if both were run (`figs_compare.png`) | 1 min |

```
>> D1_step1_rayleigh_isotropic_inv
>> D2_step2_prep_mineos_idagrn
>> D3_step2_waveform_inv                                  % plain rjMcMC (default)
>> opt.sampler = 'pt'; D3_step2_waveform_inv              % parallel tempering
$ CHAINS="1" octave-cli --eval D3_step2_waveform_inv &    % chains in parallel processes
$ ENSEMBLES="2" octave-cli --eval "opt.sampler='pt'; D3_step2_waveform_inv" &
```

With parallel processes, create the setup first (with `CHAINS="[]"` or
`ENSEMBLES="[]"`: Gauss–Newton fits of the nuisance parameters at the
reference and at the true dispersion), then start the chains/ensembles.

**Naming:** `D#_*.m` are the drivers, run in that order;
`parameters_RL.m` holds the settings for all of them. Everything else is in
`./functions`: functions, and two *drivee* scripts, `sub_setup_records.m`
and `sub_plot_results.m`, which the drivers call and which are not meant to
be run on their own.

## Step 1: isotropic Rayleigh inversion

- Data: MINEOS phase velocities of PA5, fundamental-mode Rayleigh,
  12 periods at 20–100 s, uncertainty 0.5%.
- Start: PREM below the PA5 Moho (12 km below sea level) with the PA5
  water, sediments and crust above it, made isotropic (Voigt average,
  η = 1; `functions/make_ref_card.m`, `functions/isotropic_card.m`).
- Inversion: `run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m` (the engine of
  `c2`) with Rayleigh data only, ξ damped to 1 (ε_ξ = 100), η damped to 1,
  smoothing broken only at the crustal discontinuities, MINEOS phase
  velocities and kernels recalculated after each of 4 iterations.

## Step 2: multimode waveform inversion

| | |
|---|---|
| events | 15 real 2012 Pacific-rim events, Mw 6.4–7.8, 33–86° (`events_nomelt_2012.cmt`) |
| data | idagrn6 seismograms of PA5, **branches 0–15**, transverse (Love) and vertical (Rayleigh) displacement, 17–112 s, normalized by each record's peak, 3% band-limited noise |
| modeled | **branches 0–7** |
| reference spectra | idagrn6_mask single-branch seismograms of the step-1 model |
| dispersion | MINEOS phase velocities of the step-1 model + MINEOS kernels × model perturbation, with Vp and ρ scaled to Vs (`functions/mineos_lin_setup.m`) |
| windows | from the step-1 model's MINEOS group velocities: overtone window from Δ/(max U of branches 1–7 + 0.6 km/s) to Δ/(max U₀ + 0.2), fundamental window from there to Δ/(min U₀ − 0.35) + 60 s |
| model | Vsv and ξ as transdimensional cubic B-splines from the Moho (7 km below the seafloor) to 350 km, bottom tied to the reference |
| nuisance | per record and branch: smooth spectral correction B_n(ω), frequency-dependent time shift τ_n(ω) and arrival angle δθ_n(ω) (see below) |
| noise | AR(1)-correlated, σ and r per window type, sampled hierarchically |

There is no path correction: the step-1 reference already predicts the
fundamental-mode Rayleigh path within the range the time shifts absorb.

### Nuisance parameters

For each record j and modeled branch n the predicted spectrum is

U_jn(ω) = B_n(ω) · U_ref,jn(ω) · exp(−iω τ_n(ω)) · exp(−iω [Δ/c_n(ω) − Δ/c_ref,n(ω)]) (+ the δθ_n(ω) term),

with B on 3 hat functions, τ_n(ω) on `opt.spec_Kt` Legendre terms in
frequency (default 3: constant + linear + quadratic; Kt = 1 is a constant
shift) and δθ_n(ω) on 2 Legendre terms. A frequency-dependent τ_n absorbs
what a constant shift cannot: the overtone path delays relative to the
reference are frequency dependent, because the 3-D/1-D structure the
reference misses is sampled differently at each period.

| option | values | meaning |
|---|---|---|
| `opt.nuisance` | `'marginal'` (default) | B, τ and δθ integrated out (Laplace approximation about a linearization point); only Vsv, ξ and the noise are sampled. Works with both samplers. |
| | `'sample'` | τ and δθ sampled by Metropolis with the exact waveforms, B integrated out (plain sampler only) |
| `opt.relin` | 1 (default) | the linearization points of τ, δθ and B of every state are moved to their conditional posterior mean (one Gauss–Newton step, clamped to ±2 s / ±3°) every `relin` iterations, i.e. before every proposal by default; larger values are cheaper and less accurate |
| `opt.spec_Kt` | 3 | number of τ_n(ω) terms |
| `opt.tau_prior`, `opt.dtau_prior` | 15 s, 10 s | prior std of the constant term and of the higher τ terms |

The initial linearization points come from a Gauss–Newton fit at the
reference dispersion with *continuation* (`functions/spec_refine.m`): a
constant τ_n first, then all terms. Freeing all terms at once lets
overlapping branches trade delays and ends in poorer minima (at the true
structure, median χ² of the Love fundamental windows 2.0 without and 1.13
with the continuation).

Why marginalize: sampling τ and δθ (9 × 8 × 30 extra parameters) mixes
slowly and couples to the structure, so different chains settle in
different (τ, structure) combinations; integrating them out leaves a
low-dimensional, tempering-friendly posterior. The Laplace approximation is
accurate when it is taken about the current conditional mode, hence the
relinearization.

RESULTS_PLACEHOLDER

## Files in `./functions`

| File | Purpose |
|---|---|
| `sub_setup_records.m`, `sub_plot_results.m` | drivee scripts: records, data and noise; figures and summary |
| `prep_mineos_idagrn.m` | MINEOS and idagrn6 runs of step 2 (cached) |
| `mineos_modes.m`, `idagrn_synth.m` | MINEOS tables/dispersion/kernels of any branch; idagrn6 full or single-branch seismograms |
| `make_ref_card.m`, `isotropic_card.m`, `card_at_depth.m` | cards: PREM under the PA5 crust, isotropic version, values at depth |
| `mineos_lin_setup.m`, `mineos_dispersion.m` | MINEOS dispersion of a model, Vp and ρ scaled to Vs |
| `run_pt_LR.m` (plain and PT, marginalized nuisance), `run_chain_LR.m` (sampled nuisance), `lr_common.m` | samplers |
| `spec_refine.m`, `nuis_relin.m`, `nuis_zero.m`, `spec_layout.m`, `spec_basis.m` | nuisance linearization points, Gauss–Newton / relinearization step, parameter layout and frequency bases |
| `event_G_spec.m`, `event_loglike_ar.m`, … | waveform operator and likelihood (from `nomelt_LR4`) |
| `fit_nuisance_LR.m` | best nuisance fit for a fixed structure (log L of the truth) |
| `update_reference.m` | optional further iteration: reference card from a posterior mean, MINEOS + idagrn6 rerun |
| `compare_runs.m` | posteriors of several runs in one figure |
| `read_cmtsolution.m`, `gc_dist_az.m`, `tukeywin_.m`, … | helpers |
