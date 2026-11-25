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
   isotropic PA5 without its lid/LVZ that has PA5's Vp/Vs, ρ/Vs and η (or
   an isotropic PREM). The result has a lid and a low-velocity zone.
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
| `D2_step2_prep_mineos_idagrn.m` | MINEOS tables, dispersion and kernels; idagrn6 data (PA5, branches 0–15), reference spectra (step-1 model, branches 0–7) and seismograms of the step-1 starting model (all branches, for the figures); cached in `prep_step2.mat`; record sections of the full mode-sum synthetics (PA5 and the step-1 start, T and Z, with the window bounds) in `figs_step2/` | ~10 min |
| `D3_step2_waveform_inv.m` | the waveform inversion with the parallel-tempering (`opt.sampler = 'pt'`, default) or plain (`'plain'`) sampler | hours |
| `D4_step2_plot_results.m` | figures and `summary.txt` in `figs_<sampler>/`; compares the two samplers if both were run (`figs_compare.png`) | 1 min |
| `D5_step2_check_linearization.m` | full MINEOS vs the linearized dispersion at the posterior mean and at the (projected) truth: linearization error vs what the parameterization cannot represent (`d5_linearization.txt`, `fig10_linearization.png`) | ~5 min |
| `D6_step2_update_reference_and_rerun.m` | optional outer iteration: new reference from the posterior mean of a run (MINEOS + idagrn6), then D3/D4 again about it (tag `<tag>_it2`) | ~10 min + D3 |

```
>> D1_step1_rayleigh_isotropic_inv
>> D2_step2_prep_mineos_idagrn
>> D3_step2_waveform_inv                                  % parallel tempering (default)
>> opt.sampler = 'plain'; D3_step2_waveform_inv           % plain rjMcMC chains
$ ENSEMBLES="1" octave-cli --eval D3_step2_waveform_inv & % ensembles in parallel processes
$ CHAINS="2" octave-cli --eval "opt.sampler='plain'; D3_step2_waveform_inv" &
```

With parallel processes, create the setup first (with `ENSEMBLES="[]"`, or
`CHAINS="[]"` for the plain sampler: Gauss–Newton fits of the nuisance parameters at the
reference and at the true dispersion), then start the chains/ensembles.
Start all of them at once: with the end-of-burn-in reference update (the
default) each process waits after its burn-in until all chains/ensembles
have finished theirs; the first process to see all burn-in files runs
MINEOS + idagrn6 for the new reference (~10 min) and the others then
continue with it. MINEOS/idagrn6 runs are serialized through the lock
directory `../run_MINEOS/.prep_lock` (remove it by hand if a process died
while holding it).

**Naming:** `D#_*.m` are the drivers, run in that order;
`parameters_RL.m` holds the settings for all of them. Everything else is in
`./functions`: functions, and two *drivee* scripts, `sub_setup_records.m`
and `sub_plot_results.m`, which the drivers call and which are not meant to
be run on their own.

## Step 1: isotropic Rayleigh inversion

- Data: MINEOS phase velocities of PA5, fundamental-mode Rayleigh,
  12 periods at 20–100 s, uncertainty 0.5%.
- Start (`opt.s1_start_model`):
  - `'pa5_iso'` (default): an isotropic PA5 without its lid/LVZ structure,
    built as the starting model of
    `SURF96_Rayl_1Dinv/test_kernel_inversion_vsvh_vpvh_rho_eta_PA5_twostage.m`
    (`functions/pa5_iso_ref_card.m`). Between the Moho (12 km below sea
    level) and `opt.s1_zlin` = 250 km, Vsv = Vsh is a straight line from
    the isotropic (Voigt) Vs of PA5 at the Moho (4.72 km/s) to that at
    250 km (4.47 km/s): no lid, no LVZ, no Vs discontinuities. There Vpv,
    Vph and ρ are this Vs times PA5's Vpv/Vsv, Vph/Vsh and ρ/Vsv, so these
    ratios, η and Q (`opt.ref_qmod`) are PA5's at every depth; PA5 itself
    above the Moho and below 250 km. Because Vs has no drop at 68.5 km while
    PA5's Vp/Vs and ρ/Vs jump up there, the starting Vp and ρ jump up at
    68.5 km (ρ by ~0.2 g/cm³) and are higher than PA5's down to 250 km
    (`figs_step2/pa5_iso_reference.png`).
  - `'prem'`: PREM below the PA5 Moho with the PA5 water, sediments and
    crust above it, made isotropic (Voigt average, η = 1;
    `functions/make_ref_card.m`, `functions/isotropic_card.m`).
- Inversion: `run_mineos_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m` (the engine of
  `c2`) with Rayleigh data only, ξ damped to 1 (ε_ξ = 100), η damped to the start,
  smoothing broken only at the crustal discontinuities, MINEOS phase
  velocities and kernels recalculated after each of 4 iterations. With
  `'pa5_iso'` the starting (PA5) Vp/Vs, ρ/Vs and η are imposed exactly on
  the result, so the step-2 reference (the step-1 model, which recovers a
  lid and an LVZ) has PA5's Vp/Vs, ρ/Vs and η: with the true Vsv and Vsh it
  is PA5, and η, the Vp/ρ scaling and the structure below 350 km are no
  longer sources of error in step 2. From the `'pa5_iso'` start, step 1
  fits the data (reduced χ² 47.7 → 0.01) with a smooth lid (Vs ≈ 4.66 km/s
  at ~35 km) and LVZ (≈ 4.30 km/s at ~110 km): `results/step1_pa5iso_start.png`.

## Step 2: multimode waveform inversion

A full description with equations and notation: `docs/step2_waveform_inversion.tex` (PDF alongside).

| | |
|---|---|
| events | 15 real 2012 Pacific-rim events, Mw 6.4–7.8, 33–86° (`events_nomelt_2012.cmt`) |
| data | idagrn6 seismograms of PA5, **branches 0–15**, transverse (Love) and vertical (Rayleigh) displacement, 17–112 s, normalized by each record's peak, 3% band-limited noise |
| modeled | **branches 0–7** |
| reference spectra | idagrn6_mask single-branch seismograms of the step-1 model |
| dispersion | MINEOS phase velocities of the reference + Σ over Vsv, Vsh, Vpv, Vph, ρ, η of each parameter's MINEOS kernel × its perturbation; Vpv, Vph and ρ follow Vsv, Vsh with the reference ratios, η is the reference's unless `opt.eta_free` (`functions/mineos_lin_setup.m`, `mineos_lin_disp.m`) |
| reference update | by default the reference (spectra, phase velocities, kernels) is recomputed once at the end of burn-in from the mean model of all chains, then held fixed (`opt.update_ref`) |
| windows | from the step-1 model's MINEOS group velocities: by default (`opt.windows = 'single'`) one surface-wave window per wave type from Δ/U_start to Δ/(min U₀ − 0.35) + 60 s, containing the overtones and the fundamental mode, with one noise level and correlation per wave type; U_start = max U of branches 1–7 + 0.6 km/s, or `opt.U_start` (km/s; e.g. 6 to leave out the early, body-wave-like energy). `'split'`: an overtone window (to Δ/(max U₀ + 0.2)) and a fundamental-mode window with separate noise parameters. The windows stay those of the step-1 reference when the reference is updated |
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
| `opt.relin` | 1 (default) | the linearization points of τ, δθ and B of every state are moved to their conditional posterior mean (one Gauss–Newton step, clamped to ±2 s / ±3°) every `relin` iterations, i.e. before every proposal by default; larger values are cheaper and less accurate. Cost: one rebuild of the waveform derivatives and one likelihood per record, about the cost of a structural proposal, so `relin` = 1 roughly doubles the time per iteration (the posterior mean used for the update is cached from the last likelihood evaluation) |
| `opt.spec_Kt` | 3 | number of τ_n(ω) terms |
| `opt.windows` | `'single'` (default) / `'split'` | one surface-wave window per wave type (one σ, r each), or separate overtone and fundamental-mode windows (σ, r each); separate noise levels down-weight whichever part the model fits worse (e.g. the overtone window with unmodeled branches 8–15) |
| `opt.U_start` | `[]` (default) or km/s | start of the (overtone part of the) window; `[]`: max group velocity of branches 1–7 + 0.6 km/s (8.3 km/s Love, 9.2 km/s Rayleigh for the step-1 model) |
| `opt.update_ref` | `'burnin'` (default) / `'none'` | recompute the reference model (idagrn6_mask spectra, MINEOS phase velocities and kernels) once at the end of burn-in, from the mean Vsv, Vsh (and η) of the last `opt.upd_frac` (0.2) of the burn-in of all chains, and hold it fixed afterwards; `'none'`: the step-1 reference throughout. The data, windows and normalizations never change |
| `opt.eta_free` | `false` (default) / `true` | η held at the reference model, or a third transdimensional B-spline profile (k_η ≤ k_V, prior U[`opt.etamin`, `opt.etamax`] = U[0.80, 1.05], bottom tied to the reference, a fraction `opt.p_eta` = 0.2 of the moves). Only Rayleigh waves are sensitive to η |
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

## Results (PA5 recovery test)

Figures and summaries in `results/` (`plain/`: the run below;
`baseline_plain/`: the earlier constant-τ run with sampled nuisance
parameters). Both runs: 4 independent chains × 3000 iterations (1500 burn-in),
split windows, step-1 reference throughout, η at the reference. The runs
with the new defaults (single window, end-of-burn-in reference update) and
parallel tempering are in the next subsection.

| | baseline (τ_n constant, τ and δθ sampled) | marginalized nuisance, τ_n(ω) (Kt = 3), relinearized every iteration |
|---|---|---|
| R-hat median Vsv / Vsh / ξ | 327 / 8.2 / 3.3 | **4.3 / 4.2 / 2.0** |
| log L of the chains (each vs log L of PA5 with fitted nuisance parameters) | 386k–392k vs 368k | 400.4k–400.8k vs 385k |
| truth inside 5–95%: Vsv / Vsh / ξ | 65 / 50 / 50 % | 73 / 53 / 38 % |
| mean abs. error Vsv / Vsh (km/s), ξ | 0.19 / 0.27, 0.051 | 0.15 / 0.21, 0.072 |
| σ / σ_true (L ovt, L fund, R ovt, R fund) | 1.47 1.57 1.19 1.42 | 1.41 1.14 1.15 1.02 |
| ξ 7–60 / 60–160 / 160–250 / 250–350 km (true 1.08 / 1.06 / 1.00 / 1.00) | 1.14 / 1.09 / 1.02 / 0.92 | 1.01 / 1.14 / 0.98 / 0.91 |
| time per iteration (4 chains on 4 cores) | ~1.8 s | ~3.3 s (4.3 s with other jobs sharing the CPU); relinearization is about half of it |

What the second run shows:

- **Convergence.** Marginalizing τ_n(ω), δθ_n(ω) and B_n(ω) with
  relinearization improved the mixing a lot (the four chains now reach the
  same log-likelihood, 400.4k–400.8k, and R-hat drops by one to two orders
  of magnitude), but it is still not converged (R-hat 2–4).
- **Fit vs truth.** The waveforms are fit closely
  (`results/plain/fig4_waveform_fits_event1.png`: data, PREM, step-1
  reference and final model) and the noise levels are close to the true
  ones, yet the chains reach a log-likelihood ~15k **above** that of the
  true structure. Within this forward model the truth is not the best fit:
  the structure absorbs forward-model errors. The posterior Vsv oscillates
  with depth and ξ sits on the prior bounds (~1.2 at 60–160 km, 0.90 below
  250 km), and the truth is outside much of the 5–95% range.
- **Why (D5, `results/plain/d5_linearization.txt`).** Full MINEOS vs the
  linearized dispersion, rms per branch (km/s):

  | | Love n = 0 | Love n ≥ 1 | Rayleigh n = 0 | Rayleigh n ≥ 1 |
  |---|---|---|---|---|
  | linearization error at the posterior mean | 0.050 | 0.008–0.040 | 0.014 | 0.013–0.030 |
  | linearization error at the (projected) truth | 0.010 | 0.004–0.040 | 0.020 | 0.006–0.008 |
  | not representable (projected truth − PA5) | 0.001 | 0.017–0.030 | 0.027 | 0.051–0.065 |
  | … after adding PA5's η | 0.001 | 0.017–0.030 | 0.025 | 0.022–0.042 |
  | … and PA5's Vp and ρ (no ratio scaling) | 0.003 | 0.011–0.023 | 0.003 | 0.009–0.017 |
  | … and PA5 below 350 km | 0.002 | 0.003–0.007 | 0.003 | 0.002–0.005 |

  The largest errors are not linearization errors but structure the model
  cannot represent: **η** for the Rayleigh overtones, the **Vp/Vs and ρ/Vs
  scaling** for the Rayleigh fundamental mode, and **PA5 below 350 km** for
  the overtones of both. These are comparable to the phase-velocity signal
  of the radial anisotropy itself, so the sampler bends Vsv and ξ to
  compensate (e.g. the deep ξ at the lower bound).
- **D6 not run.** An outer iteration would remove the linearization error
  at the posterior mean (up to 0.05 km/s for the Love fundamental mode) but
  not the unrepresentable part, and the posterior is not converged. More
  useful next steps: `opt.eta_free = true` (η is the largest missing term
  for the Rayleigh overtones), a deeper model (`cfg.zmax`) and, for the
  Rayleigh fundamental mode, Vp and ρ not tied to Vs. The new defaults
  (single surface-wave window, end-of-burn-in reference update) also apply
  to the next run.


### Runs with the new defaults: plain vs parallel tempering

`results/plainnew/` (4 plain chains × 3000 iterations, 1500 burn-in) and
`results/ptnew/` (2 parallel-tempering ensembles of K = 4 temperatures up to
T = 10, 1200 iterations, 600 burn-in), both with the single surface-wave
window, the end-of-burn-in reference update and η fixed at the reference.
`results/compare_plainnew_ptnew.png` compares the two posteriors. The
log-likelihoods are not comparable with the split-window runs above: with
the single window, the true structure has log L = 387.2k (385.1k with split
windows).

| | plain, split windows (above) | plain, new defaults | parallel tempering, new defaults |
|---|---|---|---|
| saved samples | 1200 | 1200 | 240 (2 ensembles, T = 1) |
| R-hat median Vsv / Vsh / ξ | 4.3 / 4.2 / 2.0 | 4.6 / 3.7 / 3.4 | **1.27 / 1.36 / 1.15** |
| R-hat max Vsv / Vsh / ξ | 162 / 12.6 / 9.0 | 362 / 29.7 / 28.9 | 10.5 / 6.8 / 2.4 |
| log L of the chains vs log L of PA5 | 400.4k–400.8k vs 385.1k | 397.2k–399.0k vs 387.2k | 398.6k–398.9k vs 387.2k |
| truth inside 5–95%: Vsv / Vsh / ξ | 73 / 53 / 38 % | 74 / 52 / 58 % | 70 / 54 / 73 % |
| mean abs. error Vsv / Vsh (km/s), ξ | 0.15 / 0.21, 0.072 | 0.11 / 0.19, 0.058 | 0.10 / 0.19, 0.062 |
| σ / σ_true (L ovt, L fund, R ovt, R fund) | 1.41 1.14 1.15 1.02 | 1.34 1.34 1.08 1.08 | 1.34 1.34 1.08 1.08 |
| ξ 7–60 / 60–160 / 160–250 / 250–350 km (true 1.08 / 1.06 / 1.00 / 1.00) | 1.01 / 1.14 / 0.98 / 0.91 | 1.04 / 1.11 / 0.97 / 0.91 | 1.00 / 1.10 / 0.99 / 0.92 |
| end-of-burn-in reference update | – | rejected (dL −5.6k … −8.8k) | rejected (dL −17.0k, −18.8k) |
| time per iteration (6 processes on 4 cores) | – | ~6 s | ~15 s (4 temperatures) |

- **Parallel tempering mixes much better.** Its two ensembles agree
  (R-hat median 1.2–1.4, and both reach log L ≈ 398.7k), while the plain
  chains still sit in different modes (R-hat 3–5; one chain has only 5
  interior control points). Its posterior is also the best calibrated for ξ
  (73% of the true profile inside the 5–95% range). But only 2 ensembles of
  600 post-burn-in iterations were run, so its R-hat is a rough estimate;
  swap acceptance between neighbouring temperatures was low (0–6%).
- **Single window.** With the single window the fundamental and overtone
  noise levels are tied per wave type (σ/σ_true 1.34 for Love, 1.08 for
  Rayleigh). The models and errors are similar to the split-window run, a
  little closer to the truth (mean error Vsv 0.10–0.11 vs 0.15 km/s).
- **Both runs still fit the data ~11k better than the true structure**, with
  ξ too high at 60–160 km and too low (0.91–0.92) below 250 km. D5 on both
  runs (`results/*new/d5_linearization.txt`) gives the same split of the
  phase-velocity errors as for the first run: the unrepresentable part
  (η, Vp/ρ scaling, PA5 below 350 km; 0.02–0.065 km/s for the Rayleigh
  overtones) is larger than the linearization error at the posterior mean
  (0.009–0.06 km/s).
- **Why the reference update was rejected** (`results/refupd_diagnostic.txt`).
  The end-of-burn-in states fit 6–19k worse with the new reference. About
  1.5–3k of that is a bias of the check itself (the same states re-scored
  with the *old* reference through the same path lose that much, because
  the nuisance linearization point is refitted from scratch); the rest is
  real. With the new reference, even the true structure fits worse
  (log L 374.7k–377.5k instead of 387.2k), and the mean model is further
  from PA5 than the step-1 reference is (rms Vsh 0.16–0.18 vs 0.10 km/s).
  The late burn-in mean is a model that partly compensates for the
  unrepresentable structure, and reference spectra computed from it are a
  worse description of the data than those of the step-1 model. The
  safeguard worked as intended: both runs sampled with the step-1 reference.
- **D6 not run**, for the same reason: it would rebuild the reference about
  the same kind of compensating mean model. The useful next steps remain
  `opt.eta_free = true`, a deeper model (`cfg.zmax`), and Vp and ρ not tied
  to Vs, then more and longer parallel-tempering ensembles.

## Files in `./functions`

| File | Purpose |
|---|---|
| `sub_setup_records.m`, `sub_plot_results.m` | drivee scripts: records, data and noise; figures and summary |
| `prep_mineos_idagrn.m` | MINEOS and idagrn6 runs of step 2 (cached) |
| `mineos_modes.m`, `idagrn_synth.m` | MINEOS tables/dispersion/kernels of any branch; idagrn6 full or single-branch seismograms |
| `make_ref_card.m`, `isotropic_card.m`, `card_at_depth.m` | cards: PREM under the PA5 crust, isotropic version, values at depth |
| `mineos_lin_setup.m`, `mineos_lin_disp.m`, `mineos_dispersion.m` | linearized MINEOS dispersion: per-parameter kernels; c for given Vsv, Vsh, η node profiles (Vp, ρ scaled to Vs); c of a model |
| `apply_reference.m`, `ref_card_update.m`, `mineos_lock.m` | install a reference model (spectra, kernels) with the data and windows unchanged; card from Vsv/Vsh(/η) profiles with Vp and ρ scaled; serialize MINEOS runs |
| `run_pt_LR.m` (plain and PT, marginalized nuisance), `run_chain_LR.m` (sampled nuisance), `lr_common.m` | samplers |
| `spec_refine.m`, `nuis_relin.m`, `nuis_zero.m`, `spec_layout.m`, `spec_basis.m` | nuisance linearization points, Gauss–Newton / relinearization step, parameter layout and frequency bases |
| `event_G_spec.m`, `event_loglike_ar.m`, … | waveform operator and likelihood (from `nomelt_LR4`) |
| `fit_nuisance_LR.m` | best nuisance fit for a fixed structure (log L of the truth) |
| `update_reference.m` | new reference card from Vsv/Vsh(/η) profiles (e.g. a posterior mean), MINEOS + idagrn6 rerun (end-of-burn-in update, D6) |
| `compare_runs.m` | posteriors of several runs in one figure |
| `plot_record_sections.m` | record sections of the full mode-sum synthetics (D2) |
| `read_cmtsolution.m`, `gc_dist_az.m`, `tukeywin_.m`, … | helpers |
