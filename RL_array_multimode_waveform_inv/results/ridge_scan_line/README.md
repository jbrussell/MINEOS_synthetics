# Ridge scan, line tests 2b / 2f

Does the line-test posterior have separate modes, or one ridge/basin? Log-likelihood of the line data
along straight lines between the posterior-mean models of the ensembles of tests 2b (random starts) and
2f (`refT1_TNjitter` start + joint dispersion term), and the true model. It is evaluated in the 2b setup:
- step-1 isotropic waveform reference;
- τ fixed at the truth;
- B free and Laplace-marginalized (`nuis_relin`);
- dθ about 0.

The noise σ and ρ are the medians of all posterior samples.

Point t of a line: C(t) = (1−t) C_A + t C_B for Love and Rayleigh, all branches. The dispersion is
linearized about the reference card, so this equals the line between the mean profiles.

Files:
- `ridge_scan.m`, `ridge_wave_L.m`: run from a work copy holding the setup and chain files of 2b and 2f.
- `ridge_plot.m` makes `ridge_scan.png`.
- `ridge_scan_log.txt` has the numbers. "barrier" = the lowest L for 0<t<1 minus the lower end point;
  a negative value means a barrier.

## Results

**No barrier on any of the 18 lines.** L is concave along every line and stays above the lower end
point. The posterior means of the ensembles lie on the flanks of one connected basin, not in separate
modes. The midpoints between pairs of ensembles fit better than either end.

**The truth is not at the top.** The mean model of 2f ensemble 1 fits about 235 better than the truth;
the maximum along truth → 2f E1 is at t ≈ 0.9.

| Model | ξ: lid | ξ: 160–250 km |
|---|---|---|
| Truth | 1.077 | 1.001 |
| 2f E1 | 1.039 | 1.037 |

This is with the posterior σ, which is about 3.4× the true noise. At the true noise the difference
would be about 10 times larger.

The bias of 2f is therefore in the likelihood of this setup, not a sampling failure. The test 2b
ensembles scatter around this basin without mixing; their mean models have lower L than the
combinations between them. Both points fit the explanation that the isotropic step-1 waveform
reference, with B free to absorb its error, moves the likelihood maximum. Test 2g (2b with
`update_ref 'burnin'`) tests this directly.
