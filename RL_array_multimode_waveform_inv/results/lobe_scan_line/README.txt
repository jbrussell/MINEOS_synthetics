Side-lobe scan of the single-event line data (201210280304A, 41 stations on a 2000 km line along the great circle
through NoMelt, noise 0.005): logL vs a uniform relative shift (-3..+3 %, all periods) of ONE branch's phase
velocity, the other branches at the true dispersion (lobe_scan.m in the work directory; lobe_plot.m plots it).
  oracle : setup_line1.mat  (reference = true model, B/tau/dth fixed at the truth)
  2b     : setup_line2b.mat (step-1 isotropic reference, tau fixed at the truth, B free, dth held at ~0)
  with sigma = sigma_true and 3 sigma_true (about the posterior noise level), AR(1) rho from the data.
Result: every branch has a SINGLE maximum within +-3 % (no cycle-skip side lobes along these directions): the
broad band (17-112 s) and the moveout across the line remove the 1-cycle lobes of a narrow band. Love branches and
the Rayleigh fundamental are sharp (logL -3000 at +-0.1 % (R0) to +-1 % (L0) with sigma_true); the Rayleigh overtones
are nearly flat (R2, R3: < 150 logL over +-3 %). With the step-1 reference and free B (2b) the maxima are biased from
the truth: L0 +0.23 % (+46), L3 +0.23 % (+127), R0 +0.03 % (+346), R2/R3 up to +-0.7 % (flat).
