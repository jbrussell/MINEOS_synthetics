Single-event line test 0: dispersion-only inversion (D7_dispersion_only_test.m) of the TRUE dispersion,
the upper bound on the xi resolution of the waveform tests line1/line2.

Geometry (shared with line1, line2): event 201210280304A (Haida Gwaii, Mw 7.8), 41 stations every 50 km on a
2000 km line along its great circle through the NoMelt centre (35.7-53.7 deg; stations_line_201210280304A_2000km.txt,
A00b_make_station_line.m).

Data: Love + Rayleigh phase velocities of branches 0-3 at the 10 dispersion periods (17-112 s), relative std 0.25%.
  lin    - the linearized dispersion of the true model about the step-1 reference (exactly representable)
  mineos - the exact MINEOS dispersion of the true model (adds the linearization error about the step-1 model;
           D7 uses one linearization, without the burn-in updates of D3)
Model space, priors (xi strength sampled, hierarchical), proposals and linearization: those of run line2
(setup_line2.mat: step-1 reference, current defaults), opt.disp_tag = '_line2'. 4 chains x 20000 iterations.

Result (summary.txt): xi comes out as a smooth gradient; the 160 km step is not resolved by branches 0-3 at 17-112 s
  xi band averages      lin                     mineos                  true
     7- 60 km   1.066 [1.056-1.082]     1.059 [1.053-1.069]     1.077
    60-160 km   1.050 [1.045-1.056]     1.043 [1.040-1.046]     1.064
   160-250 km   1.018 [1.012-1.024]     1.025 [1.019-1.028]     1.001
   250-350 km   0.999 [0.993-1.003]     1.009 [1.003-1.011]     1.001
