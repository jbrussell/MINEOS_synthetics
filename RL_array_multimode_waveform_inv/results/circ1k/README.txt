Project circ1k
  run tag _circ1k, sampler pt, archived 2026-11-28 05:29 from /tmp/claude-0/-home-user-MINEOS-synthetics/578903d7-e992-5265-9feb-4720fd8705a2/scratchpad/work11/RL_array_multimode_waveform_inv

Key settings (all of them in settings.txt and parameters/opt_used.mat):
  trueid               = 'pa5_5km'
  aswms_dir            = 'ASWMS_circ1k/'
  nbranch_data         = 4
  nm                   = 4
  sampler              = 'pt'
  nens                 = 4
  K                    = 6
  Tmax                 = 20
  pt_nit               = 3000
  pt_burn              = 750
  nchains              = 4
  nit                  = 3000
  burn                 = 1500
  noise_model          = 'ar1'
  noise_levels         = 'pre'
  windows              = 'auto'
  nuisance             = 'marginal'
  spec_K               = 4
  spec_Ka              = 2  (mc)
  spec_Kt              = 3
  update_ref           = 'lin'
  lin_updates          = 3
  smooth_sig           = [0.0002 0.0002]
  smooth_hier          = false
  vsv_ref_sig          = 'step1'
  xi_ref_hier          = true
  synth_arrival_shift  = []
  noise                = 0.03
  s1_cstd              = 0.0025
  prepfile             = 'prep_circ1k.mat'

Synthetic data (A0, ASWMS_circ1k/synth_truth.mat): true model 'pa5_5km', data branches 4, noise 0.03 x peak, seed 3, station noise spread 0

Results: figures/summary.txt (statistics), figures/posterior_summary_circ1k.mat
  (pooled posterior and its statistics), figures/*.png, step2/chain_circ1k_XX.mat (samples).

Folders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/, data/ (waveform data and model cards)

Rerun (from RL_array_multimode_waveform_inv/):
  >> opt = load_project('circ1k');          % the exact options of this run
  >> opt = load_project('circ1k', true);    % the same, with the archived waveform data (data/)
  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,
  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data
  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.
  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.

--- Repository copy: figures, summary, settings and parameters only (chains, prep and data files stay in the
working directory; *.mat files are not tracked). ---
Filled-circle synthetic test, 1000 km aperture (stations_circle_1000km_N16_fill.txt, 44 stations, 14 events, 28
records), current defaults (spec_K 4, spec_Kt 3, update_ref 'lin', random prior starts) with pt_nit 3000.
Result: xi 1.041/1.072/1.060/0.995 (true 1.077/1.064/1.001/1.001); truth in the xi 5-95 % band for 75 % of the
profile; R-hat xi median 2.07 (not converged). Same pattern as the line tests with the step-1 isotropic reference
(2c/2d/2f): lid anisotropy too weak, excess anisotropy at 160-250 km. Compared with x1k12 (1000 km, 12 events;
1.000/1.082/1.056/1.013, coverage 56 %): better lid and coverage, the same deep excess.
Fundamental-mode nuisance: the recovered tau of 15 of 28 records is shifted by 1-3 whole cycles, and B differs from
its fit at the true structure by rms |B| 3.2 / phase 60 deg. The logL of the true structure (68556) is not
comparable with the chains (565-566k): its own nuisance fit is cycle-skipped in 13 branch-0 records.
Ran with SIGSTOP pauses for line tests 1b-2g (wall time not representative).
