Project line2b
  run tag _line2b, sampler pt, archived 2026-11-27 19:32 from /tmp/claude-0/-home-user-MINEOS-synthetics/578903d7-e992-5265-9feb-4720fd8705a2/scratchpad/work11/RL_array_multimode_waveform_inv

Key settings (all of them in settings.txt and parameters/opt_used.mat):
  trueid               = 'pa5_5km'
  aswms_dir            = 'ASWMS_line/'
  nbranch_data         = 4
  nm                   = 4
  sampler              = 'pt'
  nens                 = 4
  K                    = 6
  Tmax                 = 20
  pt_nit               = 1500
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
  noise                = 0.005
  s1_cstd              = 0.0025
  prepfile             = 'prep_line2.mat'

Synthetic data (A0, ASWMS_line/synth_truth.mat): true model 'pa5_5km', data branches 4, noise 0.005 x peak, seed 3, station noise spread 0

Results: figures/summary.txt (statistics), figures/posterior_summary_line2b.mat
  (pooled posterior and its statistics), figures/*.png, step2/chain_line2b_XX.mat (samples).

Folders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/, data/ (waveform data and model cards)

Rerun (from RL_array_multimode_waveform_inv/):
  >> opt = load_project('line2b');          % the exact options of this run
  >> opt = load_project('line2b', true);    % the same, with the archived waveform data (data/)
  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,
  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data
  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.
  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.

--- Repository copy: figures, summary, settings and parameters only (chains, prep and data files stay in the
working directory; *.mat files are not tracked). ---
Single-event line test 2b: as test 2 (results/line2: event 201210280304A, 41 stations every 50 km on a 2000 km line
along its great circle through the NoMelt centre; noise 0.005 of the peak; step-1 (isotropic) reference, prep_line2.mat;
tau held at the true path delays, fix_tau; B free; current defaults otherwise: spec_K 4, spec_Kt 3, PT 4 x 6, pt_nit
1500, update_ref 'lin'), but the arrival angles held at ~0 (dth_prior = dth1_prior = 0.01 deg).
Result: the radial anisotropy is recovered (xi 1.070/1.047/1.016/0.993 in 7-60/60-160/160-250/250-350 km, true
1.077/1.064/1.001/1.001; test 2 with free dth gave ~1.000 everywhere): on a line along the great circle a free
arrival angle trades off with the Love phase velocity (cos dth / c) and absorbed the anisotropy in test 2.
R-hat xi median 3.2 (ensembles not mixed), median logL 61.9-62.9k vs 62.95k at the truth.
