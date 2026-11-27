Project line2d
  run tag _line2d, sampler pt, archived 2026-11-27 22:51 from /tmp/claude-0/-home-user-MINEOS-synthetics/578903d7-e992-5265-9feb-4720fd8705a2/scratchpad/work11/RL_array_multimode_waveform_inv

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

Results: figures/summary.txt (statistics), figures/posterior_summary_line2d.mat
  (pooled posterior and its statistics), figures/*.png, step2/chain_line2d_XX.mat (samples).

Folders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/, data/ (waveform data and model cards)

Rerun (from RL_array_multimode_waveform_inv/):
  >> opt = load_project('line2d');          % the exact options of this run
  >> opt = load_project('line2d', true);    % the same, with the archived waveform data (data/)
  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,
  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data
  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.
  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.

--- Repository copy: figures, summary, settings and parameters only (chains, prep and data files stay in the
working directory; *.mat files are not tracked). ---
Single-event line test 2d: as test 2c (results/line2c: step-1 isotropic reference, tau held at the truth, B free,
arrival angles held at ~0) but pt_start 'refT1_TNjitter': the T = 1 state of every ensemble starts at the step-1 model,
the hotter ones at it with control-point values perturbed by N(0, [0.05 km/s, 0.02]).
Result: no better than 2b (random starts) or 2c: R-hat xi median 3.3, xi 1.018/1.062/1.045/1.003 (true
1.077/1.064/1.001/1.001): the lid anisotropy is pushed below 160 km (a depth trade-off along which the chains from the
isotropic start stay), truth in the xi 5-95 % band for only 34 % of the profile; median logL 62.2-62.7k (truth 62.95k).
