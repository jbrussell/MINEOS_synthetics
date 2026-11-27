Project line1c
  run tag _line1c, sampler pt, archived 2026-11-27 19:59 from /tmp/claude-0/-home-user-MINEOS-synthetics/578903d7-e992-5265-9feb-4720fd8705a2/scratchpad/work11/RL_array_multimode_waveform_inv

Key settings (all of them in settings.txt and parameters/opt_used.mat):
  trueid               = 'pa5_5km'
  aswms_dir            = 'ASWMS_line/'
  nbranch_data         = 4
  nm                   = 4
  sampler              = 'pt'
  nens                 = 4
  K                    = 10
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
  update_ref           = 'none'
  lin_updates          = 3
  smooth_sig           = [0.0002 0.0002]
  smooth_hier          = false
  vsv_ref_sig          = 'step1'
  xi_ref_hier          = true
  synth_arrival_shift  = []
  noise                = 0.005
  s1_cstd              = 0.0025
  prepfile             = 'prep_line1.mat'

Synthetic data (A0, ASWMS_line/synth_truth.mat): true model 'pa5_5km', data branches 4, noise 0.005 x peak, seed 3, station noise spread 0

Results: figures/summary.txt (statistics), figures/posterior_summary_line1c.mat
  (pooled posterior and its statistics), figures/*.png, step2/chain_line1c_XX.mat (samples).

Folders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/, data/ (waveform data and model cards)

Rerun (from RL_array_multimode_waveform_inv/):
  >> opt = load_project('line1c');          % the exact options of this run
  >> opt = load_project('line1c', true);    % the same, with the archived waveform data (data/)
  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,
  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data
  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.
  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.

--- Repository copy: figures, summary, settings and parameters only (chains, prep and data files stay in the
working directory; *.mat files are not tracked). ---
Single-event line test 1c (oracle): as test 1b (results/line1b: true reference pa5_5km, fix_nuisance, jump_every 0,
update_ref 'none', PT 4 ensembles x 10 temperatures up to 20, pt_nit 1500, prep_line1.mat), but every PT state starts
at the reference = true model (pt_start 'ref': least-squares B-spline fit with 8 interior points; the smooth profile
cannot reproduce the true steps at 60 and 160 km: rms 0.031 km/s in Vsv, 0.005 in xi).
Result: the true mode is stable. All ensembles stay in it (median logL 65.3-65.7k, 1.2-1.6k below the true stepped
model, 66.9k), R-hat xi median 1.56 (1b: 2.6), xi 1.076/1.058/1.011/1.005 (true 1.077/1.064/1.001/1.001), phase-velocity
tracks within +-0.3 % of the truth (1b: up to 1.2 %). Test 1b's failure is finding the mode from random prior starts,
not the posterior itself.
