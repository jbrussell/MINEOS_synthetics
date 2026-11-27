Project line1b
  run tag _line1b, sampler pt, archived 2026-11-27 16:18 from /tmp/claude-0/-home-user-MINEOS-synthetics/578903d7-e992-5265-9feb-4720fd8705a2/scratchpad/work11/RL_array_multimode_waveform_inv

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

Results: figures/summary.txt (statistics), figures/posterior_summary_line1b.mat
  (pooled posterior and its statistics), figures/*.png, step2/chain_line1b_XX.mat (samples).

Folders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/, data/ (waveform data and model cards)

Rerun (from RL_array_multimode_waveform_inv/):
  >> opt = load_project('line1b');          % the exact options of this run
  >> opt = load_project('line1b', true);    % the same, with the archived waveform data (data/)
  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,
  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data
  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.
  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.

--- Repository copy: figures, summary, settings and parameters only (chains, prep and data files stay in the
working directory; *.mat files are not tracked). ---
Single-event line test 1b (oracle): as test 1 (results/line1: event 201210280304A, 41 stations every 50 km on a
2000 km line along its great circle through the NoMelt centre; noise 0.005 of the peak; reference = the TRUE model,
opt.ref_card_id = 'pa5_5km'; B, tau, dth held at their true values, fix_nuisance, jump_every 0; prep_line1.mat),
but update_ref 'none' (the linearization stays at the true model; no linearization updates) and a denser PT ladder:
K 10 temperatures (1, 1.39, 1.95, ..., 20) instead of 6. pt_nit 1500, 4 ensembles.
Result: closer to the truth than test 1 (median logL 63.0-64.7k vs 61.7-64.0k; truth 66.9k) and better swap rates,
but still not converged (R-hat xi median 2.6) and still separate dispersion tracks in fig11: the T = 1 states stay in
modes below the true one, so the multimodality is not (only) the linearization drift of test 1.
