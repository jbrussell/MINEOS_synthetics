Project line2g
  run tag _line2g, sampler pt, archived 2026-11-28 02:51 from /tmp/claude-0/-home-user-MINEOS-synthetics/578903d7-e992-5265-9feb-4720fd8705a2/scratchpad/work11/RL_array_multimode_waveform_inv

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
  update_ref           = 'burnin'
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

Results: figures/summary.txt (statistics), figures/posterior_summary_line2g.mat
  (pooled posterior and its statistics), figures/*.png, step2/chain_line2g_XX.mat (samples).

Folders: parameters/ (parameters_RL.m of the run, opt_used.mat), figures/, step1/, step2/, data/ (waveform data and model cards)

Rerun (from RL_array_multimode_waveform_inv/):
  >> opt = load_project('line2g');          % the exact options of this run
  >> opt = load_project('line2g', true);    % the same, with the archived waveform data (data/)
  then, as for a new run: D1_step1_Bayesian_disp (step 1), D2_step2_prep_mineos_idagrn,
  D3_step2_waveform_inv (D4_step2_plot_results is called at the end). For synthetic data
  made anew, A0_synth_make_sac_data, A1_sac2eventmat, A2_cleanup_events, A3_run_autowinpick first.
  To change a setting, edit opt (e.g. opt.spec_K = 12) and give the run a new opt.project_tag.

--- Repository copy: figures, summary, settings and parameters only (chains, prep and data files stay in the
working directory; *.mat files are not tracked). ---
Single-event line test 2g: as test 2b (results/line2b: step-1 isotropic reference, tau held at the truth, B free,
arrival angles held at ~0, random prior starts) but update_ref 'burnin': the whole waveform reference (MINEOS +
idagrn6) recomputed at the end of burn-in from the mean of 120 late burn-in samples of all ensembles, the true
path delays re-applied relative to it (refix_tau.m).
Reference update: logL change of the end-of-burn-in states [926 -92 785 638] (accepted); logL of the true
structure 62953.9 -> 64750.4 (+1796): the isotropic step-1 reference was costing the truth ~1800 in logL.
Result: xi 1.068/1.057/1.028/1.005 (true 1.077/1.064/1.001/1.001; 2b 1.070/1.047/1.016/0.993, 2f
1.029/1.061/1.045/1.008): the lid anisotropy is recovered and the deep excess roughly halved vs 2f; truth in the xi
5-95 % band for 74 % of the profile (2f 30 %); R-hat xi median 1.91 (2b 3.2, 2f 1.34). Median logL per ensemble
62.3-64.4k vs truth 64.75k: the truth is now above every ensemble (the sampler has not reached the maximum;
not a biased likelihood). Love fundamental: two ensemble groups, one on the truth, the other ~0.9 % fast at 15 s
(the 2f bias). The tau "rms error" lines of the summary compare the held tau (0, relative to the shifted
reference) with the delay of the truth relative to the new reference: a reporting artefact of fix_tau + burnin.
