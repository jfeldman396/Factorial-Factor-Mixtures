# Sample-Size Simulation Scripts

This folder contains the runnable scripts for the fixed-DGP IFEval-like
simulation study.

The paper-facing study is the three-arm fixed-DGP recovery comparison launched
by `run_three_arm_recovery_study.sh`.

## Canonical Workflow

1. Generate the true loading matrices and heatmaps:

   ```bash
   Rscript scripts/sample_size/plot_fixed_ifeval_lambda_heatmaps.R
   ```

2. Launch or resume the complete recovery simulation:

   ```bash
   zsh scripts/sample_size/run_three_arm_recovery_study.sh
   ```

   The default grid is `n = 100, 200, 400`, Product MAP `p = 500, 1000,
   1500, 2000`, Gibbs `p = 500, 1000`, `H = 5, 10`, all-2/all-3 marginal
   component counts, 25 replications, IFEval-like block sizes with at least
   30 items in the smallest block, `Uniform(1, 2)` nonzero loadings,
   cross-loading probability `0.05`, mixture separation `2`, shared Product
   MAP/Gibbs loading penalties `lambda = 3` at `n = 100` and `lambda = 5` at
   `n = 200, 400`, and EM-SVD pretraining stopped when the left singular
   subspace change is below `2e-3`.

   The launcher runs three DGP arms with one common loading/intercept design:
   separated symmetric mixtures, asymmetric mixture probabilities, and
   moderately overlapping mixtures. Product MAP is run at every `p`; Gibbs
   with the matched Laplace prior is run at `p = 500, 1000`. Set
   `RUN_VIROLI_GAUSSIAN=TRUE` to add the diffuse-Gaussian Gibbs baseline.

   Product MAP is scored after final canonical normalization: each fitted
   factor's marginal mixture has mean zero and variance one, with `alpha` and
   `Lambda` transformed to preserve the probit linear predictor. Viroli Gibbs
   uses the same canonical parameterization by normalizing each posterior draw
   before averaging with the same `canonical_normalize_factor_parameters()`
   helper and the same `CANONICAL_MIN_SCALE` constant. Raw non-scaled
   parameters are diagnostic only.

   Keep `STABLE_SCENARIO_SEEDS=TRUE` for method comparisons; this makes the
   generated dataset depend on the scenario itself rather than the position of
   that scenario inside a launcher chunk.

   Useful controls are:

   ```bash
   # Recreate figures from completed results without fitting.
   RUN_FITS=FALSE RUN_PLOTS=TRUE \
   zsh scripts/sample_size/run_three_arm_recovery_study.sh

   # Run one arm only.
   ARM_FILTER=moderate_overlap \
   zsh scripts/sample_size/run_three_arm_recovery_study.sh
   ```

   Historical partial-run utilities remain available:
   `run_n100_lambda3_study.sh` reproduces the `n = 100` block, and
   `run_product_map_robustness_completion.sh` fills Product MAP cells for
   `n = 200, 400`. They are recovery tools, not the preferred fresh-run entry
   point.

3. The authoritative launcher creates all eight presentation boxplots for
   every arm after fitting: factor scores, probabilities, loadings,
   intercepts, mixture means, mixture variances, mixture weights, and runtime.
   To plot one run directly:

   ```bash
   RUN_LABEL=fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_separated \
   G_COMPONENT_FILTER=2,3 \
   Rscript scripts/sample_size/plot_fixed_ifeval_grouped_boxplot_panels.R
   ```

4. Recreate an example true/estimated Lambda recovery panel, factor-score
   scatter panel, and fitted marginal-mixture panel for one design cell:

   ```bash
   N_VALUE=200 P_VALUE=500 H_TRUE=5 G_TRUE=3 REP_VALUE=1 \
   LASSO_PENALTY=5 SEPARATIONS=2 \
   Rscript scripts/sample_size/plot_example_lambda_recovery.R
   ```

## Script Roles

- `run_fixed_ifeval_lambda_simulation.R`: high-level launcher for the fixed
  IFEval-like design. It sets the simulation grid and worker allocation.
- `run_three_arm_recovery_study.sh`: authoritative, resumable full-study
  launcher. It fixes the three DGP arms, canonicalization, penalties,
  convergence criteria, worker counts, output paths, and final plotting.
- `validate_three_arm_recovery_study.R`: verifies expected row coverage,
  duplicate-free scientific keys, and matched DGP seeds across methods.
- `compare_original_simulation_joint_mfa_gibbs.R`: low-level fitting engine
  called by the launcher for each chunk.
- `plot_fixed_ifeval_lambda_heatmaps.R`: deterministic DGP visualization and
  Lambda table export.
- `plot_fixed_ifeval_lambda_progress.R`: reads checkpoint files from a running
  simulation and creates line/boxplot summaries. It also writes a cell-level
  summary table, a Product MAP versus Gibbs overlap table once Gibbs rows are
  available, a matched Product MAP versus Viroli Laplace scatter plot, and
  fixed-`H,G` recovery panels showing how RMSE changes with `p`.
- `plot_fixed_ifeval_grouped_boxplot_panels.R`: creates the presentation-style
  grouped boxplots by metric, faceted by `p` and `H/G`.
- `plot_example_lambda_recovery.R`: recreates one simulation cell and writes
  side-by-side true/Product MAP/Viroli Laplace heatmaps, aligned factor-score
  plots, and fitted marginal-mixture overlays.
- `run_n100_lambda3_study.sh`: partial-run recovery utility for the `n = 100`
  lambda-3 block of all three mixture settings.
- `run_product_map_robustness_completion.sh`: partial-run recovery utility for
  the missing `n = 200, 400` Product MAP robustness cells using a common
  `p_max = 2000` master loading matrix.
- `run_lambda5_sensitivity_asym_overlap.sh`: superseded five-replication pilot
  for the asymmetric and overlap arms; retained only for audit history.
- `run_rotation_ablation_three_arms.sh`: runs the paired rotation ablation over
  the full `n`, `p`, `H`, and `G` grid for all three mixture settings. Each
  dataset is evaluated with FastICA only, FastICA followed by mixture rotation,
  and identity/SVD initialization followed by mixture rotation, using both the
  estimated and oracle latent-response signals. It also fits an overcomplete
  rank-15 estimated signal and records eigengap-based estimates of `H`.
- `compare_rotation_fastica_diagnostic.R`: fitting engine for the rotation
  ablation. It records pre-refinement and post-refinement recovery, end-to-end
  runtime, estimated-versus-oracle signal error, and the full eigengap profile.
