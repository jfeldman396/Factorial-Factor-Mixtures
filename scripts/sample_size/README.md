# Sample-Size Simulation Scripts

This folder contains the runnable scripts for the fixed-DGP IFEval-like
simulation study.

The current paper-facing run label is:

```text
fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_subspace2e3_seedfix_full_grid
```

## Canonical Workflow

1. Generate the true loading matrices and heatmaps:

   ```bash
   Rscript scripts/sample_size/plot_fixed_ifeval_lambda_heatmaps.R
   ```

2. Launch or resume the simulation:

   ```bash
   zsh scripts/sample_size/run_full_grid_lambda5_subspace2e3_simulation.sh
   ```

   The default grid is `n = 100, 200, 400`, Product MAP `p = 500, 1000,
   1500, 2000`, Gibbs `p = 500, 1000`, `H = 5, 10`, all-2/all-3 marginal
   component counts, 25 replications, IFEval-like block sizes with at least
   30 items in the smallest block, `Uniform(1, 2)` nonzero loadings,
   cross-loading probability `0.05`, mixture separation `2`, shared
   Product MAP/Viroli-Laplace penalty `lambda = 5`, and EM-SVD pretraining
   stopped when the left singular subspace change is below `2e-3`.

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

   To run the two light robustness arms, use:

   ```bash
   zsh scripts/sample_size/run_lambda5_sensitivity_asym_overlap.sh
   ```

   This keeps the main loading design and estimator settings, then checks
   moderate asymmetric mixture weights and moderate extra component overlap.

3. Plot interim or final results from the chunk checkpoints:

   ```bash
   Rscript scripts/sample_size/plot_fixed_ifeval_lambda_progress.R
   ```

   To plot only the Product MAP estimator:

   ```bash
   METHOD_FILTER=independent_marginal_mixture \
   OUTPUT_TAG=product_map_only \
   PLOT_DIR=results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_subspace2e3_seedfix_full_grid/product_map_only \
   Rscript scripts/sample_size/plot_fixed_ifeval_lambda_progress.R
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
- `run_lambda5_sensitivity_asym_overlap.sh`: runs the asymmetric-probability
  and moderate-overlap robustness checks using the same core grid as the main
  simulation, excluding only the dominated diffuse-Gaussian Gibbs baseline.
