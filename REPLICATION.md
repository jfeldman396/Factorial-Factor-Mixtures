# Replication Instructions

This guide explains how to replicate the two analysis tracks in this repository:

1. the IFEval empirical analysis;
2. the sample-size simulation comparing the proposed factorial factor-mixture method with Viroli-style Gibbs samplers.

Run these commands in a local macOS Terminal, not in the GitHub web interface.
All commands assume the repository root is the working directory:

```sh
cd "/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
```

## Software

The code is written primarily in R. The main scripts use:

- base R and `parallel`;
- `MASS`;
- `truncnorm`;
- `glmnet`;
- `coda`;
- `ggplot2`;
- `reshape2`;
- `mclust`.

Optional 3D HTML visualizations use Python packages:

- `numpy`;
- `pandas`;
- `matplotlib`;
- `plotly`.

The scripts were developed on macOS. The sample-size simulation uses Unix-style multicore parallelism where available.

## Repository Layout

Important paths:

- `R/`: shared copies of the proposed binary probit pretraining and refinement algorithms.
- `scripts/ifeval/`: IFEval fitting, tuning, comparison, and visualization scripts.
- `scripts/sample_size/`: simulation, Gibbs comparator, and checkpoint plotting scripts.
- `data/ifeval/`: cleaned IFEval binary matrix and item/model metadata.
- `configs/sample_size_intercepts_centered.env`: the main simulation grid/settings.
- `results/selected_tables/`: selected CSV outputs committed for review.
- `results/selected_plots/`: selected figures committed for review.
- `writeup/`: IFEval writeups and rendered PDFs.
- `CODE_AUDIT.md`: current static audit notes, parse checks, and caveats.

Full regenerated outputs should be written under `results/full/`, which is ignored by git.

## IFEval Analysis

The cleaned input files are:

- `data/ifeval/openeval_ifeval_only_binary_matrix.csv`;
- `data/ifeval/openeval_item_metadata.csv`;
- `data/ifeval/openeval_model_metadata.csv`.

The current IFEval workflow fits only the proposed independent-mixture probit
factor model.  Rank `H`, component counts `G_h`, and the entrywise sparse
loading penalty are tuned by held-out predictive log likelihood on randomly
removed model-by-item cells.  The current component-wise grid allows
`G_h in {1,2,3}` with at most one Gaussian coordinate (`G_h = 1`).

### Run The Current IFEval Pipeline

The easiest way to reproduce the current analysis is:

```sh
bash scripts/ifeval/run_full_analysis.sh
```

This runs held-out-likelihood CV, selected mixture refit, Lambda heatmaps,
marginal mixture plots, loading summaries, and factor-score visualizations.

### Tune H, Component-Wise G, And Sparse Loading Penalty

```sh
OUT_DIR=results/full/ifeval_columnwise_G_cv_atmost1_gaussian \
H_GRID=1:5 \
G_MODE=column_grid \
G_GRID=1,2,3 \
G_COMPONENT_VALUES=1,2,3 \
MAX_GAUSSIAN_COORDS=1 \
LAMBDA_L1_GRID=0,1,2,4,8,12 \
K_FOLDS=3 \
WORKERS=6 \
PRETRAIN_AUG_ITER=20 \
REFINE_ITER=20 \
MIXTURE_MAX_ITER=200 \
FIT_SELECTED_AFTER_CV=TRUE \
SAVE_FITS=FALSE \
RESUME_EXISTING=TRUE \
REFRESH_PLOTS=TRUE \
Rscript scripts/ifeval/cv_ifeval_rank_lambda_models.R
```

Primary outputs:

- `ifeval_rank_lambda_cv_fold_scores.csv`;
- `ifeval_rank_lambda_cv_histories.csv`;
- `ifeval_rank_lambda_cv_summary.csv`;
- `ifeval_rank_lambda_selected_by_heldout_ll.csv`;
- `rank_lambda_top_candidates_by_heldout_ll.csv`;
- `rank_lambda_top_candidates_by_heldout_ll.png`;
- selected full-data mixture refit folder when `FIT_SELECTED_AFTER_CV=TRUE`.

### Refit The Interpretable Component-Wise Model

The current component-wise interpretation uses:

- `H = 4`;
- `G = (3,3,1,3)`, so the third coordinate is Gaussian;
- sparse-loading MAP refinement with `lambda_l1_penalty = 4`;
- `20` maximum pretraining iterations and `20` maximum refinement iterations.

Run:

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
ITEM_METADATA_PATH=data/ifeval/openeval_item_metadata.csv \
OUT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
H_FIXED=4 \
G_FIXED=3,3,1,3 \
WORKERS=8 \
PRETRAIN_AUG_ITER=20 \
REFINE_ITER=20 \
MIXTURE_MAX_ITER=200 \
REQUIRE_MIXTURE_CONVERGENCE=TRUE \
REFINEMENT_LAMBDA_L1_PENALTY=4 \
Rscript scripts/ifeval/fit_interpret_ifeval_mixture.R
```

Primary outputs:

- `openeval_item_intercepts_loadings_metadata.csv`;
- `openeval_model_factor_scores_profiles.csv`;
- `openeval_factor_interpretation_summary.csv`;
- `openeval_top_loading_item_examples.csv`;
- `openeval_factor_mixture_groups.csv`;
- loading, marginal-mixture, factor-score, and profile plots.

### Summarize Loadings And Cross-Loadings

```sh
LOADINGS_PATH=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4/openeval_item_intercepts_loadings_metadata.csv \
OUT_DIR=results/full/ifeval/loadings_crossloadings \
Rscript scripts/ifeval/summarize_ifeval_loadings_crossloadings.R
```

Primary outputs:

- `ifeval_sparse_learned_lambda_all_items_ordered.csv`;
- `ifeval_sparse_lambda_matrix_threshold_0p5.csv`;
- `ifeval_G3_factor_only_exact_items.csv`;
- `ifeval_G3_cross_loading_items.csv`;
- heatmaps for sparse and cross-loading items.

### Regenerate Selected-Fit Visualizations

```sh
FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/selected_visualizations \
Rscript scripts/ifeval/plot_openeval_mixture_lambda_heatmaps.R

FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/selected_visualizations \
Rscript scripts/ifeval/plot_factor_marginal_mixtures.R

FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/ifeval_3d_factor_visualizations \
Rscript scripts/ifeval/plot_ifeval_learned_factors_baseR.R
```

Optional Python/Plotly plots:

```sh
FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/ifeval_3d_factor_visualizations \
python3 scripts/ifeval/plot_ifeval_learned_factors.py
```

The rendered component-wise writeup is:

```text
writeup/ifeval_componentwise_G3313.pdf
```

It summarizes the selected fit, compact cross-loading examples, and LLM
ability profiles by factor.

## Sample-Size Simulation

The current paper-facing simulation is the fixed-DGP IFEval-like Lambda
simulation.  It asks whether Product MAP can recover latent factors and model
parameters competitively with Viroli-style probit Gibbs samplers when the item
loading structure resembles the IFEval analysis.

```text
docs/fixed_ifeval_lambda_simulation_design.md
docs/fixed_ifeval_lambda_simulation_design.pdf
```

Run or resume the full simulation from the repository root:

```sh
Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

By default this runs:

- `n in {100, 200, 400}`;
- `p in {500, 1000, 1500, 2000}` for Product MAP;
- `p in {500, 1000}` for the Viroli baselines;
- `H in {5, 10}`;
- `G_h = 2` for every factor, or `G_h = 3` for every factor;
- separation `2`;
- IFEval-like unbalanced item blocks, with at least 30 primary items in the
  smallest block;
- nonzero loading magnitudes sampled from `Uniform(2, 3)`;
- cross-loading probability `0.05`;
- randomly signed cross-loadings and block-level primary-loading signs;
- 25 Monte Carlo replications per setting;
- item intercepts using the IFEval-like intercept design;
- loading-based sign/permutation alignment for all recovery metrics;
- Product MAP with EM-SVD likelihood pretraining, Riemannian rotation, and MAP
  refinement;
- Viroli Gibbs with a Laplace loading prior using the same n-dependent loading
  penalty schedule as Product MAP;
- Viroli Gibbs with a diffuse Gaussian loading prior.

The main output folder is:

```text
results/full/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10
```

The run is resumable at the task-chunk level.  Each chunk writes its own
`comparison_results_checkpoint.csv` under
`results/full/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/chunks`.
The launcher combines completed chunks into:

```text
results/full/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/comparison_results.csv
```

The completed main run has 2400 result rows: 1200 Product MAP rows, 600 Viroli
Laplace Gibbs rows, and 600 Viroli Gaussian Gibbs rows.  Full outputs under
`results/full/` are ignored by git because they contain logs and chunk-level
artifacts.  Selected CSV snapshots and plots are committed under
`results/selected_tables/sample_size/` and
`results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/`.

The default penalty schedule is:

```text
LAMBDA_L1_PENALTY_BY_N=100=5,200=5,400=8
```

For Product MAP this sets `PRETRAIN_LOADING_PENALTY`,
`ROTATION_LOADING_L1_PENALTY`, and `LAMBDA_L1_PENALTY`.  For Viroli Laplace it
sets `VIROLI_LAMBDA_L1_PENALTY`.  Viroli Gaussian sets the Laplace penalty to
zero.

### Final Simulation Smoke Tests

For a fast Product MAP code-path check:

```sh
N_VALUES=20 \
P_VALUES_PRODUCT=40 \
P_VALUES_GIBBS=40 \
H_VALUES=2 \
G_CONFIG_TYPES=all2 \
REP_VALUES=1 \
TASK_WORKERS_PRODUCT=1 \
PRODUCT_INTERNAL_WORKERS=2 \
VIROLI_ITER=6 \
VIROLI_BURN=3 \
VIROLI_COMPUTE_PARAMETER_ESS=FALSE \
RUN_LABEL=fixed_ifeval_lambda_smoke \
Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

To smoke-test only Product MAP, disable Gibbs by giving an empty Gibbs grid:

```sh
N_VALUES=20 \
P_VALUES_PRODUCT=40 \
P_VALUES_GIBBS= \
H_VALUES=2 \
G_CONFIG_TYPES=all2 \
REP_VALUES=1 \
TASK_WORKERS_PRODUCT=1 \
PRODUCT_INTERNAL_WORKERS=2 \
RUN_LABEL=fixed_ifeval_lambda_smoke_product_only \
Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

These smoke tests only verify code paths.  They are not scientific simulations.

### Parallelization

The launcher uses two levels of parallelism.

For Product MAP, the default is one task chunk at a time with 18 internal
workers:

```text
TASK_WORKERS_PRODUCT=1
PRODUCT_INTERNAL_WORKERS=18
```

Within Product MAP, internal workers are used for independent marginal mixture
fits, itemwise loading regressions, and subject-wise factor-score updates.
The outer EM-SVD, rotation, and refinement sweeps remain sequential.

For Viroli Gibbs, the current launcher runs one task chunk at a time and uses
four internal workers inside each Gibbs fit:

```text
TASK_WORKERS_GIBBS=1
GIBBS_INTERNAL_WORKERS_SERIAL=4
GIBBS_INTERNAL_WORKERS_PARALLEL=4
```

This keeps replication scheduling simple while still parallelizing the
computationally heavy Gibbs conditionals where the implementation supports it.

### Regenerate Simulation Figures

Representative DGP heatmaps for the fixed IFEval-like loading design:

```sh
Rscript scripts/sample_size/plot_fixed_ifeval_lambda_heatmaps.R
```

This writes PNG heatmaps to:

```text
results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/true_lambda_heatmaps
```

and matching loading matrices to:

```text
results/selected_tables/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/true_lambda
```

Progress plots from completed chunks:

```sh
RUN_LABEL=fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10 \
Rscript scripts/sample_size/plot_fixed_ifeval_lambda_progress.R
```

This reads the corresponding `results/full/` directory, writes a completed-results
snapshot to:

```text
results/selected_tables/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10_completed_results.csv
```

and writes line/boxplot summaries under:

```text
results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10
```

To regenerate an example loading-recovery panel for one cell:

```sh
N_VALUE=200 P_VALUE=500 H_TRUE=5 G_TRUE=3 REP_VALUE=1 \
LASSO_PENALTY=5 SEPARATIONS=2 \
Rscript scripts/sample_size/plot_example_lambda_recovery.R
```

This writes true/Product MAP/Viroli Laplace loading heatmaps, factor-score
scatter panels, fitted factor-marginal overlays, and aligned matrices under:

```text
results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/lambda_recovery_examples
```

### Interpret Simulation Metrics

The main recovery metrics are computed after aligning estimated factors and
loadings to the data-generating factors and loadings by best signed
permutation.

Important fields:

- `factor_score_rmse`: RMSE between aligned estimated and true factor scores.
- `lambda_rmse`: RMSE between aligned estimated and true loading entries.
- `alpha_rmse`: RMSE between estimated and true item intercepts.
- `marginal_mu_rmse`: RMSE between estimated and true marginal mixture means.
- `marginal_var_rmse`: RMSE between estimated and true marginal mixture variances.
- `marginal_weight_rmse`: RMSE between estimated and true marginal mixture weights.
- `stage1_signal_relative_frobenius_error`: relative error in the first-stage
  centered low-rank probit signal for Product MAP.
- `stage1_sinTheta_op`: operator-norm subspace angle error for the first-stage
  signal estimate.
- `seconds`: wall-clock runtime for the method in that repetition.
- `gibbs_min_parameter_ess`, `gibbs_median_parameter_ess`,
  `gibbs_mean_parameter_ess`: Gibbs effective sample size summaries when ESS
  calculation is enabled.

The DGP columns `dgp_min_total_nonzero_loadings_by_factor`,
`dgp_mean_cross_loadings_per_item`, `dgp_min_loading_l2_by_factor`, and related
fields connect recovery to effective signal strength and loading support.

## Committed Selected Results

Selected outputs are committed for immediate inspection:

- fixed IFEval-like simulation DGP heatmaps and recovery plots:
  `results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10`;
- fixed IFEval-like simulation summary CSV snapshots:
  `results/selected_tables/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10_*.csv`;
- IFEval plots: `results/selected_plots/ifeval`;
- IFEval tables: `results/selected_tables/ifeval`;
- IFEval writeup: `writeup/ifeval_componentwise_G3313.pdf`;
- simulation design writeup:
  `docs/fixed_ifeval_lambda_simulation_design.pdf`.

The latest static audit notes are in `CODE_AUDIT.md`. They record which R files
were parsed, which input files were checked, which PDF was rendered, and which
caveats remain before claiming a full fresh end-to-end rerun.

These selected outputs are snapshots. Full reruns should write to `results/full/`.

## Git Notes

Generated `.rds` files, logs, and full output directories are ignored by `.gitignore`. Commit source code, replication instructions, selected summary tables, selected plots, and final writeups.
