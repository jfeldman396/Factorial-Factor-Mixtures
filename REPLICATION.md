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

The current final simulation asks when the first-stage binary signal/subspace
estimate is accurate enough for rotation and MAP refinement.  The design and
algorithms are summarized in:

```text
writeup/final_simulation_design/final_simulation_design_algorithms.pdf
```

Run the launcher from the repository root:

```sh
Rscript scripts/sample_size/run_final_product_viroli_simulation.R
```

By default this runs:

- `n in {100, 200}`;
- `p in {500, 1000, 1500, 2000}` for product MAP;
- `p in {500, 1000}` for the Viroli baselines;
- `H in {5, 10, 15, 20}`;
- `G in {2, 3}`, expanded to equal component counts across factors;
- separation `1`;
- loading magnitudes `Uniform(1.25, 1.75)` and `Uniform(2.50, 3.00)`;
- cross-loading probabilities `0.075` and `0.20`;
- randomly signed cross-loadings;
- balanced and IFEval-like unbalanced item blocks;
- 25 Monte Carlo repetitions per setting;
- item intercepts using the IFEval-like intercept design;
- loading-based sign/permutation alignment for all recovery metrics;
- Product MAP with EM-SVD likelihood pretraining, sparse rotation, and MAP
  refinement;
- Viroli Gibbs with a Laplace loading prior and penalty 10;
- Viroli Gibbs with a diffuse Gaussian loading prior.

The main output folder is:

```text
results/full/signal_support_grid
```

The run is resumable at the task-chunk level.  Each chunk writes its own
`comparison_results_checkpoint.csv` under
`results/full/signal_support_grid/chunks`.  The launcher combines completed
chunks into:

```text
results/full/signal_support_grid/comparison_results.csv
results/full/signal_support_grid/comparison_summary.csv
```

These full outputs are ignored by git because they can become large.  The
checkpoint records the full DGP setting, penalties, mixture priors,
convergence diagnostics, elapsed seconds, Gibbs ESS summaries, item-prevalence
diagnostics, loading support diagnostics, and product-MAP first-stage
signal/subspace diagnostics.

### Final Simulation Smoke Tests

For a fast Product MAP code-path check:

```sh
N_VALUES=20 \
P_VALUES_PRODUCT=40 \
P_VALUES_GIBBS=99999 \
H_VALUES=2 \
G_VALUES=2 \
LOADING_STRENGTHS=weak \
CROSS_LOADING_PROBS=0.075 \
BLOCK_SIZE_MODES=balanced \
REP_VALUES=1 \
TASK_WORKERS_PRODUCT=1 \
PRODUCT_INTERNAL_WORKERS=2 \
RUN_LABEL=signal_support_grid_smoke \
Rscript scripts/sample_size/run_final_product_viroli_simulation.R
```

To smoke-test all three method paths, use a tiny Gibbs run:

```sh
N_VALUES=12 \
P_VALUES_PRODUCT=20 \
P_VALUES_GIBBS=20 \
H_VALUES=2 \
G_VALUES=2 \
LOADING_STRENGTHS=weak \
CROSS_LOADING_PROBS=0.075 \
BLOCK_SIZE_MODES=balanced \
REP_VALUES=1 \
TASK_WORKERS_PRODUCT=1 \
PRODUCT_INTERNAL_WORKERS=2 \
TASK_WORKERS_GIBBS=2 \
VIROLI_ITER=6 \
VIROLI_BURN=3 \
VIROLI_COMPUTE_PARAMETER_ESS=FALSE \
RUN_LABEL=signal_support_grid_smoke_gibbs \
Rscript scripts/sample_size/run_final_product_viroli_simulation.R
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

For Viroli Gibbs, the default is six independent task chunks at a time and one
internal Gibbs worker per chunk:

```text
TASK_WORKERS_GIBBS=6
GIBBS_INTERNAL_WORKERS=1
```

This parallelizes Gibbs across independent replications/configurations, which
is usually more efficient than trying to parallelize every MCMC transition
inside a single chain.  Increase `TASK_WORKERS_GIBBS` if memory permits.

### Regenerate Simulation Figures

Representative DGP heatmaps for the final loading designs:

```sh
Rscript scripts/sample_size/plot_dgp_loading_heatmaps.R
```

This writes PNG heatmaps to:

```text
results/selected_plots/sample_size/signal_support_grid/dgp_heatmaps
```

and matching loading matrices to:

```text
results/selected_tables/sample_size/signal_support_grid_dgp
```

Progress plots from completed chunks:

```sh
Rscript scripts/sample_size/plot_signal_support_simulation_progress.R
```

This reads `results/full/signal_support_grid`, writes a completed-results
snapshot to:

```text
results/selected_tables/sample_size/signal_support_grid_completed_results.csv
```

and writes line/boxplot summaries under:

```text
results/selected_plots/sample_size/signal_support_grid
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
- `ess_min`, `ess_median`, `ess_mean`: Gibbs effective sample size summaries
  when ESS calculation is enabled.

The DGP columns `dgp_min_total_nonzero_loadings_by_factor`,
`dgp_mean_cross_loadings_per_item`, `dgp_min_loading_l2_by_factor`, and related
fields connect recovery to effective signal strength and loading support.

## Committed Selected Results

Selected outputs are committed for immediate inspection:

- final simulation DGP heatmaps:
  `results/selected_plots/sample_size/signal_support_grid/dgp_heatmaps`;
- final simulation representative Lambda matrices:
  `results/selected_tables/sample_size/signal_support_grid_dgp`;
- IFEval plots: `results/selected_plots/ifeval`;
- IFEval tables: `results/selected_tables/ifeval`;
- IFEval writeup: `writeup/ifeval_componentwise_G3313.pdf`;
- simulation design writeup:
  `writeup/final_simulation_design/final_simulation_design_algorithms.pdf`.

The latest static audit notes are in `CODE_AUDIT.md`. They record which R files
were parsed, which input files were checked, which PDF was rendered, and which
caveats remain before claiming a full fresh end-to-end rerun.

These selected outputs are snapshots. Full reruns should write to `results/full/`.

## Git Notes

Generated `.rds` files, logs, and full output directories are ignored by `.gitignore`. Commit source code, replication instructions, selected summary tables, selected plots, and final writeups.
