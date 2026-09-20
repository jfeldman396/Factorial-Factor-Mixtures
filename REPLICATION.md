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
```

Run or resume the complete three-arm recovery study from the repository root:

```sh
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

This authoritative launcher fixes all scientific and computational settings
needed to reproduce the paper-facing recovery figures. It runs:

- `n in {100, 200, 400}`;
- `p in {500, 1000, 1500, 2000}` for Product MAP;
- `p in {500, 1000}` for the Viroli baselines;
- `H in {5, 10}`;
- `G_h = 2` for every factor, or `G_h = 3` for every factor;
- separation `2`;
- IFEval-like unbalanced item blocks, with at least 30 primary items in the
  smallest block;
- nonzero primary and cross-loading magnitudes sampled from `Uniform(1, 2)`;
- cross-loading probability `0.05`;
- randomly signed cross-loadings and block-level primary-loading signs;
- 25 Monte Carlo replications per setting;
- item intercepts using the IFEval-like intercept design;
- loading-based sign/permutation alignment for all recovery metrics;
- Product MAP with EM-SVD likelihood pretraining, Riemannian rotation, and MAP
  refinement;
- Viroli Gibbs with a Laplace loading prior using the same n-dependent loading
  penalty schedule as Product MAP;
- Viroli Gibbs with a diffuse Gaussian loading prior when
  `RUN_VIROLI_GAUSSIAN=TRUE` (off by default in the three-arm comparison).

The shared loading-penalty schedule is:

```text
LAMBDA_L1_PENALTY_BY_N=100=3,200=5,400=5
```

For Product MAP, the per-`n` value is used in pretraining, mixture rotation,
and MAP refinement. For Gibbs with a Laplace loading prior, the same value is
used in the Gaussian scale-mixture representation of the Laplace prior.

All method comparisons use the same canonical parameterization. Each factor
has marginal mean zero and variance one. Product MAP is canonicalized after
fitting, with the corresponding transformations of `alpha`, `Lambda`, and the
mixture moments preserving the probit linear predictor. Every retained Gibbs
draw is canonicalized by the same shared helper before posterior averaging.

### Mixture Arms

The multiplier below acts on the base means `(-2, 2)` for `G=2` and
`(-2, 0, 2)` for `G=3`, before canonicalization.

| Arm | G | Weights | Mean multiplier | Standard deviations |
|---|---:|---|---:|---|
| separated | 2 | `(0.50, 0.50)` | `1.35` | `(0.45, 0.45)` |
| separated | 3 | `(0.30, 0.40, 0.30)` | `1.35` | `(0.45, 0.65, 0.45)` |
| asymmetric | 2 | `(0.65, 0.35)` | `1.35` | `(0.45, 0.45)` |
| asymmetric | 3 | `(0.20, 0.50, 0.30)` | `1.35` | `(0.45, 0.65, 0.45)` |
| overlap | 2 | `(0.50, 0.50)` | `1.10` | `(0.60, 0.60)` |
| overlap | 3 | `(0.30, 0.40, 0.30)` | `1.10` | `(0.60, 0.75, 0.60)` |

### Output And Row Counts

The three default output folders are:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_separated/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_asymmetric_pi/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_moderate_overlap/
```

With the default two-method comparison, each arm has `1800` unique result
rows: `1200` Product MAP rows and `600` Laplace-Gibbs rows. Enabling the
Gaussian Gibbs baseline adds `600` rows per arm. Each chunk writes a checkpoint
below `<run>/chunks/`; rerunning the launcher skips completed chunks and
rebuilds `<run>/comparison_results.csv`.

Each output root also contains `replication_manifest.txt`, which records the
actual grid, DGP, methods, penalties, and worker allocation used for that run.
The launcher finishes by running `validate_three_arm_recovery_study.R`, which
checks expected row counts, duplicate-free cell/replication keys, and matched
DGP seeds between Product MAP and Gibbs.
Full outputs under `results/full/` are ignored by git; selected figures and
summary tables are written below:

```text
results/selected_plots/sample_size/three_mixture_settings/
results/selected_tables/sample_size/three_mixture_settings/
```

### Replot Without Refitting

```sh
RUN_FITS=FALSE RUN_PLOTS=TRUE \
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

The plotting phase creates one four-by-four faceted figure for each of eight
outcomes: factor-score RMSE, probability RMSE, loading RMSE, intercept RMSE,
mixture-mean RMSE, mixture-variance RMSE, mixture-weight RMSE, and runtime.
Columns are `p`, rows are the four `H/G` combinations, and boxes are grouped by
sample size and method.

### Targeted And Smoke Runs

Use `ARM_FILTER` for one or more comma-separated arms:

```sh
ARM_FILTER=asymmetric_pi,moderate_overlap \
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

For a short end-to-end smoke test, use a new run-label prefix so the diagnostic
does not share chunks with the production study:

```sh
RUN_LABEL_PREFIX=fixed_ifeval_three_arm_smoke \
ARM_FILTER=separated \
N_VALUES=100 \
P_VALUES_PRODUCT=500 \
P_VALUES_GIBBS=500 \
H_VALUES=5 \
G_CONFIG_TYPES=all2 \
REP_VALUES=1,2 \
PRODUCT_INTERNAL_WORKERS=2 \
TASK_WORKERS_GIBBS=1 \
GIBBS_INTERNAL_WORKERS_PARALLEL=2 \
VIROLI_ITER=20 \
VIROLI_BURN=10 \
VIROLI_COMPUTE_PARAMETER_ESS=FALSE \
RUN_PLOTS=FALSE \
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

These smoke settings only verify code paths and must not be used for scientific
comparisons.

### Partial-Run Recovery Utilities

The following launchers reproduce historical blocks or finish interrupted
runs, but are not the preferred entry point for a fresh replication:

- `run_n100_lambda3_study.sh`: all three arms at `n = 100`, with 25 reps.
- `run_product_map_robustness_completion.sh`: Product MAP at `n = 200, 400`
  for the asymmetric and overlap arms, using one `p_max = 2000` master loading
  matrix and 25 reps.
- `run_lambda5_sensitivity_asym_overlap.sh`: superseded five-replication pilot.

The full launcher should be used for a clean study because it guarantees that
all nested `p` subsets in an arm come from the same `p_max = 2000` loading
matrix and that the DGP is identical across matched methods.

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

For Viroli Gibbs, the production launcher permits four independent task chunks
and four workers inside each fit:

```text
TASK_WORKERS_GIBBS=4
GIBBS_INTERNAL_WORKERS_SERIAL=4
GIBBS_INTERNAL_WORKERS_PARALLEL=4
```

These are the recorded production settings, not a guarantee that 16 cores are
busy continuously: the Gibbs blocks have different parallel fractions and R's
fork scheduling adds overhead. Runtime comparisons must therefore report the
worker policy alongside wall time. Override these variables to match a
different machine; do not change them silently within a reported run.

### Regenerate Simulation Figures

The full launcher plots all three arms automatically. For one completed arm,
the equivalent direct command is:

```sh
RUN_LABEL=fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_separated \
RESULTS_FILE=results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_separated/comparison_results.csv \
PLOT_DIR=results/selected_plots/sample_size/three_mixture_settings/separated \
TABLE_DIR=results/selected_tables/sample_size/three_mixture_settings \
METHOD_FILTER=independent_marginal_mixture,viroli_laplace_gibbs \
G_COMPONENT_FILTER=2,3 \
P_FILTER=500,1000,1500,2000 \
OUTPUT_TAG=separated_product_map_vs_gibbs \
Rscript scripts/sample_size/plot_fixed_ifeval_grouped_boxplot_panels.R
```

The same command with the other run label and output directory regenerates the
asymmetric and overlap figures. Product MAP occupies all four `p` columns;
Gibbs appears only at `p = 500, 1000`, by design.

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

- three-arm recovery plots:
  `results/selected_plots/sample_size/three_mixture_settings/`;
- three-arm cell means and integrity checks:
  `results/selected_tables/sample_size/three_mixture_settings/`;
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
