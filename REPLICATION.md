# Replication Instructions

This guide explains how to replicate the two analysis tracks in this repository:

1. the IFEval empirical analysis;
2. the sample-size simulation comparing the proposed factorial factor-mixture method with a joint-mixture Gibbs sampler.

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

This runs the singular-value shelf diagnostic, held-out-likelihood CV, selected
mixture refit, Lambda heatmaps, marginal mixture plots, loading summaries, and
factor-score visualizations.

### Singular-Value Shelf Diagnostic

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval/reproduced_rank_diagnostics \
H_MAX=10 \
Rscript scripts/ifeval/plot_ifeval_rank_diagnostics.R
```

Primary outputs:

- `ifeval_singular_value_shelf.csv`;
- `ifeval_singular_value_shelf.png`;
- `ifeval_intercept_only_alpha_for_rank_diagnostics.csv`.

This diagnostic is descriptive.  Model selection is done by held-out predictive
log likelihood.

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

The main simulation compares:

- `independent_marginal_mixture`: the proposed method with augmented-probit SVD pretraining, independent marginal mixture rotation, item intercepts, and MAP refinement;
- `joint_mixture_factor_gibbs`: a correctly specified joint-mixture factor model with `G^H` latent profiles and diagonal within-profile covariance, fitted by full Gibbs sampling with probit augmentation.

The main grid is:

- `H in {3, 4}`;
- `G in {2, 3}`;
- `n in {100, 500, 1000, 2000}`;
- `p = 500`;
- `25` Monte Carlo repetitions;
- loading designs `balanced_moderate_few_positive_cross` and `balanced_moderate_dense_signed_cross`;
- item intercept mode `ifeval_like`;
- unequal mixture variances;
- proposed-method pretraining max iterations `10`;
- proposed-method refinement max iterations `10`;
- Gibbs iterations `2000`, burn-in `1000`, thin `1`.

These settings are recorded in:

```sh
configs/sample_size_intercepts_centered.env
```

The unbalanced-block extension uses the same `H`, `G`, `n`, `p`, repetitions,
loading designs, intercepts, and MAP/Gibbs settings, but changes
`BLOCK_SIZE_MODE` to `moderate_ifeval_like`.  Its launcher is:

```sh
Rscript scripts/sample_size/run_ifeval_blocksize_crossloading_sample_size_MAP_intercepts.R
```

By default, outputs go to:

```text
results/full/moderate_ifeval_blocksize_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered
```

### Run The Full Simulation

Run:

```sh
Rscript scripts/sample_size/run_moderate_crossloading_sample_size_MAP_intercepts.R
```

By default, outputs go to:

```text
results/full/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered
```

The launcher has `RESUME_EXISTING=TRUE`, so rerunning the command resumes from existing checkpoint files rather than starting over.

Primary outputs:

- `comparison_results_checkpoint.csv`;
- `comparison_results.csv`;
- `comparison_summary.csv`;
- `comparison_convergence.csv`;
- per-repetition parameter-recovery CSV files;
- Gibbs history CSV/PNG files;
- checkpoint/final line plots.

### Run A Smoke Test

For a fast check of the code path:

```sh
cd "/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"

OUT_DIR=results/full/smoke_sample_size \
H_VALUES=3 \
G_VALUES=2 \
NP_GRID=n100p500:100:500 \
LOADING_DESIGNS=balanced_moderate_few_positive_cross \
REP_VALUES=1 \
PRETRAIN_AUG_ITER=3 \
REFINE_ITER=3 \
MFA_ITER=50 \
MFA_BURN=25 \
RESUME_EXISTING=FALSE \
Rscript scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R
```

This is only a code-path check, not a scientific simulation.

### Refresh Checkpoint Plots

For the main full run:

```sh
OUT_DIR=results/full/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered \
Rscript scripts/sample_size/plot_moderate_crossloading_checkpoint.R
```

Or use the refresher loop:

```sh
OUT_DIR=results/full/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered \
bash scripts/sample_size/refresh_moderate_crossloading_plots.sh
```

The plotting script produces:

- flat all-parameter correlation lines;
- lambda correlation lines;
- joint-mixture variance correlation lines;
- timing lines;
- progress heatmap;
- all-parameter panel plots by `H, G, loading_design`.

### Regenerate Paper-Facing Sample-Size Figures

The selected paper-facing simulation plots are generated from the checkpoint
CSV and the exact DGP loading generator. Set `OUT_DIR` to the full simulation
output directory:

```sh
OUT_DIR=results/full/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered
```

To regenerate the representative DGP loading heatmaps for
`Loadings = "Sparse"` and `Loadings = "Cross"` at `H = 3, 4`:

```sh
OUT_DIR=$OUT_DIR \
Rscript scripts/sample_size/plot_dgp_loading_heatmaps.R
```

This writes:

- `dgp_lambda_heatmap_sparse_H3.png`;
- `dgp_lambda_heatmap_sparse_H4.png`;
- `dgp_lambda_heatmap_cross_H3.png`;
- `dgp_lambda_heatmap_cross_H4.png`;
- matching `dgp_lambda_matrix_*.csv` files containing the exact representative
  loading matrices shown in the heatmaps.

To regenerate the six-panel recovery figures for each loading/H/G setting:

```sh
OUT_DIR=$OUT_DIR \
Rscript scripts/sample_size/plot_sample_size_rmse_panels.R
```

These figures report RMSE for item intercepts, loadings, mixture means,
mixture variances, and mixture weights, plus the flattened factor-score
correlation. The mixture-weight RMSE is computed by vectorizing the aligned
joint-profile weights and taking the RMSE against the true profile weights.

To regenerate the runtime figures:

```sh
OUT_DIR=$OUT_DIR \
Rscript scripts/sample_size/plot_sample_size_timing_lines.R
```

These figures plot mean `log(seconds)` with `+/- 2` standard deviations across
Monte Carlo repetitions. The `seconds` field is total elapsed wall-clock time
for each method fit; for the Gibbs comparator it includes the configured 2000
Gibbs iterations and 1000 burn-in iterations.

To refresh the committed selected copies after regenerating full outputs:

```sh
cp "$OUT_DIR"/dgp_lambda_heatmap_*.png results/selected_plots/sample_size/
cp "$OUT_DIR"/checkpoint_parameter_recovery_panel_*.png results/selected_plots/sample_size/
cp "$OUT_DIR"/checkpoint_timing_log_seconds_*.png results/selected_plots/sample_size/
cp "$OUT_DIR"/dgp_lambda_matrix_*.csv results/selected_tables/sample_size/
cp "$OUT_DIR"/checkpoint_parameter_recovery_summary.csv results/selected_tables/sample_size/
cp "$OUT_DIR"/checkpoint_timing_log_seconds_summary.csv results/selected_tables/sample_size/
```

For the unbalanced-block simulation, use the unbalanced output directory above
and copy the files with an `unbalanced_blocks_` prefix when preparing selected
paper-facing snapshots.

### Interpret Simulation Metrics

The main correlation metrics are computed after aligning estimated factors/loadings to the data-generating factors/loadings by best signed permutation.

Important fields:

- `factor_corr`: correlation between aligned estimated and true factor scores.
- `lambda_corr`: correlation between aligned estimated and true loading entries.
- `alpha_corr`: correlation between estimated and true item intercepts.
- `mixture_mu_corr`: correlation between estimated and true marginal mixture means.
- `mixture_var_corr`: correlation between estimated and true marginal mixture variances.
- `flat_parameter_corr`: correlation after flattening all comparable parameter blocks into one vector.
- `seconds`: wall-clock runtime for the method in that repetition.

For `G = 3`, variance recovery can be sensitive to saturated binary item blocks, so inspect both `mixture_var_corr` and variance RMSE/parameter recovery tables.

## Committed Selected Results

Selected outputs from prior runs are committed for immediate inspection:

- sample-size plots: `results/selected_plots/sample_size`;
- sample-size tables: `results/selected_tables/sample_size`;
- IFEval plots: `results/selected_plots/ifeval`;
- IFEval tables: `results/selected_tables/ifeval`;
- IFEval writeups: `writeup/ifeval_analysis_writeup.pdf` and
  `writeup/ifeval_componentwise_G3313.pdf`.

The latest static audit notes are in `CODE_AUDIT.md`. They record which R files
were parsed, which input files were checked, which PDF was rendered, and which
caveats remain before claiming a full fresh end-to-end rerun.

These selected outputs are snapshots. Full reruns should write to `results/full/`.

## Git Notes

Generated `.rds` files, logs, and full output directories are ignored by `.gitignore`. Commit source code, replication instructions, selected summary tables, selected plots, and final writeups.
