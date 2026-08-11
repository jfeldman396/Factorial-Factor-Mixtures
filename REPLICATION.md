# Replication Instructions

This guide explains how to replicate the two analysis tracks in this repository:

1. the IFEval empirical analysis;
2. the sample-size simulation comparing the proposed factorial factor-mixture method with a joint-mixture Gibbs sampler.

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
- `writeup/`: IFEval LaTeX source and rendered PDF.

Full regenerated outputs should be written under `results/full/`, which is ignored by git.

## IFEval Analysis

The cleaned input files are:

- `data/ifeval/openeval_ifeval_only_binary_matrix.csv`;
- `data/ifeval/openeval_item_metadata.csv`;
- `data/ifeval/openeval_model_metadata.csv`.

The main selected model is the independent-mixture binary probit factor model with `H = 3`, `G = 3`, item intercepts `alpha_j`, and sparse loadings selected by the lambda penalty tuning step.

### Refit The Selected Mixture Model

Run:

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
ITEM_METADATA_PATH=data/ifeval/openeval_item_metadata.csv \
OUT_DIR=results/full/ifeval/reproduced_openeval_ifeval_H3_G3_interpretation \
H_FIXED=3 \
G_FIXED=3 \
WORKERS=8 \
REFINEMENT_LAMBDA_L1_PENALTY=10 \
Rscript scripts/ifeval/fit_interpret_ifeval_H3_G3.R
```

Primary outputs:

- `openeval_item_intercepts_loadings_metadata.csv`;
- `openeval_model_factor_scores_profiles.csv`;
- `openeval_factor_interpretation_summary.csv`;
- `openeval_top_loading_item_examples.csv`;
- `openeval_factor_mixture_groups.csv`;
- loading and factor-score plots.

### Tune The Sparse Loading Penalty

Run:

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval/reproduced_openeval_ifeval_lambda_sparsity_tuning \
H_FIXED=3 \
G_FIXED=3 \
WORKERS=8 \
Rscript scripts/ifeval/tune_ifeval_lambda_sparsity.R
```

Primary outputs:

- `ifeval_lambda_sparsity_tuning_summary.csv`;
- `ifeval_selected_lambda_sparsity_penalty.csv`;
- `ifeval_lambda_sparsity_tuning_path.png`;
- fitted `.rds` files for each lambda penalty.

The committed selected analysis used lambda penalty `10`.

### Fit Ordinary Binary Probit Comparator

Run:

```sh
MMLU_MATRIX=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval/reproduced_openeval_ifeval_ordinary_probit_H3_visualization \
H=3 \
WORKERS=8 \
AUG_ITER=4 \
REFINE_ITER=10 \
LAMBDA_L1=2 \
Rscript scripts/ifeval/fit_visualize_ordinary_probit_factors.R
```

Primary outputs:

- `ordinary_probit_lambda.csv`;
- `ordinary_probit_factor_scores.csv`;
- `ordinary_probit_lambda_heatmap.png`;
- ordinary-factor score plots.

### Compare Ordinary And Mixture Factors

After refitting both the selected mixture model and ordinary comparator, run:

```sh
ORDINARY_DIR=results/full/ifeval/reproduced_openeval_ifeval_ordinary_probit_H3_visualization \
MIXTURE_DIR=results/full/ifeval/reproduced_openeval_ifeval_H3_G3_interpretation \
OUT_DIR=results/full/ifeval/reproduced_openeval_ordinary_vs_mixture_H3 \
Rscript scripts/ifeval/summarize_openeval_ordinary_probit_factors.R
```

Then regenerate the side-by-side loading plots:

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
MIXTURE_LOADINGS_PATH=results/full/ifeval/reproduced_openeval_ifeval_H3_G3_interpretation/openeval_item_intercepts_loadings_metadata.csv \
ORDINARY_LAMBDA_PATH=results/full/ifeval/reproduced_openeval_ifeval_ordinary_probit_H3_visualization/ordinary_probit_lambda.csv \
ORDINARY_MIXTURE_COR_PATH=results/full/ifeval/reproduced_openeval_ordinary_vs_mixture_H3/ordinary_mixture_factor_correlation_matrix.csv \
OUT_DIR=results/full/ifeval/reproduced_openeval_ordinary_vs_mixture_H3 \
Rscript scripts/ifeval/plot_openeval_side_by_side_lambdas.R
```

### Summarize Loadings And Cross-Loadings

Run:

```sh
LOADINGS_PATH=results/full/ifeval/reproduced_openeval_ifeval_H3_G3_interpretation/openeval_item_intercepts_loadings_metadata.csv \
OUT_DIR=results/full/ifeval/loadings_crossloadings \
Rscript scripts/ifeval/summarize_ifeval_loadings_crossloadings.R
```

Primary outputs:

- `ifeval_sparse_learned_lambda_all_items_ordered.csv`;
- `ifeval_sparse_lambda_matrix_threshold_0p5.csv`;
- `ifeval_G3_factor_only_exact_items.csv`;
- `ifeval_G3_cross_loading_items.csv`;
- heatmaps for sparse and cross-loading items.

### Regenerate 3D IFEval Visualizations

Base R plots:

```sh
FIT_DIR=results/full/ifeval/reproduced_openeval_ifeval_H3_G3_interpretation \
OUT_DIR=results/full/ifeval/ifeval_3d_factor_visualizations \
Rscript scripts/ifeval/plot_ifeval_learned_factors_baseR.R
```

Optional Python/Plotly plots:

```sh
FIT_DIR=results/full/ifeval/reproduced_openeval_ifeval_H3_G3_interpretation \
OUT_DIR=results/full/ifeval/ifeval_3d_factor_visualizations \
python3 scripts/ifeval/plot_ifeval_learned_factors.py
```

### Reproduce Rank-Selection Comparison Plots

The repository includes lightweight CV summary inputs for the mixture rank-selection fits under:

- `results/selected_tables/ifeval/cv_inputs/mixture_G2_H_summary.csv`;
- `results/selected_tables/ifeval/cv_inputs/mixture_G2_fold_scores.csv`;
- `results/selected_tables/ifeval/cv_inputs/mixture_G3_H_summary.csv`;
- `results/selected_tables/ifeval/cv_inputs/mixture_G3_fold_scores.csv`.

To rerun the ordinary-probit CV from the raw IFEval matrix:

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval/openeval_ifeval_cv_H1_10_ordinary \
H_GRID=1:10 \
K_FOLDS=3 \
WORKERS=8 \
AUG_ITER=4 \
REFINE_ITER=5 \
LAMBDA_L1=2 \
Rscript scripts/ifeval/openeval_ordinary_probit_cv_H_selection.R
```

Then regenerate the rank-selection comparison figures:

```sh
G2_H_SUMMARY=results/selected_tables/ifeval/cv_inputs/mixture_G2_H_summary.csv \
G2_FOLD_SCORES=results/selected_tables/ifeval/cv_inputs/mixture_G2_fold_scores.csv \
G3_H_SUMMARY=results/selected_tables/ifeval/cv_inputs/mixture_G3_H_summary.csv \
G3_FOLD_SCORES=results/selected_tables/ifeval/cv_inputs/mixture_G3_fold_scores.csv \
ORDINARY_H_SUMMARY=results/full/ifeval/openeval_ifeval_cv_H1_10_ordinary/ordinary_probit_H_summary.csv \
ORDINARY_FOLD_SCORES=results/full/ifeval/openeval_ifeval_cv_H1_10_ordinary/ordinary_probit_factor_fold_scores.csv \
OUT_DIR=results/full/ifeval/reproduced_openeval_ifeval_rank_selection_comparison \
Rscript scripts/ifeval/summarize_openeval_rank_selection_cv.R
```

The heavy per-fold mixture-CV `.rds` files are not committed. The committed summaries are sufficient to regenerate the paper-facing rank-selection comparison tables and plots.

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
- IFEval writeup: `writeup/ifeval_analysis_writeup.pdf`.

These selected outputs are snapshots. Full reruns should write to `results/full/`.

## Git Notes

Generated `.rds` files, logs, and full output directories are ignored by `.gitignore`. Commit source code, replication instructions, selected summary tables, selected plots, and final writeups.
